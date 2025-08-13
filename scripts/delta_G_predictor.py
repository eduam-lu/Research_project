#%% Import and install libraries, load models
import os,time,subprocess,re,sys
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def format_pytorch_version(version):
  return version.split('+')[0]

def format_cuda_version(version):
  return 'cu' + version.replace('.', '')

TORCH_version = torch.__version__
print(TORCH_version)
TORCH = format_pytorch_version(TORCH_version)
CUDA_version = torch.version.cuda
print(CUDA_version)
CUDA = format_cuda_version(CUDA_version)

IF_model_name = "~/scripts/esm_if1_gvp4_t16_142M_UR50.pt"

import numpy as np
import argparse
import json
from pathlib import Path
import esm
import biotite.structure
# alias the new filter function so the old name still works
biotite.structure.filter_backbone = biotite.structure.filter_peptide_backbone
from esm.inverse_folding.util import load_structure, extract_coords_from_structure,CoordBatchConverter
from esm.inverse_folding.multichain_util import extract_coords_from_complex,_concatenate_coords,load_complex_coords

# Tell PyTorch it’s OK to unpickle argparse.Namespace
torch.serialization.add_safe_globals([argparse.Namespace])
print("importing the model")
# Get the folder where your script lives:
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the absolute path to your .pt file
model_path = os.path.join(script_dir, "esm_if1_gvp4_t16_142M_UR50.pt")
model, alphabet = esm.pretrained.load_model_and_alphabet(model_path)
model.eval().cuda().requires_grad_(False)

print("--> Installations succeeded")

#%% FUNCTIONS

def run_model(coords,sequence,model,cmplx=False,chain_target='A'):

    device = next(model.parameters()).device

    batch_converter = CoordBatchConverter(alphabet)
    batch = [(coords, None, sequence)]
    coords, confidence, strs, tokens, padding_mask = batch_converter(
        batch, device=device)

    prev_output_tokens = tokens[:, :-1].to(device)
    target = tokens[:, 1:]
    target_padding_mask = (target == alphabet.padding_idx)

    logits, _ = model.forward(coords, padding_mask, confidence, prev_output_tokens)

    logits_swapped=torch.swapaxes(logits,1,2)
    token_probs = torch.softmax(logits_swapped, dim=-1)

    return token_probs

def score_variants(sequence,token_probs,alphabet):

    aa_list=[]
    wt_scores=[]
    skip_pos=0

    alphabetAA_L_D={'-':0,'_' :0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}
    alphabetAA_D_L={v: k for k, v in alphabetAA_L_D.items()}

    for i,n in enumerate(sequence):
      aa_list.append(n+str(i+1))
      score_pos=[]
      for j in range(1,21):
          score_pos.append(masked_absolute(alphabetAA_D_L[j],i, token_probs, alphabet))
          if n == alphabetAA_D_L[j]:
            WT_score_pos=score_pos[-1]

      wt_scores.append(WT_score_pos)

    return aa_list, wt_scores

def masked_absolute(mut, idx, token_probs, alphabet):

    mt_encoded = alphabet.get_idx(mut)

    score = token_probs[0,idx, mt_encoded]
    return score.item()

#%% Input check
parser = argparse.ArgumentParser(description="Generate summary dataframes for a project structures")
parser.add_argument('--input_folder', help="Folder that contains all the input structures", type=str, required= True)
parser.add_argument('--output_folder', help="Folder that contains all the output structures", type=str, required= True)
args = parser.parse_args()
input_path = Path(args.input_folder)
output_path = Path(args.output_folder)

#%% MAIN

a=0.10413378327743603 ## fitting param from the manuscript to convert IF score scale to kcal/mol
b=0.6162549378400894 ## fitting param from the manuscript to convert IF score scale to kcal/mol

#Initialise dictionary for storing the delta g
delta_g_dict={}

for file in list(input_path.glob("*.pdb")) + list(input_path.glob("*.cif")):  # loop over all .pdb and .cif files
    #print(f"Processing {file.name}...")

    # Load structure
    structure = load_structure(str(file), 'A')  # load current file
    coords_structure, sequence_structure = extract_coords_from_structure(structure)

    # Compute probabilities
    prob_tokens = run_model(coords_structure, sequence_structure, model, chain_target='A')
    aa_list, wt_scores = score_variants(sequence_structure, prob_tokens, alphabet)

    # Compute delta G in kcal/mol
    dg_IF = np.nansum(wt_scores)
    dg_kcalmol = a * dg_IF + b

    delta_g_dict[file.name] = dg_kcalmol  # use file name as the key

# Save delta_g_dict as JSON
with open(f"{output_path}/delta_g_dict.json", "w") as f:
    json.dump(delta_g_dict, f)