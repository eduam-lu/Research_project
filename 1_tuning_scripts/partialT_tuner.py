# -*- coding: utf-8 -*-
"""
partialT_tuner.py input_folder input_csv_file

This script takes a set of selected structures (in an input csv file) from an input folder and generates a dataframe to assess what's the best partial T for a project

1. For each file, runs RF diffusion 3 times 
2. For each RF diffusion prediction, computes 20 seqs with ProteinMPNN and selects the one with the lowest score
3. Predicts the structure with Alphafold  3
4. RMSD calculation with biopython
5. Appends all the info to a final dataframe with the following columnns
    [file_ID] [original_sequence] [MPNN_seq] [sequence_length] [DSSP_original] [H%_original] [E%_original] [DSSP_final] [H%_final] [E%_final] [category] [partial_T] [pLDDT] [RMSD]
    
"""
#%% import modules
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess
import re
import json
from statistics import mean, stdev
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
#%%Functions
def parse_MPNN_file (file):
    # Given an MPNN formatted file it parses it sequence by sequence, stores each sequence in a list and each score in another list
    final_scores = []
    final_sequences = []
    with open(file, "r") as faa_file :
        for line in faa_file:
            if line.startswith(">"):
                match = re.search(r"score=([\d.]+)", line)
                if match:
                    score = float(match.group(1))
                    final_scores.append(score)
            else:
                final_sequences.append(line.strip())
    return final_sequences,final_scores

def best_seq_extractor (file, seq_info, score_info):
    # Given a file name,  a list with all the sequences and a list with all the scores. It returns a list with the sequence with the lowest score and the score
    index_of_lowest = score_info.index(min(score_info))
    return [str(file), seq_info[index_of_lowest], score_info[index_of_lowest]]

def parse_MPNN_folder(folder):
    # Given a folder with MPNN output files, it takes each file and
    # a. stores all the scores and seqs in a list
    # b. find the tuple with the minimum score and appends it to the global list
    # c. returns a list with all the best sequences and a list with the scores
    input_folder = Path(folder)
    best_sequences = list()
    for file in input_folder.iterdir():
        if file.is_file():
            extracted_seqs,extracted_scores = parse_MPNN_file(file)
            best_seq_info = best_seq_extractor(file,extracted_seqs,extracted_scores)
            best_sequences.append(best_seq_info)            
            
    return best_sequences

def extract_plddt_from_cif(cif_path):
    plddt_scores = []
    with open(cif_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                fields = line.strip().split()  # or split() if it's space-delimited
                try:
                    score = float(fields[14])  # 15th field (index 14)
                    plddt_scores.append(score)
                except (ValueError, IndexError):
                    continue
    return plddt_scores

def pLDDT_summary_generator(scores):
    summary_stats = [round(mean(scores), 2),round(stdev(scores), 2),round(min(scores), 2), round(max(scores), 2)]
    return summary_stats

def load_structure(path, struct_id):
        ext = Path(path).suffix.lower()
        if ext == '.pdb':
            parser = PDBParser(QUIET=True)
        elif ext == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            raise ValueError("Unsupported file format. Use .pdb or .cif")
        return parser.get_structure(struct_id, path)

def calculate_rmsd(structure_path1, structure_path2, chain_id='A'):
    """
    Calculates RMSD between two structures (.pdb or .cif) using CA atoms from a specified chain.

    Parameters:
    - structure_path1 (str): Path to the first structure file (.pdb or .cif).
    - structure_path2 (str): Path to the second structure file.
    - chain_id (str): Chain identifier (default: 'A').

    Returns:
    - float: RMSD value.
    """

    structure1 = load_structure(structure_path1, "struct1")
    structure2 = load_structure(structure_path2, "struct2")

    atoms1 = [atom for atom in structure1[0][chain_id].get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom for atom in structure2[0][chain_id].get_atoms() if atom.get_id() == 'CA']

    if len(atoms1) != len(atoms2):
        raise ValueError("Structures do not have the same number of CA atoms.")

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(structure2.get_atoms())

    return sup.rms

def seq_comparison(seq1, seq2):
    """
    Compare two sequences and count the number of positions where they differ.

    Parameters:
    - seq1 (str): First sequence
    - seq2 (str): Second sequence

    Returns:
    - int: Number of differing positions
    """
    # Compare up to the length of the shortest sequence
    min_len = min(len(seq1), len(seq2))
    differences = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    
    # Optionally count differences due to length mismatch
    differences += abs(len(seq1) - len(seq2))

    return differences

#%% Input check

parser = argparse.ArgumentParser(description="Generate summary dataframes for a project structures")
parser.add_argument('--folder', help="Folder that contains all the input structures", type=str, required= True)
parser.add_argument('--sample', help="File that contains the selected structures from the folder", type=str, required= True)
args = parser.parse_args()
folder = Path(args.folder)
sample = Path(args.sample)

#%% Main execution
#### IMPORTANT #####
# This script must be run in madoka and in my home directory

# Make sure we are in Madoka

# Generate necessary output folders
Path("partial_T_tuning/rf_outputs").mkdir(parents=True, exist_ok=True)
Path("partial_T_tuning/MPNN_outputs").mkdir(parents=True, exist_ok=True)
Path("partial_T_tuning/AF3_outputs").mkdir(parents=True, exist_ok=True)
Path("partial_T_tuning/AF3_outputs/json_inputs").mkdir(parents=True, exist_ok=True)
# Store the samples in a dataframe
sample_df = pd.read_csv(sample)

# Parse the file IDs and run 3 different RF diffusions for each
file_list = list(sample_df['file_ID'])
length_list = list(sample_df['length'])
partial_Ts = [10,20,30]

for T in partial_Ts:
    output_path = f"./partial_T_tuning/rf_outputs/{T}_"
    for file,length in zip(file_list, length_list):
        path_to_pdb = folder / file
        formated_length = f"[{length},{length}]"
        rf_command = [
            "export MKL_THREADING_LAYER=GNU && conda", "run", "-n", "SE3nv", "python", 
            "~/RFdiffusion/scripts/run_inference.py",
            f"inference.output_prefix={output_path}{file}",
            f"inference.input_pdb={path_to_pdb}",
            f"contigmap.contigs=\'[{length}-{length}]\'",
            "inference.num_designs=1",
            f"diffuser.partial_T={T}"
        ]
        
        # Join manually to preserve quotes
        full_command = " ".join(rf_command)
        #subprocess.run(full_command, shell=True)
# Move all the pdbs to a single folder
# Define source and destination folders
src_folder = Path("./partial_T_tuning/rf_outputs")
dst_folder = Path("./partial_T_tuning/rf_outputs/rf_pdbs")

# Make sure the destination folder exists
dst_folder.mkdir(parents=True, exist_ok=True)

# Move each .pdb file
for pdb_file in src_folder.glob("*.pdb"):
    new_location = dst_folder / pdb_file.name
    pdb_file.rename(new_location)
# Parse the output folder for RFdiffussion, generate 20 seqs for each file

rf_pdb_folder = "./partial_T_tuning/rf_outputs/rf_pdbs"
MPNN_output_path = "partial_T_tuning/MPNN_outputs"

json_command = (
    "python ~/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path "
    + rf_pdb_folder
    + " --output_path "
    + rf_pdb_folder
    + "/test.jsonl"
)

subprocess.run(json_command, shell=True)

MPNN_command = (
    "python ~/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target=20 --batch_size=10 --out_folder=" #ask about batch_size
    + MPNN_output_path
    + " --use_soluble_model --sampling_temp 0.1 --save_scores --save_probs MPPN_probs --jsonl_path="
    + rf_pdb_folder
    + "/test.jsonl" # ask about sampling 
)

subprocess.run(MPNN_command, shell=True)

# Parse the MPNN output to select the sequence with the least score and append the sequence and the score to a list
sequences = parse_MPNN_folder("partial_T_tuning/MPNN_outputs/seqs") # Contains a list of lists with the format [[file, sequence, score], ...]

### Generate alphafold predictions with the sequence list
# Set up variables

output_dir = "partial_T_tuning/AF3_outputs"
env_path = "/home/ingemar/anaconda3/envs/alphafold3"
json_path = Path("partial_T_tuning/AF3_outputs/json_inputs")

# Generate json inputs for each sequence
    
for sequence_info in sequences:
    file_name = Path(sequence_info[0]).stem
    seq = sequence_info[1]
    json_data = {
        "name": file_name,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }
    
#    with open(json_path / f"{file_name}.json", "w") as f:
 #       json.dump(json_data, f, indent=2)

# Parse the json files and generate predictions


for file in json_path.iterdir():
    af_command = (
    f"conda run -p {env_path} python /mnt/data/alphafold3/run_alphafold.py "
    f"--json_path={file} "
    f"--output_dir={output_dir} "
    f"--db_dir=/mnt/data/alphafold3/alphafold_databases/ "
    f"--model_dir=/mnt/data/alphafold3/models/"
    )
    #af_process = subprocess.run(af_command, shell = True, text=True)

# Parse pLDDT and compute RMSD with the original structure
# I wil parse the sequences, updating the list from [file, seq, score] to [file, seq, score, pLDDT_mean, "stdev_pLDDT", "min_pLDDT","max_pLDDT", "RMSD"]
sequences_plus_pLDDT = []
for sequence_info in sequences:
    cif_name = str(Path(sequence_info[0]).stem)[:-6]
    pdb_name = cif_name[3:]
    print(pdb_name)
    AF3_cif_path = f"partial_T_tuning/AF3_outputs/{cif_name.lower()}.pdb_0/seed-1_sample-0/model.cif"
    original_pdb_path = f"input_pdbs/{pdb_name}.pdb"

    pLDDT_scores = extract_plddt_from_cif(AF3_cif_path)
    pLDDT_summary = pLDDT_summary_generator(pLDDT_scores)
    RMSD = calculate_rmsd(AF3_cif_path, original_pdb_path, "A")

    updated_information = list(sequence_info)
    
    # Flatten the list and convert np.float64 to regular float
    updated_information.extend([float(x) for x in pLDDT_summary])
    updated_information.append(float(RMSD))
    
    sequences_plus_pLDDT.append(updated_information)

print(sequences_plus_pLDDT[1])

# Update the original dataframe with all the new information
# First I split the partial T used from the name of the file and clean it
final_sequences =[]
for seq_info in sequences_plus_pLDDT :
    cleaned_info = list(seq_info)
    cleaned_info[0] = cleaned_info[0][35:-5]
    partial_T_used = cleaned_info[0][:2]
    cleaned_info[0] = cleaned_info[0][3:]
    cleaned_info.append(partial_T_used)
    final_sequences.append(cleaned_info)
print(final_sequences[1])
print(final_sequences[5])

# Final sequences format ["file_name", "MPNN_seq", "pLDDT_mean", "pLDDT_mean", "pLDDT_mean", "pLDDT_mean", "RMSD", "partial_T"]
# Now i use the first element of each final sequence to search for the proper row in the original df, and concat those rows together
final_df = pd.DataFrame(columns=list(sample_df.columns)+[ "MPNN_seq", "MPNN_score", "pLDDT_mean", "pLDDT_stdev", "pLDDT_min", "pLDDT_max", "RMSD", "partial_T"])
for sequence_info in final_sequences:
    entry_id = sequence_info[0]
    new_info = sequence_info[1:]
    row = sample_df[sample_df['file_ID'] == entry_id].values.flatten().tolist()
    row.extend(new_info)
    small_df = pd.DataFrame([row], columns=final_df.columns)
    final_df = pd.concat([final_df,small_df])

# Now I wanna compute the number of changes in the original sequence and the MPNN sequence
for index, row in final_df.iterrows():
    n_changes = seq_comparison(str(row['sequence']), str(row['MPNN_seq']))
    final_df.loc[index, 'changes'] = n_changes
# Then I reorder the columns properly
desired_order = ["file_ID","category", "partial_T", "pLDDT_mean", "pLDDT_stdev", "pLDDT_min", "pLDDT_max", "RMSD","sequence","MPNN_seq","length", "changes", "MPNN_score","DSSP","H%","E%"]
final_df = final_df[desired_order]
final_df.to_csv('partial_T_sampling.csv')
