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

#%%Functions
def parse_MPNN_file ():
    # Given an MPNN formatted file it parses it sequence by sequence, storing in a list a tuple of sequence and score
    
    return
def parse_MPNN_folder():
    # Given a folder with MPNN output files, it takes each file and
    # a. stores all the scores and seqs in a list
    # b. find the tuple with the minimum score and appends it to the global list
    # c. returns a list with all the best sequences and a list with the scores
    return

def extract_partialT_df ():
    # Given a data frame it takes the first two digits of the columnn 'file_ID' and stores them in a new columm 'partial_T'
    # Then it changes the file name removing the first 3 digits
    # For example: 10_edu.pdb would have a file ID of edu.pdb and a partial_T of 10
    
    return
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
#        subprocess.run(full_command, shell=True)
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

#subprocess.run(json_command, shell=True)

MPNN_command = (
    "conda run -n SE3nv python ~/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target=20 --batch_size=10 --out_folder=" #ask about batch_size
    + MPNN_output_path
    + " --use_soluble_model --sampling_temp 0.1 --jsonl_path="
    + rf_pdb_folder
    + "/test.jsonl" # ask about sampling 
)

subprocess.run(MPNN_command, shell=True)

# Parse the MPNN output to select the sequence with the least score and append the sequence and the score to a list

# Generate alphafold predictions with the sequence list

# RMSD with the original structure file

# Update the original dataframe with all the new information
