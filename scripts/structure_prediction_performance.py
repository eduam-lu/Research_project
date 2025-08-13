# -*- coding: utf-8 -*-
"""
structure_prediction_performance.py input_sequence_file

This script takes an input sequence and predicts its structure with 3 different methods
 - Alpha Fold 3 without an MSA template
 - ESM 
 - Chai 1

Apart from the structures, it outputs summaries of the quality of the prediction and the time it took to predict them
    
"""

#%% Modules
import pandas as pd
import subprocess
import time
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import json
from statistics import mean


#%% Functions
def run_ESM2 (row, output_path, original_structure_path):
    sequence = row['MPNN_seq']
    file_name = row['file_ID']
    partial_T = row['partial_T']
    # Launch sequence to the server
    esm_command = f"curl -X POST --data \"{sequence}\" https://api.esmatlas.com/foldSequence/v1/pdb/ > {output_path}/{partial_T}_{file_name.lower()}"
    start = time.time()
    esm_process = subprocess.run(esm_command, shell = True, text=True)
    end = time.time()
    run_time = end - start
    # Get pLDDT
    new_structure_path = f"{output_path}/{partial_T}_{file_name.lower()}"
    pLDDT = extract_plddt(new_structure_path,10)
    # Get RMSD
    RMSD = calculate_rmsd(original_structure_path, new_structure_path)
    return run_time, 100*pLDDT, RMSD

def run_AF3_wo_MSA(row, output_path, original_structure_path):
    # Create temporary JSON file
    json_generator(row, output_path)
    # Run AF3
    af_command = (
    f"conda run -p /home/ingemar/anaconda3/envs/alphafold3 python /mnt/data/alphafold3/run_alphafold.py "
    f"--json_path={output_path}/temp.json "
    f"--output_dir={output_path} "
    f"--db_dir=/mnt/data/alphafold3/alphafold_databases/ "
    f"--model_dir=/mnt/data/alphafold3/models/ "
    f"--num_diffusion_samples=1"
    )
    start = time.time()
    af_process = subprocess.run(af_command, shell = True, text=True)
    end = time.time()
    run_time = end - start
    # Get pLDDT
    new_structure_path = f"{output_path}/{row['partial_T']}_{file_name.lower()}/seed-1_sample-0/model.cif"
    pLDDT = extract_plddt(new_structure_path,14)
    # Get RMSD
    RMSD = calculate_rmsd(original_structure_path, new_structure_path)
    return run_time, pLDDT, RMSD

def run_chai(row, output_path, original_structure_path):
    sequence = row['MPNN_seq']
    file_name = row['file_ID']
    partial_T = row['partial_T']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{partial_T}{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
    # Run chai
    chai_command =f" conda run -n chai2 chai-lab fold {fasta_path} {output_path}/{partial_T}_{file_name.lower()}" 
    start = time.time()
    chai_process = subprocess.run(chai_command, shell=True)
    end= time.time()
    run_time = end - start

    # Get pLDDT # I need to know the format for the pLDDT
    new_structure_path = f"{output_path}/{partial_T}_{file_name.lower()}/pred.model_idx_0.cif"
    pLDDT = extract_plddt(new_structure_path,17)

    # Get RMSD
    RMSD = calculate_rmsd(original_structure_path, new_structure_path)
    return run_time, pLDDT, RMSD

### AUXILIARY FUNCTIONS

def json_generator(row, json_path):
    file_name = row['file_ID']
    seq = row['MPNN_seq']
    partial_T = row['partial_T']

    json_data = {
        "name": f"{partial_T}_{file_name}",
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq,
                    "unpairedMsa": f">dummy\n{seq}\n",  # <--- String directamente, no dict
                    "pairedMsa": "",                    # <--- también string vacío
                    "templates": []
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }

    with open(f"{json_path}/temp.json", "w") as f:
        json.dump(json_data, f, indent=2)



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

def extract_plddt(cif_path,position = 14):
    plddt_scores = []
    with open(cif_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                fields = line.strip().split()  # or split() if it's space-delimited
                try:
                    score = float(fields[position])  # 15th field (index 14)
                    plddt_scores.append(score)
                except (ValueError, IndexError):
                    continue
    return mean(plddt_scores)



#%% MAIN

# Extract IDs and sequences from a data frame
original_df = pd.read_csv("partial_T_sampling.csv")

sequence_df = original_df[['file_ID', 'MPNN_seq', 'partial_T']]
# Feed each sequence iteratively to each method
column_comp = ["file_ID","partial_T","sequence","af_time","esm_time","chai_time","af_pLDDT","esm_pLDDT","chai_pLDDT","af_RMSD","esm_RMSD","chai_RMSD"]
comparison_df = pd.DataFrame(columns=column_comp) #Initialise storage dataframe

#Generate output folders
Path("prediction_methods_comparison/af_outputs").mkdir(parents=True, exist_ok=True)
Path("prediction_methods_comparison/chai_outputs").mkdir(parents=True, exist_ok=True)
Path("prediction_methods_comparison/esm_outputs").mkdir(parents=True, exist_ok=True)


for index, row in sequence_df.iterrows():
    sequence = row['MPNN_seq']
    file_name = row['file_ID']
    partial_T = row['partial_T']
    original_structure_path = f'input_pdbs/{file_name}'
    # Run Alphafold3 wo MSA
    run_time_af, pLDDT_af, RMSD_af = run_AF3_wo_MSA(row, "prediction_methods_comparison/af_outputs",original_structure_path)

    # Run ESM
    run_time_esm, pLDDT_esm, RMSD_esm = run_ESM2(row, "prediction_methods_comparison/esm_outputs",original_structure_path)

    # Run chai
    run_time_chai, pLDDT_chai, RMSD_chai = run_chai(row, "prediction_methods_comparison/chai_outputs",original_structure_path)

    # New row as a Series
    elements = [file_name,partial_T,sequence,run_time_af,run_time_esm,run_time_chai,pLDDT_af,pLDDT_esm,pLDDT_chai,RMSD_af,RMSD_esm,RMSD_chai]
    new_row = pd.Series(elements, index=column_comp)
    # Add the Series as a new row
    comparison_df.loc[len(comparison_df)] = new_row

comparison_df.to_csv('prediction_methods_comparison.csv')




