"""
This script serves as a control to run shapedesign without an RF diffusion step
"""
#IMPORT MODULES ################################################################################################################################
from pathlib import Path
import argparse
import sys,os
import logging
import subprocess
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, Superimposer, NeighborSearch, PDBIO
import pandas as pd
import pymol
import tempfile
import json
from pymol import cmd
pymol.finish_launching(['pymol', '-cq'])
### FUNCTIONS ##################################################################################################################################
def help():
    return

def rename_files_to_lowercase(folder_path="."):
    for filename in os.listdir(folder_path):
        src = os.path.join(folder_path, filename)
        if os.path.isfile(src):
            new_name = filename.lower()
            dst = os.path.join(folder_path, new_name)
            if src != dst:
                os.rename(src, dst)

def run_shapedesign(file_path, output_path):
    pure_name = str(Path(file_path).name).split("_relaxed")[0]
    sym_def_path = f"symdefs/{pure_name}.symm"
    shapedesign_command = f"conda run -n zernike_pyrosetta python proteinshapedesign/run_proteinshapedesign.py --struct {file_path} --symdef {sym_def_path} --out_dir {output_path}  --popsize 4 --generations 1 --cycles 0"
    shapedesign_process = subprocess.run(shapedesign_command,shell=True)
    return

def minimise_structure(folder,output):
    minimise_command = f" conda run -p /home/ingemar/miniforge3/envs/pyrosetta2 python scripts/structure_minimisation.py --input_structure {folder} --output_folder {output}"
    minimise_process = subprocess.run(minimise_command, shell=True, text=True)
    return

def capsid_metrics(df, name, design_path):
    # Load design csv from shapedesign folder
    metric = pd.read_csv(design_path)

    # Add name column
    metric['name'] = name

    # If df is empty, initialize it with the same columns as metric
    if df.empty:
        df = pd.DataFrame(columns=metric.columns)
    else:
        df.columns = metric.columns  # ensures compatibility if needed

    # Concat both dataframes
    final_df = pd.concat([df, metric], ignore_index=True)
    return final_df

def filter_capsids(df,input_dir,output_dir):
    # Given a df w/ capsid metrics, returns a 1-row df with the best capsid
    # According to: sasa, ddg and sc

    # Columns to normalize
    columns = [
        "ddg.5fold", "ddg.3fold", "ddg.2_1fold", "ddg.2_2fold",
        "sasa.5fold", "sasa.3fold", "sasa.2_1fold", "sasa.2_2fold",
        "sasa_polar.5fold", "sasa_polar.3fold", "sasa_polar.2_1fold", "sasa_polar.2_2fold",
        "sasa_hydrophobic.5fold", "sasa_hydrophobic.3fold", "sasa_hydrophobic.2_1fold", "sasa_hydrophobic.2_2fold",
        "sc.5fold", "sc.3fold", "sc.2_1fold", "sc.2_2fold"
    ]

    norm_columns = []

    # Normalize desired columns
    for col in columns:
        norm_col = f"{col}_norm"
        norm_columns.append(norm_col)
        df[norm_col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())

    # Generate score as the mean of normalized columns
    df["score"] = df[norm_columns].sum(axis=1) / len(norm_columns)

    # Select row with the best score
    best_row = df.loc[df["score"].idxmax()]
    best_row = best_row.to_frame().T
    
    #Move file to the definitive folder
    best_file = best_row['name'].values[0]
    best_file = f"{best_file.split(".remove")[0]}_design_full.cif"
    source_path = f"{input_dir}/structures/{best_file}"
    dest_path = f"{output_dir}/{best_file}"
    os.rename(source_path, dest_path) 

    return best_row  # Return as a 1-row DataFrame

def extraction(folder,output_path):
    # Generate subfolders
    monomer_path = f"{output_path}/monomers"
    dimer_path = f"{output_path}/dimers"
    trimer_path = f"{output_path}/trimers"
    pentamer_path = f"{output_path}/pentamers"
    Path(monomer_path).mkdir(exist_ok=True, parents =True)
    Path(dimer_path).mkdir(exist_ok=True, parents =True)
    Path(trimer_path).mkdir(exist_ok=True, parents =True)
    Path(pentamer_path).mkdir(exist_ok=True, parents =True)
    # Determine chains to extract
    monomer_chains =['A']
    dimer_chains =['A','B']
    trimer_chains =['A','J','K']
    pentamer_chains =['A','B','C','D','E']


    # Parse the input
    for file in Path(folder).iterdir():
        # Monomers
        extract_chains(file, monomer_path, monomer_chains,"monomer")
        # Dimers
        extract_chains(file, dimer_path, dimer_chains,"dimer")
        # Trimers
        extract_chains(file, trimer_path, trimer_chains, "trimer")
        # Pentamers
        extract_chains(file, pentamer_path, pentamer_chains, "pentamer")
    return

def extract_chains(pdb_path, output_path, chains_to_extract, label):
    """
    Use PyMOL to extract specified chains from a PDB and save to a new file.

    Parameters:
    - pdb_path: Path object
    - output_path: Path to folder
    - chains_to_extract: list of str (e.g., ['A', 'B'])
    - label: str, a label to append to the output filename (e.g., 'dimer')
    """
    base_name = pdb_path.stem
    output_file = Path(output_path) / f"{base_name}_{label}.pdb"

    cmd.load(str(pdb_path), 'structure')
    selection_str = " or ".join([f"chain {c}" for c in chains_to_extract])
    cmd.select("extracted", selection_str)
    cmd.save(str(output_file), "extracted")
    cmd.delete("all")

def extract_seq_pdb(file_path):
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"{file_path} does not exist.")
    
    if file_path.suffix.lower() == '.pdb':
        parser = PDBParser(QUIET=True)
    elif file_path.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
        #raise ValueError("Unsupported file type. Only .pdb and .cif are supported.")
    
    structure = parser.get_structure("structure", str(file_path))
    ppb = PPBuilder()
    
    # Extract first chain with a sequence
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()
        return str(seq), len(seq)

    raise ValueError("No polypeptide chain with a recognizable sequence was found.")

def predict_solubility(seq_ID, sequence):
    # Create temp fasta file with the sequence
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_file:
        temp_file.write(f">{seq_ID}\n")
        temp_file.write(sequence)
        temp_path = temp_file.name  # Save the path for later

    try:
        # Predict solubility
        command = f"protein-sol-sequence-prediction-software/multiple_prediction_wrapper_export.sh {temp_path}"
        subprocess.run(command, shell=True, check=True)

        # Retrieve estimation from seq_prediction.txt
        with open("seq_prediction.txt", "r") as file:
            for line in file:
                if line.startswith("SEQUENCE PREDICTIONS"):
                    solubility_estimate = line.strip().split(",")[3]
                    return solubility_estimate

        raise ValueError("Solubility estimate not found in output.")

    finally:
        # Clean up temp file
        if os.path.exists(temp_path):
            os.remove(temp_path)

# Kyte-Doolittle hydrophobicity scale
hydrophobicity_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def estimate_hydrophobicity(sequence):
    estimate = 0
    for aa in sequence:
        estimate += hydrophobicity_scale[aa]
    estimate = estimate/len(sequence)
    return estimate

def calculate_stability(folder_path,df, dict_output):
    # Better to process in folders, because that way, the parameters are loaded only once
    stability_command = f" conda run -n esm_stability python scripts/delta_G_predictor.py --input_folder {folder_path} --output_folder {dict_output}"
    stability_process = subprocess.run(stability_command, shell=True, text=True)

    # Load the generated dictionary and save it as a dataframe
    with open(f"{dict_output}/delta_g_dict.json", 'r') as file:
        data = json.load(file)
    delta_G_df = pd.DataFrame(list(data.items()), columns=["file_ID", "stab_shape"])
    # Join both dataframes so that the stability column is added
    print(df['file_ID'])
    print (delta_G_df['file_ID'])
    merged_df = pd.merge(df, delta_G_df, on='file_ID', how='inner')

    return merged_df

def convert_all_cif_to_pdb(input_folder, output_folder):
    command = f"mamba run -n benchmark python scripts/pdb_converter.py {input_folder} {output_folder}"
    process = subprocess.run(command, shell=True)

def move_files(source_folder, destination_folder, ending):
    # Ensure destination folder exists
    os.makedirs(destination_folder, exist_ok=True)

    # Iterate through files in source folder
    for filename in os.listdir(source_folder):
        if filename.endswith(ending):
            source_path = os.path.join(source_folder, filename)
            destination_path = os.path.join(destination_folder, filename)
            os.rename(source_path, destination_path)
### INPUT CHECK ################################################################################################################################
parser = argparse.ArgumentParser(
    description="This script runs the pipeline shapedesign with an additional RF_diffusion step"
)
parser.add_argument('--folder', help="Folder that contains all the input structures", type=str)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# If --detailed-help is provided, print custom help and exit
if args.detailed_help:
    help()
    sys.exit()

# If --folder was not provided, show error and exit
if not args.folder:
    parser.error("--folder is required unless --detailed-help is used")

# Convert to Path object
input_folder = Path(args.folder)

# If input folder is not an existing directory, show error and exit
if not input_folder.exists():
    parser.error(f"The folder '{input_folder}' does not exist.")
if not input_folder.is_dir():
    parser.error(f"The path '{input_folder}' is not a directory.")


### GENERATE FOLDER STRUCTURE ###################################################################################################################
# Generate the global folder, making sure that it doesn't exist already
global_output_name=Path("control_shapedesign_output_29-06")
#count=1
#while global_output_name.exists():
#    global_output_name = Path(f"{global_output_name}_{count}")
#    count += 1
global_output_name.mkdir(parents=True, exist_ok=True)

### SET UP LOG FILE #############################################################################################################################
logging.basicConfig(
    filename=f"{str(global_output_name)}/improved_shapedesign.log",       # Log file path
    filemode='a',                   # Append mode
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)

def save_to_log(message):
    logging.info(message)

### MAIN EXECUTION################################################################################################################################

#save_to_log("Starting execution")
# 0. Format the names of all the input_files, so that they are all lowercase ###########################################
#save_to_log("Processing input files")
#rename_files_to_lowercase(input_folder)

### 1. Relax structures
minimised_folder = f"{str(global_output_name)}/minimised_structures"
Path(minimised_folder).mkdir(parents=True, exist_ok=True)

minimise_structure(input_folder,minimised_folder)
### 2. Shapedesign iterations

#save_to_log("Starting shapedesign iterations")
save_to_log("Starting shapedesign iterations")
shapedesign_output = f"{str(global_output_name)}/shapedesign"
Path(shapedesign_output).mkdir(parents=True, exist_ok=True)
final_capsids = f"{str(global_output_name)}/final_capsids"
Path(final_capsids).mkdir(parents=True, exist_ok=True)
indiv_capsid_metrics = pd.DataFrame()
global_capsid_metrics = pd.DataFrame()

for file in Path(minimised_folder).iterdir():
    # Assuming file is a Path object pointing to the original input structure
    new_file_name = file

    for i in range(5):
        old_file_name = new_file_name
        new_file_name = file.parent / f"{file.stem}_capsid_{i}.remove"

        # Rename old_file_name to new_file_name (if needed)
        Path(old_file_name).rename(new_file_name)

        # Run shape design
        run_shapedesign(new_file_name, shapedesign_output)

        # Append metrics for each generated capsid
        indiv_capsid_metrics = capsid_metrics(
            indiv_capsid_metrics,
            new_file_name.name,
            f"{str(shapedesign_output)}/design.csv"
        )

    # Filter the best capsid from those generated for this file
    selected_capsid = filter_capsids(indiv_capsid_metrics, shapedesign_output,final_capsids)

    # Add selected capsid to the global dataframe
    global_capsid_metrics = pd.concat([global_capsid_metrics, selected_capsid], ignore_index=True)

    # Reinitialize individual metrics DataFrame
    indiv_capsid_metrics = pd.DataFrame()



global_capsid_metrics.to_csv(f"{str(global_output_name)}/capsid_metrics.csv")
#move_files(f"{str(global_output_name)}/shapedesign/structures", "96_capsids_wo","full.cif")

### Extract monomers from the selected capsids
extracted_folder = f"{global_output_name}/extracted_pdbs"
Path(extracted_folder).mkdir(exist_ok=True,parents=True)
extraction(final_capsids,extracted_folder)
### Compute metrics
monomer_folder = f"{global_output_name}/extracted_pdbs/monomers"

monomer_df = pd.DataFrame(columns=["file_ID","seq","length","sol_shape","hydro_shape"])

for file in Path(monomer_folder).iterdir():
   file_name = file.name
   seq,length = extract_seq_pdb(file)
   sol = predict_solubility(file_name,seq)
   hydro = estimate_hydrophobicity(seq)
   row = pd.Series([file_name,seq,length,sol,hydro],index=["file_ID","seq","length","sol_shape","hydro_shape"])
   monomer_df.loc[len(monomer_df)] = row

monomer_df=calculate_stability(input_folder,monomer_df,f"{global_output_name}")
monomer_df.to_csv(f"{global_output_name}/monomer_metrics.csv")

