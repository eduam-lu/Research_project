#
"""

""" 
import subprocess
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, Superimposer, NeighborSearch, PDBIO
from pathlib import Path
from statistics import mean, stdev
import random
import re,os
import json
import tempfile

### MODEL FUNCTIONS ###################################################################################################################

def run_rf_partial(row, original_path, output_path):
    T = 20
    inferences = 10
    length = row["length"]  # Fix: access dict key with quotes
    file_id = row["file_ID"]

    # Proper formatting: Use f-strings, escape nested quotes
    rf_command = (
        f'export MKL_THREADING_LAYER=GNU && conda run -n SE3nv python '
        f'~/RFdiffusion/scripts/run_inference.py '
        f'inference.output_prefix="{output_path}/{file_id}_partialdiff" '
        f'inference.input_pdb="{original_path}{file_id}" '
        f'contigmap.contigs="[{length}-{length}]" '
        f'inference.num_designs={inferences} '
        f'diffuser.partial_T={T}'
    )
    subprocess.run(rf_command, shell=True)
    return

def run_rf_scaffold(row, original_path, output_path):
    inferences = 10
    file_id = row["file_ID"]
    length = row["length"]

    scaffold_dir = f"{output_path}/scaffolds/{file_id}"
    orig_pdb = f"{original_path}/{file_id}"

    # Run the helper script
    scaffold_command = (
        f'export MKL_THREADING_LAYER=GNU && conda run -n SE3nv python '
        f'~/RFdiffusion/helper_scripts/make_secstruc_adj.py '
        f'--input_pdb "{orig_pdb}" --out_dir "{scaffold_dir}"'
    )
    subprocess.run(scaffold_command, shell=True)

    # Now run the inference script
    rf_command = (
        f'export MKL_THREADING_LAYER=GNU && conda run -n SE3nv python '
        f'~/RFdiffusion/scripts/run_inference.py '
        f'inference.output_prefix="{output_path}/{file_id}_conditioning" '
        f'inference.num_designs={inferences} '
        f'contigmap.contigs="[{length}-{length}]" '
        #f'scaffoldguided.scaffoldguided=True '
        f'scaffoldguided.target_pdb=False '
        f'scaffoldguided.scaffold_dir="{scaffold_dir}" '
        f'denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0'
    )
    subprocess.run(rf_command, shell=True)
    return

def run_MPNN(pdb_folder, output_path):

    json_command = (
        "python ~/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path "
        + pdb_folder
        + " --output_path "
        + pdb_folder
        + "/test.jsonl"
    )

    subprocess.run(json_command, shell=True)

    MPNN_command = (
        "export MKL_THREADING_LAYER=GNU && conda run -n SE3nv python ~/ProteinMPNN/protein_mpnn_run.py --num_seq_per_target=40 --batch_size=10 --out_folder=" #ask about batch_size
        + output_path
        + " --use_soluble_model --sampling_temp 0.1  --jsonl_path="
        + pdb_folder
        + "/test.jsonl" # ask about sampling 
    )

    subprocess.run(MPNN_command, shell=True)
    return

def run_ESM (row, output_path):
    sequence = row['seq']
    file_name = row['file_ID']

    # Launch sequence to the server
    esm_command = f"curl -X POST --data \"{sequence}\" https://api.esmatlas.com/foldSequence/v1/pdb/ > {output_path}/{file_name.lower()}.pdb"
    esm_process = subprocess.run(esm_command, shell = True, text=True)

    return 

def run_chai(row, output_path):
    sequence = row['seq']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
    # Run chai
    chai_command =f" conda run -n chai2 chai-lab fold {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return 


def run_shapedesign(file_path, output_path):
    pure_name = str(Path(file_path).name).split(".pdb")[0]
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
### 1D PARAMETERS FUNCTIONS ######################################################################################################################################################################

# As of now i only include one 1D filter, which gets rids of sequences that have more than 4 times the same aa in a row
# This gets rid of alanine and glutamate runs

# Might be interesting to think about more 1D filters in the future

def oned_sequence_filter(seq):
    """Returns True if sequence has bad features, False if it passes."""
    REPEAT_REGEX = r"(.)\1{3,}"  # 4+ repeated amino acids
    return re.search(REPEAT_REGEX, seq) is not None

def oned_dataframe_filter(df):
    """Filters out sequences that fail 1D criteria (e.g., repetitive residues)."""
    # Apply filter and keep only passing sequences
    df["filter"] = df["seq"].apply(oned_sequence_filter)
    checked_df = df[~df["filter"]].copy()  # Keep only sequences where filter is False
    return checked_df[["file_ID", "seq"]]

### 3D PARAMETERS FUNCTIONS ########################################################################################################################################################################
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
    if plddt_scores:
        return mean(plddt_scores),stdev(plddt_scores)
    else:
        return 0,0


def load_structure(path, struct_id):
    ext = Path(path).suffix.lower()
    if ext == '.pdb':
        parser = PDBParser(QUIET=True)
    elif ext == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
        #raise ValueError("Unsupported file format. Use .pdb or .cif")
    return parser.get_structure(struct_id, path)

def detect_clashes(file_path, struct_id="structure", clash_distance=2.0, bond_distance=1.2):
    """
    Detects steric clashes in a protein structure.

    Args:
        file_path (str): Path to structure file (.pdb or .cif).
        struct_id (str): Structure ID.
        clash_distance (float): Distance threshold to define clashes (Å).
        bond_distance (float): Distance threshold to consider atoms bonded (Å).

    Returns:
        tuple:
            - int: number of clashes
            - float: clashes per atom
            - list of tuples: clashing atom pairs
    """
    structure = load_structure(file_path, struct_id)
    atoms = [atom for atom in structure.get_atoms() if not atom.get_name().startswith("H")]
    ns = NeighborSearch(atoms)
    clashes = set()
    clash_pairs = []

    for atom in atoms:
        close_atoms = ns.search(atom.coord, clash_distance)
        for neighbor in close_atoms:
            if atom == neighbor:
                continue

            # Avoid double-counting pairs
            pair = frozenset((atom, neighbor))
            if pair in clashes:
                continue

            # Exclude likely bonded atoms (e.g., within 1.2 Å)
            distance = atom - neighbor
            if distance < bond_distance:
                continue

            if distance < clash_distance:
                clashes.add(pair)
                clash_pairs.append((atom, neighbor))

    num_clashes = len(clashes)
    num_atoms = len(atoms)
    clashes_per_atom = num_clashes / num_atoms if num_atoms > 0 else 0.0

    return num_clashes, clashes_per_atom

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
        return 1000

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(structure2.get_atoms())

    return sup.rms

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
    delta_G_df = pd.DataFrame(list(data.items()), columns=["file_ID", "delta_G"])
    # Join both dataframes so that the stability column is added
    print(df['file_ID'])
    print (delta_G_df['file_ID'])
    merged_df = pd.merge(df, delta_G_df, on='file_ID', how='inner')

    return merged_df

def threed_params_df(folder, original_path, original_df, pLDDT_loc):
    """
    Given a folder with structures, returns the following df: "file_ID", "sequence", "stability", "pLDDT mean", "pLDDT std dev", "RMSD", "clashes"
    """
    # Initialise df
    elements = ["file_ID", "seq", "pLDDT_mean", "pLDDT_stdev", "RMSD", "clashes","clashes_per_atom","orig_delta_G"]
    threed_df = pd.DataFrame(columns = elements)

    # Parse folder
    folder = Path(folder)
    for file in folder.iterdir():
        file_path = str(file)
        if not file_path.endswith(".json"):
            original_name = str(file.name).split(".")[0] + ".pdb" # Remove all the locators added during the pipeline
            original_file = f"{original_path}/{original_name}"
            # pLDDT
            pLDDT_mean,pLDDT_stdev= extract_plddt(file_path,pLDDT_loc)
            # RMSD
            RMSD = calculate_rmsd(file_path,original_file)
            # Clashes
            num_clashes,clash_atom = detect_clashes(file_path)
            # Sequence
            seq = extract_seq_pdb(file_path)[0]
            # Original stab
            orig_stab = original_df[original_df['file_ID'] == original_name]["delta_G"].iloc[0]
            # Generate row
            row = pd.Series([file.name, seq, pLDDT_mean, pLDDT_stdev, RMSD, num_clashes,clash_atom, orig_stab], index=list(threed_df.columns))
            threed_df.loc[len(threed_df)] = row
    
    # Add new stability
    print(threed_df.head(5))
    threed_df = calculate_stability(folder,threed_df,folder)
    print(threed_df.head(5))
    return threed_df



def threed_filter_df(df):
    MIN_PLDDT = 0.8
    MIN_RMSD = 1
    MAX_RMSD = 8
    MAX_CLASHES = 1.01
    # Drop rows with low pLDDT
    df = df[df["pLDDT_mean"] >= MIN_PLDDT]

    # Drop rows with RMSD out of range
    df = df[(df["RMSD"] >= MIN_RMSD) & (df["RMSD"] <= MAX_RMSD)]

    # Drop rows with too many clashes
    df = df[df["clashes_per_atom"] <= MAX_CLASHES]

    return df

def final_3d_filter(df,original_df,global_output):
    MIN_PLDDT = 0.8
    MIN_RMSD = 1
    MAX_RMSD = 8
    MAX_CLASHES = 1.01
    # Drop rows with low pLDDT
    df = df[df["pLDDT_mean"] >= MIN_PLDDT]

    # Drop rows with RMSD out of range
    df = df[(df["RMSD"] >= MIN_RMSD) & (df["RMSD"] <= MAX_RMSD)]

    # Drop rows with too many clashes
    df = df[df["clashes_per_atom"] <= MAX_CLASHES]

    #Parses the file_ID column in order to identify how many different types there are
    file_list = []
    # Initialize new columns with None
    for index, row in df.iterrows():
        file_name = row['file_ID'].split(".")[0] + ".pdb"

        # Optional: track processed files
        if file_name not in file_list:
            file_list.append(file_name)

    # Initialise final df
    elements = list(df.columns) + ["delta_delta_Tm"]
    final_df = pd.DataFrame(columns=elements)
    # Parses all the possible origin inputs and stores them in a sepparate df. For that df
    #   - The delta delta G is calculated
    #   - Structures are ordered by pLDDT and delta delta G
    #   - A dataframe with the top 5 candidates is appended to the final df
    for file in file_list:
        file_df = df[df['file_ID'].str.startswith(file)]
        # Calculate delta_Tm
        # original_Tm = original_df[original_df['file_ID'] == file]['stability'].iloc[0]
        #original_Tm = 1
        #file_df['delta_Tm'] = file_df["stability"] - original_Tm
        # Normalize values between 0 and 1
        file_df["delta_G_norm"] = (file_df["delta_G"] - file_df["delta_G"].min()) / (file_df["delta_G"].max() - file_df["delta_G"].min())
        file_df["pLDDT_norm"] = (file_df["pLDDT_mean"] - file_df["pLDDT_mean"].min()) / (file_df["pLDDT_mean"].max() - file_df["pLDDT_mean"].min())
        # Combine both with equal weight (or adjust weights as needed)
        file_df["score"] = 0.5 * file_df["delta_G_norm"] + 0.5 * file_df["pLDDT_norm"]
        # Sort by the combined score
        sorted_df = file_df.sort_values(by="score", ascending=False).reset_index(drop=True)
        sorted_df = sorted_df.head(5)
        sorted_df = sorted_df.sort_values(by="delta_G", ascending=False).reset_index(drop=True)
        final_df = pd.concat([final_df,sorted_df.head(1)])
    
    final_df['seq'] = None
    final_df['sol_RF'] = None
    final_df['hydro_RF'] = None

    final_df = final_df.reset_index(drop=True)

    for i in range(len(final_df)):
        row = final_df.iloc[i]
        file_name = row['file_ID'].split(".")[0] + ".pdb"
        file_path = f"{global_output}/Chai_pdb_files/{str(row['file_ID'])}"

        seq, length = extract_seq_pdb(file_path)
        final_df.loc[i, 'seq'] = seq

        sol = predict_solubility(file_name, seq)
        final_df.loc[i, 'sol_RF'] = sol

        hydro = estimate_hydrophobicity(seq)
        final_df.loc[i, 'hydro_RF'] = hydro
    return file_list,final_df
### AUXILIARY FUNCTIONS ###########################################################################################################################################
def help():
    pass

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

def process_MPNN_folder(folder):

    folder = Path(folder)

    sequence_df = pd.DataFrame(columns=["file_ID","seq","score"]) 
    for file in folder.iterdir():
        file_df = process_MPNN_file(file)
        sequence_df = pd.concat([sequence_df, file_df], ignore_index = True)

    return sequence_df

def process_MPNN_file(file):
    # Process the file
    scores = []
    sequences = []
    with open(file, "r") as faa_file :
        for line in faa_file:
            if line.startswith(">"):
                match = re.search(r"score=([\d.]+)", line)
                if match:
                    score = float(match.group(1))
                    scores.append(score)
            else:
                sequences.append(line.strip())
    
    #Store as a dataframe
    sequence_df = pd.DataFrame(columns=["file_ID","seq","score"]) 
    file_name = str(file.name)
    file_name = file_name[:-3]
    count = 1
    for score_seq, seq in zip(scores,sequences):
        entry_name = f"{file_name}_seq_{count}"
        row = pd.Series([entry_name, seq, score_seq], index=["file_ID","seq","score"])
        sequence_df.loc[len(sequence_df)] = row
        count += 1
    
    # Sort by score descending and keep top 20
    top_20_df = sequence_df.sort_values(by="score", ascending=True).head(20).reset_index(drop=True)

    return top_20_df

def format_chai_folder(source_dir,dest_dir):
    os.makedirs(dest_dir, exist_ok=True)

    for subfolder in os.listdir(source_dir):
        subfolder_path = os.path.join(source_dir, subfolder)

        if os.path.isdir(subfolder_path):
            original_file = os.path.join(subfolder_path, "pred.model_idx_0.cif")

            if os.path.isfile(original_file):
                new_filename = f"{subfolder}.cif"
                destination_file = os.path.join(dest_dir, new_filename)

                # Copy the file content manually
                with open(original_file, 'rb') as src_file:
                    content = src_file.read()

                with open(destination_file, 'wb') as dst_file:
                    dst_file.write(content)

                print(f"Copied and renamed: {original_file} -> {destination_file}")
            else:
                print(f"File not found: {original_file}")

def folder_from_dataframe(df, original_folder, output_folder):
    """
    Given a DataFrame, this function reads the 'file_ID' column,
    finds each corresponding file in original_folder, and copies it
    to output_folder using manual read/write.
    """
    os.makedirs(output_folder, exist_ok=True)
    
    for index, row in df.iterrows():
        file_id = row['file_ID']
        original_path = os.path.join(original_folder, file_id)
        destination_path = os.path.join(output_folder, file_id)
        
        if os.path.isfile(original_path):
            # Manually copy file using read/write
            with open(original_path, 'rb') as src:
                data = src.read()
            with open(destination_path, 'wb') as dst:
                dst.write(data)
        else:
            print(f"File not found: {original_path}")

def remove_small_files(folder_path, min_size=1024):
    """
    Remove all files smaller than `min_size` bytes in the specified folder.

    Args:
        folder_path (str): Path to the folder.
        min_size (int): Minimum file size in bytes (default: 1024 = 1 KB).
    """
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        return

    removed_files = 0
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            size = os.path.getsize(file_path)
            if size < min_size:
                os.remove(file_path)
                removed_files += 1

    return removed_files

def rename_files_to_lowercase(folder_path="."):
    for filename in os.listdir(folder_path):
        src = os.path.join(folder_path, filename)
        if os.path.isfile(src):
            new_name = filename.lower()
            dst = os.path.join(folder_path, new_name)
            if src != dst:
                os.rename(src, dst)

def convert_all_cif_to_pdb(input_folder, output_folder):
    command = f"conda run -n benchmark python scripts/pdb_converter.py {input_folder} {output_folder}"
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