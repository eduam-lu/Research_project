#
"""
This script contains the functions needed to run the script benchmark.py
Function index:
A.
"""
### IMPORT MODULES ##################################################################################################
import subprocess
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, Superimposer, NeighborSearch
from pathlib import Path
from statistics import mean, stdev
import random
import re,os
import pymol
from pymol import cmd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import tmscoring
from tmscoring import TMscoring
import json
import tempfile
pymol.finish_launching(['pymol', '-cq'])

### EXTRACTION FUNCTIONS ############################################################################################

def extraction(folder,output_path):
    """
    Use PyMOL to extract specified chains from a PDB and save to a new file.

    Parameters:
    - pdb_path: Path object
    - output_path: Path to folder
    - chains_to_extract: list of str (e.g., ['A', 'B'])
    - label: str, a label to append to the output filename (e.g., 'dimer')
    """
    
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
    
    return

### PREDICTION FUNCTIONS ############################################################################################

def run_chai_mer(file, output_path, n_mer:str):
    
    mer_dict = {"monomer":1,
                "dimer": 2,
                "trimer": 3,
                "pentamer":5}
    
    sequence = extract_seq_pdb(file)[0]
    
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        for i in range(mer_dict[n_mer]):
            f.write(f">protein|name={file}-{i}\n")      # Header line
            f.write(f"{sequence}\n")         # The actual sequence

    # Run chai
    chai_command = f"export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True && conda run -n chai2 chai-lab fold {fasta_path} {output_path}/{file.name}_{n_mer}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return

def run_chai_folder(folder,output_path):

    # Mer list

    mer_list = ["monomer", "dimer", "trimer", "pentamer"]

    # Parse folder

    for file in Path(folder).iterdir():
        mer_output= f"{output_path}/{file.name}"
        Path(mer_output).mkdir(exist_ok=True, parents=True)
        for mer in mer_list:
            run_chai_mer(file,mer_output,mer)

    return


### METRICS FUNCTIONS ###############################################################################################
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



def align_structures_2(pdb1_path, pdb2_path, output_folder, image_name="alignment.png", image_width=1200, image_height=900):
    """
    Aligns two multi-chain PDB structures using PyMOL, saves a PNG of the alignment,
    and returns the RMSD value.

    Parameters:
    - pdb1_path (str or Path): path to the first PDB file
    - pdb2_path (str or Path): path to the second PDB file
    - output_folder (str or Path): where to save the alignment image
    - image_name (str): output filename for the PNG image
    - image_width (int): image width in pixels
    - image_height (int): image height in pixels

    Returns:
    - float: RMSD value from alignment
    """

    cmd.reinitialize()  # clean any previous session


    cmd.load(str(pdb1_path), "mobile")
    cmd.load(str(pdb2_path), "target")

    # Align structures (across all chains)
    rmsd = cmd.align("mobile", "target")[0]  # returns tuple, first is RMSD

    # Style and color
    cmd.color("cyan", "mobile")
    cmd.color("magenta", "target")
    cmd.show("cartoon", "all")
    cmd.bg_color("white")

    # Save image
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    image_file = output_path / image_name
    cmd.png(str(image_file), width=image_width, height=image_height, ray=1)

    # Use tmscoring if aligning monomers, MMalign if aligning multimers
    if "monomer" in str(pdb1_path.name):
        tmscore = calculate_tmscore_monomer(pdb1_path,pdb2_path)
    else:
        tmscore = calculate_tmscore_multimer(pdb1_path,pdb2_path)
    return rmsd, tmscore

def convert_cif_to_temp_pdb(cif_path):
    """
    Converts a CIF file to PDB format using PyMOL and returns the temporary PDB path.
    """
    import tempfile
    from pymol import cmd

    # Remove everything currently loaded to avoid naming conflicts
    cmd.delete("all")

    # Load CIF file
    cmd.load(str(cif_path))

    # Get the name of the object just loaded
    objects = cmd.get_object_list()
    if not objects:
        raise RuntimeError(f"Failed to load CIF file: {cif_path}")
    obj_name = objects[0]

    # Create temporary PDB file
    temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    cmd.save(temp_pdb.name, obj_name)

    return temp_pdb.name

def align_structures(pdb1_path, pdb2_path, output_folder, image_name="alignment.png", image_width=1200, image_height=900):
    """
    Aligns two multi-chain structures using PyMOL, saves a PNG of the alignment,
    and returns the RMSD and TM-score.

    Handles CIF input by converting it to PDB first.
    """
    # Convert to Path objects
    pdb1_path = Path(pdb1_path)
    pdb2_path = Path(pdb2_path)

    # Convert CIF files to PDB if needed
    temp_files = []  # Keep track of temp files for cleanup
    if pdb1_path.suffix == ".cif":
        pdb1_path = Path(convert_cif_to_temp_pdb(pdb1_path))
        temp_files.append(pdb1_path)
    if pdb2_path.suffix == ".cif":
        pdb2_path = Path(convert_cif_to_temp_pdb(pdb2_path))
        temp_files.append(pdb2_path)

    # Load structures into PyMOL
    cmd.reinitialize()
    cmd.load(str(pdb1_path), "mobile")
    cmd.load(str(pdb2_path), "target")

    # Align structures
    rmsd = cmd.align("mobile", "target")[0]

    # Style and image
    cmd.color("cyan", "mobile")
    cmd.color("magenta", "target")
    cmd.show("cartoon", "all")
    cmd.bg_color("white")

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    image_file = output_path / image_name
    cmd.png(str(image_file), width=image_width, height=image_height, ray=1)

    # Use appropriate TM-score method
    if "monomer" in str(pdb1_path.name):
        tmscore = calculate_tmscore_multimer(pdb1_path, pdb2_path)
    else:
        tmscore = calculate_tmscore_multimer(pdb1_path, pdb2_path)

    # Clean up temp files
    for f in temp_files:
        try:
            Path(f).unlink()
        except Exception as e:
            print(f"Warning: could not delete temp file {f}: {e}")

    return rmsd, tmscore

def convert_cif_to_pdb_pymol(cif_path, pdb_path):
    name = "temp_structure"
    cmd.load(cif_path, name)
    cmd.save(pdb_path, name)
    cmd.delete(name)

def calculate_tmscore_monomer(cif_file1, cif_file2):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb1 = os.path.join(tmpdir, "temp1.pdb")
        pdb2 = os.path.join(tmpdir, "temp2.pdb")

        # Convert CIF to PDB using PyMOL
        convert_cif_to_pdb_pymol(cif_file1, pdb1)
        convert_cif_to_pdb_pymol(cif_file2, pdb2)

        # Compute TM-score using tmscoring on the temporary PDB files
        alignment = TMscoring(pdb1, pdb2)
        score = alignment.tmscore(**alignment.get_current_values()) 

    # Temp files auto-removed when context exits
    return score


def calculate_tmscore_multimer_3(pdb1_path, pdb2_path, mmalign_path="MMalign"):
    """
    Run MM-align on two multimeric PDB files and return the global TM-score.

    Parameters:
    - pdb1_path (str): Path to the first PDB file (model A)
    - pdb2_path (str): Path to the second PDB file (model B)
    - mmalign_path (str): Path to the MM-align binary (default assumes it's in $PATH)

    Returns:
    - float: Global TM-score (A->B)
    """
    try:
        result = subprocess.run(
            [mmalign_path, pdb1_path, pdb2_path],
            capture_output=True,
            text=True,
            check=True
        )

        # Look for global TM-score
        match = re.search(r"TM-score=([\d\.]+) \(if normalized by length of chain A\)", result.stdout)
        if match:
            tm_score = float(match.group(1))
            return tm_score
        else:
            raise ValueError("Global TM-score not found in MM-align output.")

    except subprocess.CalledProcessError as e:
        print("MM-align failed:", e.stderr)
        return None
    
def calculate_tmscore_multimer2(pdb1_path, pdb2_path, mmalign_path="MMalign"):
    """
    Run MM-align on two multimeric PDB files and return the global TM-score.

    Parameters:
    - pdb1_path (str): Path to the first PDB file (model A)
    - pdb2_path (str): Path to the second PDB file (model B)
    - mmalign_path (str): Path to the MM-align binary (default assumes it's in $PATH)

    Returns:
    - float: Global TM-score (A->B), normalized by average length
    """
    try:
        result = subprocess.run(
            [mmalign_path, pdb1_path, pdb2_path, '-ter', '1'],
            capture_output=True,
            text=True,
            check=True
        )

        # Extract TM-score normalized by average length
        match = re.search(r"TM-score=([\d\.]+) \(if normalized by average length\)", result.stdout)
        if match:
            tm_score = float(match.group(1))
            return tm_score
        else:
            raise ValueError("Average-normalized TM-score not found in MM-align output.")

    except subprocess.CalledProcessError as e:
        print("MM-align failed:", e.stderr)
        return None

def calculate_tmscore_multimer(pdb1_path, pdb2_path, mmalign_path="MMalign"):
    try:
        result = subprocess.run(
            [mmalign_path, pdb1_path, pdb2_path, '-ter', '1'],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=True
        )

        print(result.stdout)  # Debug output

        matches = re.findall(r"TM-score=\s*([\d\.]+)", result.stdout)
        if matches:
            tm_score = float(matches[0])  # Or choose another index if needed
            return tm_score
        else:
            raise ValueError("TM-score not found in MM-align output.")

    except subprocess.CalledProcessError as e:
        print("MM-align failed:", e.stdout)
        return None

def clean_file_id(file_id):
    if "_relaxed_design_full_monomer" in file_id:
        return file_id.replace("_relaxed_design_full_monomer.pdb_monomer.pdb", "")
    else:
        return file_id.split(".pdb")[0]

def calculate_stability(folder_path,df, dict_output,col_name):
    # Better to process in folders, because that way, the parameters are loaded only once
    stability_command = f" conda run -n esm_stability python scripts/delta_G_predictor.py --input_folder {folder_path} --output_folder {dict_output}"
    stability_process = subprocess.run(stability_command, shell=True, text=True)

    # Load the generated dictionary and save it as a dataframe
    with open(f"{dict_output}/delta_g_dict.json", 'r') as file:
        data = json.load(file)
    delta_G_df = pd.DataFrame(list(data.items()), columns=["file_ID", col_name])
    delta_G_df["file_ID"] = delta_G_df["file_ID"].apply(clean_file_id)
    df["file_ID"] = df["file_ID"].apply(clean_file_id)
    print(df['file_ID'].iloc[0])
    print(delta_G_df['file_ID'].iloc[0])
    # Join both dataframes so that the stability column is added
    merged_df = pd.merge(df, delta_G_df, on='file_ID', how='inner')

    return merged_df

def calculate_stability_orig(folder_path, df, dict_output):
    # Run the external script using conda environment
    stability_command = f"conda run -n esm_stability python scripts/delta_G_predictor.py --input_folder {folder_path} --output_folder {dict_output}"
    subprocess.run(stability_command, shell=True, text=True)

    # Load the generated dictionary
    with open(f"{dict_output}/delta_g_dict.json", 'r') as file:
        data = json.load(file)

    # Prepare a list to store stability values
    stab_values = []

    # Extract the relevant stability values based on file_ID
    for idx, row in df.iterrows():
        name = row['file_ID'].split(".pdb")[0] + ".pdb"
        stab_values.append(data.get(name, None))  # Use .get() to handle missing keys gracefully

    # Add the new column to the DataFrame
    df['stab_orig'] = stab_values

    return df

def convert_cif_to_pdb(cif_path):
    tmp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    tmp_pdb.close()  # We will write to it via PyMOL
    object_name = os.path.basename(cif_path).replace('.', '_')
    cmd.load(cif_path, object_name)
    cmd.save(tmp_pdb.name, object_name)
    cmd.delete(object_name)
    return tmp_pdb.name

def parse_binding_energy_from_stdout(stdout):
    in_interaction_block = False
    for line in stdout.splitlines():
        if "interaction between" in line:
            in_interaction_block = True
        elif in_interaction_block and "Total" in line:
            try:
                return float(line.strip().split()[-1])
            except ValueError:
                pass
    raise ValueError("Could not find binding energy in FoldX output")

def interface_quality(pdb_file, output_dir):
    pdb_file = str(pdb_file)
    
    # Handle CIF conversion if needed
    if pdb_file.endswith(".cif"):
        print(f"Converting {pdb_file} to .pdb using active PyMOL session...")
        pdb_file = convert_cif_to_pdb(pdb_file)
        delete_temp = True
    else:
        delete_temp = False

    # Split path and prepare FoldX call
    path_parts = pdb_file.split("/")
    pdb_name = path_parts[-1]
    pdb_dir = "/".join(path_parts[:-1])

    cmd_args = [
        "FoldX",
        "--command=AnalyseComplex",
        f"--pdb={pdb_name}",
        "--complexWithDNA=false",
        f"--pdb-dir={pdb_dir}",
        f"--output-dir={output_dir}"
    ]
    print("Running FoldX:", " ".join(cmd_args))
    result = subprocess.run(cmd_args, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"FoldX failed:\n{result.stderr}")

    if delete_temp:
        os.remove(pdb_file)

    return parse_binding_energy_from_stdout(result.stdout)

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


### DATAFRAME FUNCTIONS #############################################################################################
def monomer_df(folder_w, folder_wo, global_output, orig_folder):
    df = pd.DataFrame()

    folder_w_path = Path(folder_w)
    folder_wo_path = Path(folder_wo)

    # Extraer claves de los archivos de folder_w
    files_w_map = {}
    for f in folder_w_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_w_map[key] = f.name
            
    # Extraer claves de los archivos de folder_wo
    files_wo_map = {}
    for f in folder_wo_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_wo_map[key] = f.name

    # Intersección de claves
    common_keys = set(files_w_map.keys()) & set(files_wo_map.keys())

    # Eliminar archivos en folder_w que no estén en la intersección
    for key, fname in files_w_map.items():
        if key not in common_keys:
            (folder_w_path / fname).unlink()

    # Eliminar archivos en folder_wo que no estén en la intersección
    for key, fname in files_wo_map.items():
        if key not in common_keys:
            (folder_wo_path / fname).unlink()
    
    print(f"📄 Number of files in {folder_w}: {len(os.listdir(folder_w))}")
    print(f"📄 Number of files in {folder_wo}: {len(os.listdir(folder_wo))}")
    # Parse folders, making sure that the files are ordered and match
    for file_w,file_wo in zip(Path(folder_w).iterdir(),Path(folder_wo).iterdir()):
        row = monomer_row(file_w,file_wo,global_output,orig_folder)
        df = pd.concat([df,row])
    # Make folders with pdbs
    Path(f"{folder_w}/pdbs").mkdir(exist_ok=True, parents=True)
    Path(f"{folder_wo}/pdbs").mkdir(exist_ok=True, parents=True)

    # Print number of files in the input folder before conversion
    print(f"📄 Number of files in {folder_w}: {len(os.listdir(folder_w))}")
    print(f"📄 Number of files in {folder_wo}: {len(os.listdir(folder_wo))}")

    # Run conversion
    convert_all_cif_to_pdb(folder_w, f"{folder_w}/pdbs")
    convert_all_cif_to_pdb(folder_wo, f"{folder_wo}/pdbs")
    # Stability w
    df = calculate_stability(f"{folder_w}/pdbs",df,global_output,'stab_w')

    # Stability wo
    df = calculate_stability(f"{folder_wo}/pdbs",df,global_output,'stab_wo')

    df.to_csv(f"{global_output}/just_in_case.csv")
    # Stability orig
    df = calculate_stability_orig(orig_folder,df,global_output)


    # Calculate delta delta G from original
    df['delta_delta_g_w'] = df['stab_w']- df['stab_orig']
    df['delta_delta_g_wo'] = df['stab_wo']- df['stab_orig']
    return df

def monomer_row(file_w, file_wo,global_output, orig_folder):
    # Find name
    file_ex = Path(file_w)
    name = str(file_ex.name).split("_monomer.pdb")[0]
    file_ex = Path(file_wo)
    name_wo = str(file_ex.name).split("_monomer.pdb")[0]
    # Find original monomer
    file_orig = f"{orig_folder}/{name}.pdb"
    file_orig = file_orig.split(".pdb")[0] + ".pdb"
    # Find extracted
    extracted_file_name_w = f"{global_output}/extracted_w/monomers/{name}_monomer.pdb"
    extracted_file_name_wo = f"{global_output}/extracted_wo/monomers/{name_wo}_monomer.pdb"
    extracted_file_name_wo =extracted_file_name_wo.split(".pdb")[0] + f".pdb_relaxed_capsid_{extracted_file_name_wo.split("_capsid_")[1][0]}_design_full_monomer.pdb"
    # Alignment folders 
    image_output_w = f"{global_output}/alignments/predicted_v_extracted_w/monomers/"
    image_output_wo = f"{global_output}/alignments/predicted_v_extracted_wo/monomers/"
    image_output_w_wo = f"{global_output}/alignments/predicted_w_v_wo/monomers/"
    image_output_w_orig = f"{global_output}/alignments/monomers_orig_v_w/"
    image_output_wo_orig = f"{global_output}/alignments/monomers_orig_v_w/"
    Path(image_output_w_orig).mkdir(parents=True, exist_ok=True)
    Path(image_output_wo_orig).mkdir(parents=True, exist_ok=True)
    print("RMSD ok")
    print(str(file_wo))
    print(str(file_w))
    # RMSD predicted vs extracted w
    rmsd_ex_w,tm_ex_w = align_structures(file_w, extracted_file_name_w,image_output_w , image_name=f"{name}.png")
    # RMSD predicted vs extracted wo
    rmsd_ex_wo,tm_ex_wo = align_structures(file_wo, extracted_file_name_wo,image_output_wo , image_name=f"{name}.png")
    # RMSD predicted w vs predicted wo
    rmsd_w_wo,tm_w_wo = align_structures(file_w, file_wo,image_output_w_wo , image_name=f"{name}.png")
    # RMSD predicted vs original w
    rmsd_og_w, tm_og_w = align_structures(file_w, file_orig,image_output_w_orig , image_name=f"{name}.png")
    # RMSD predicted vs original wo
    rmsd_og_wo,tm_og_wo = align_structures(file_wo, file_orig, image_output_wo_orig , image_name=f"{name}.png")
    # Solubility w
    seq_w = extract_seq_pdb(file_w)[0]
    seq_wo = extract_seq_pdb(file_wo)[0]
    sol_w = predict_solubility(name,seq_w)
    # Solubility wo
    sol_wo = predict_solubility(name,seq_wo)
    # Hydrophobicity w
    hydro_w = estimate_hydrophobicity(seq_w)
    # Hydrophobicity wo
    hydro_wo = estimate_hydrophobicity(seq_wo)
    # pLDDT from predicted w
    pLDDT_w = extract_plddt(file_w,17)
    # pLDDT from predicted wo
    pLDDT_wo = extract_plddt(file_wo,17)
    # Format row
    keys = ['file_ID','RMSD_pred-ex_w','Tm_pred-ex_w', 'RMSD_pred-ex_wo','Tm_pred-ex_wo','RMSD_w_wo','Tm_w_wo', 'RMSD_pred-og_w', 'Tm_pred-og_w', 'RMSD_pred-og_wo','Tm_pred-og_wo','solubility_w','solubility_wo','hydrophob_w','hydrophob_wo','pLDDT_w','pLDDT_wo']
    elements = [name + ".pdb", rmsd_ex_w,tm_ex_w,rmsd_ex_wo,tm_ex_wo, rmsd_w_wo,tm_w_wo, rmsd_og_w, tm_og_w, rmsd_og_wo, tm_og_wo,sol_w,sol_wo,hydro_w,hydro_wo, pLDDT_w[0],pLDDT_wo[0]]
    row = pd.Series(elements,index = keys)

    return pd.DataFrame([row])

def dimer_df(folder_w, folder_wo, global_output):
    df = pd.DataFrame()
    Path(f"{global_output}/foldx_out").mkdir(exist_ok=True, parents=True)
    folder_w_path = Path(folder_w)
    folder_wo_path = Path(folder_wo)

    # Extraer claves de los archivos de folder_w
    files_w_map = {}
    for f in folder_w_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_w_map[key] = f.name

    # Extraer claves de los archivos de folder_wo
    files_wo_map = {}
    for f in folder_wo_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_wo_map[key] = f.name

    # Intersección de claves
    common_keys = set(files_w_map.keys()) & set(files_wo_map.keys())

    # Eliminar archivos en folder_w que no estén en la intersección
    for key, fname in files_w_map.items():
        if key not in common_keys:
            (folder_w_path / fname).unlink()

    # Eliminar archivos en folder_wo que no estén en la intersección
    for key, fname in files_wo_map.items():
        if key not in common_keys:
            (folder_wo_path / fname).unlink()

    # Now iterate over the common files and process
    for file_w,file_wo in zip(Path(folder_w).iterdir(),Path(folder_wo).iterdir()):
        row = dimer_row(file_w, file_wo, global_output)
        df = pd.concat([df, row])

    return df

def dimer_row(file_w, file_wo,global_output):
    # Check that files exist
    keys = ['file_ID','interface_w','interface_wo','RMSD_pred-ex_w','Tm_pred-ex_w', 'RMSD_pred-ex_wo','Tm_pred-ex_wo', 'RMSD_w_wo','Tm_w_wo','pLDDT_w','pLDDT_wo']
    if not Path(file_w).exists() or not Path(file_wo).exists():
        return pd.Series([np.nan]*len(keys), index=keys)
    # Find name
    file_ex = Path(file_w)
    name = str(file_ex.name).split("_monomer.pdb")[0]
    file_ex = Path(file_wo)
    name_wo = str(file_ex.name).split("_monomer.pdb")[0]
    # Find extracted
    extracted_file_name_w = f"{global_output}/extracted_w/dimers/{name}_dimer.pdb"
    extracted_file_name_wo = f"{global_output}/extracted_wo/dimers/{name_wo}_dimer.pdb"
    extracted_file_name_wo =extracted_file_name_wo.split(".pdb")[0] + f".pdb_relaxed_capsid_{extracted_file_name_wo.split("_capsid_")[1][0]}_design_full_dimer.pdb"
    # Alignment folders 
    image_output_w = f"{global_output}/alignments/predicted_v_extracted_w/dimers/"
    image_output_wo = f"{global_output}/alignments/predicted_v_extracted_wo/dimers/"
    image_output_w_wo = f"{global_output}/alignments/predicted_w_v_wo/dimers/"
    # Interface w
    interface_w = interface_quality(file_w,f"{global_output}/foldx_out")
    # Interface wo
    interface_wo = interface_quality(file_wo,f"{global_output}/foldx_out")
    # RMSD predicted vs extracted w
    rmsd_ex_w, tm_ex_w = align_structures(file_w, extracted_file_name_w,image_output_w , image_name=f"{name}.png")
    # RMSD predicted vs extracted wo
    rmsd_ex_wo, tm_ex_wo = align_structures(file_wo, extracted_file_name_wo,image_output_wo , image_name=f"{name}.png")
    # RMSD predicted w vs predicted wo
    rmsd_w_wo, tm_w_wo = align_structures(file_w, file_wo,image_output_w_wo , image_name=f"{name}.png")
    # pLDDT from predicted w
    pLDDT_w = extract_plddt(file_w,17)
    # pLDDT from predicted wo
    pLDDT_wo = extract_plddt(file_wo,17)
    # Format row
    elements = [name + ".pdb", interface_w, interface_wo, rmsd_ex_w, tm_ex_w ,rmsd_ex_wo, tm_ex_wo , rmsd_w_wo, tm_w_wo, pLDDT_w[0],pLDDT_wo[0]]
    row = pd.Series(elements,index = keys)

    return pd.DataFrame([row])

def mer_df(folder_w, folder_wo,global_output,mode):
    df = pd.DataFrame()

    folder_w_path = Path(folder_w)
    folder_wo_path = Path(folder_wo)

    # Extraer claves de los archivos de folder_w
    files_w_map = {}
    for f in folder_w_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_w_map[key] = f.name

    # Extraer claves de los archivos de folder_wo
    files_wo_map = {}
    for f in folder_wo_path.iterdir():
        if f.is_file():
            key = f.name.split(".pdb")[0]
            files_wo_map[key] = f.name

    # Intersección de claves
    common_keys = set(files_w_map.keys()) & set(files_wo_map.keys())

    # Eliminar archivos en folder_w que no estén en la intersección
    for key, fname in files_w_map.items():
        if key not in common_keys:
            (folder_w_path / fname).unlink()

    # Eliminar archivos en folder_wo que no estén en la intersección
    for key, fname in files_wo_map.items():
        if key not in common_keys:
            (folder_wo_path / fname).unlink()

    for file_w,file_wo in zip(Path(folder_w).iterdir(),Path(folder_wo).iterdir()):
        row = mer_row(file_w,file_wo,global_output,mode)
        df = pd.concat([df,row])
    return df

def mer_row(file_w, file_wo, global_output, mode):
    # Find name
    file_ex = Path(file_w)
    name = str(file_ex.name).split("_monomer.pdb")[0]
    file_ex = Path(file_wo)
    name_wo = str(file_ex.name).split("_monomer.pdb")[0]
    # Find extracted
    extracted_file_name_w = f"{global_output}/extracted_w/{mode}s/{name}_{mode}.pdb"
    extracted_file_name_wo = f"{global_output}/extracted_wo/{mode}s/{name_wo}_{mode}.pdb"
    extracted_file_name_wo =extracted_file_name_wo.split(".pdb")[0] + f".pdb_relaxed_capsid_{extracted_file_name_wo.split("_capsid_")[1][0]}_design_full_{mode}.pdb"
    # Alignment folders 
    image_output_w = f"{global_output}/alignments/predicted_v_extracted_w/{mode}/"
    image_output_wo = f"{global_output}/alignments/predicted_v_extracted_wo/{mode}/"
    image_output_w_wo = f"{global_output}/alignments/predicted_w_v_wo/{mode}/"

    # RMSD predicted vs extracted w
    rmsd_ex_w,tm_ex_w = align_structures(file_w, extracted_file_name_w,image_output_w , image_name=f"{name}.png")
    # RMSD predicted vs extracted wo
    rmsd_ex_wo, tm_ex_wo = align_structures(file_wo, extracted_file_name_wo,image_output_wo , image_name=f"{name}.png")
    # RMSD predicted w vs predicted wo
    rmsd_w_wo,tm_w_wo = align_structures(file_w, file_wo,image_output_w_wo , image_name=f"{name}.png")
    # pLDDT from predicted w
    pLDDT_w = extract_plddt(file_w,17)
    # pLDDT from predicted wo
    pLDDT_wo = extract_plddt(file_wo,17)
    # Format row
    keys = ['file_ID','RMSD_pred-ex_w','Tm_pred-ex_w', 'RMSD_pred-ex_wo','Tm_pred-ex_wo', 'RMSD_w_wo','Tm_w_wo','pLDDT_w','pLDDT_wo']
    elements = [name + ".pdb", rmsd_ex_w,tm_ex_w,rmsd_ex_wo,tm_ex_wo,rmsd_w_wo,tm_w_wo, pLDDT_w[0],pLDDT_wo[0]]
    row = pd.Series(elements,index = keys)

    return pd.DataFrame([row])


### PLOTTING FUNCTIONS ##############################################################################################
monomer_keys = ['file_ID','stab_w','stab_wo','stab_og','RMSD_pred-ex_w', 'RMSD_pred-ex_wo','RMSD_w_wo', 'RMSD_pred-og_w', 'RMSD_pred-og_wo','pLDDT_w','pLDDT_wo']
dimer_keys = ['file_ID','interface_w','interface_wo','RMSD_pred-ex_w', 'RMSD_pred-ex_wo', 'RMSD_w_wo','pLDDT_w','pLDDT_wo']
mer_keys = ['file_ID','RMSD_pred-ex_w', 'RMSD_pred-ex_wo', 'RMSD_w_wo','pLDDT_w','pLDDT_wo']

def monomer_plots(monomer_df, plot_output):
    stab_palette = {
    "stab_w": "#a1dab4", 
    "stab_wo": "#41b6c4"
    }
    delta_palette = {
    'delta_delta_g_w': "#a1dab4", 
    'delta_delta_g_wo': "#41b6c4"
    }
    RMSD_palette = {
    'RMSD_pred-ex_w': "#66c2a5",  # teal-green
    'RMSD_pred-ex_wo': "#b2df8a", # light green
    'RMSD_w_wo': "#3288bd",       # stronger blue
    'RMSD_pred-og_w': "#abdda4"   # soft green
    }
    tm_palette = {
    'Tm_pred-ex_w': "#66c2a5",  # teal-green
    'Tm_pred-ex_wo': "#b2df8a", # light green
    'Tm_w_wo': "#3288bd",       # stronger blue
    'Tm_pred-og_w': "#abdda4"   # soft green
    }
    pLDDT_palette = {
    'pLDDT_w': "#fc8d62",   # soft orange
    'pLDDT_wo': "#8da0cb"   # lavender-blue
    }

    # Stability boxplot
    boxplot(monomer_df,
        ["stab_w","stab_wo"],"Stability (ΔG)",
        stab_palette,
        "Stability Comparison monomers","ΔG (kcal/mol)",
        f"{plot_output}/monomer_stab_boxplot.png",(6, 5))
    # Delta delta boxplot
    boxplot(monomer_df,
        ["delta_delta_g_w","delta_delta_g_wo"],"Stability (ΔG)",
        delta_palette,
        "Increase in stability (ΔΔG) monomers","ΔΔG (kcal/mol)",
        f"{plot_output}/monomer_delta_delta_boxplot.png",(6, 5))
    #Stability paired lines
    plot_paired(monomer_df,
                "stab_wo","stab_w",
                "Stability (ΔG) in monomers",
                "ΔG (kcal/mol))",
                f"{plot_output}/monomer_stab_paired_lines.png")
    # delta delta paired lines
    plot_paired(monomer_df,
                "delta_delta_g_wo","delta_delta_g_w",
                "Increase in stability (ΔΔG) monomers",
                "ΔG (kcal/mol))",
                f"{plot_output}/monomer_delta_delta_paired_lines.png")
    # RMSD boxplot
    boxplot(monomer_df,
        ['RMSD_pred-ex_w', 'RMSD_pred-ex_wo','RMSD_w_wo', 'RMSD_pred-og_w'],"RMSD",
        RMSD_palette,
        "RMSD comparison","Arsmtrongs",
        f"{plot_output}/monomer_RMSD_boxplot.png",(20, 8))
    #RMSD paired lines
    plot_paired(monomer_df,
                'RMSD_pred-ex_wo','RMSD_pred-ex_w',
                "RMSD monomers",
                "Armstrongs",
                f"{plot_output}/monomer_RMSD_paired_lines.png")
    # Tm boxplot
    boxplot(monomer_df,
        ['Tm_pred-ex_w', 'Tm_pred-ex_wo','Tm_w_wo', 'Tm_pred-og_w'],"Tm-score",
        tm_palette,
        "Tm score alignment comparison","Tm score",
        f"{plot_output}/monomer_tm_boxplot.png",(20, 8))
    #Tm paired lines
    plot_paired(monomer_df,
                'Tm_pred-ex_wo','Tm_pred-ex_w',
                "Tm score monomers",
                "Tm score",
                f"{plot_output}/monomer_tm_paired_lines.png")
    # pLDDT boxplot
    boxplot(monomer_df,
        ['pLDDT_w','pLDDT_wo'],"pLDDT",
        pLDDT_palette,
        "Stability Comparison (stab_w vs. stab_wo)","ΔG (kcal/mol)",
        f"{plot_output}/monomer_pLDDT_boxplot.png",(6, 5))
    
    # pLDDT paired lines
    plot_paired(monomer_df,
                "pLDDT_wo","pLDDT_w",
                "pLDDT in monomers",
                "ΔG (kcal/mol))",
                f"{plot_output}/monomer_pLDDT_paired_lines.png")
    return

def dimer_plots(dimer_df):
    # Interface baxplot
    return

def mer_plots(df,mode):
    return

def combined_plots(mon_df,dimer_df,trimer_df,pentamer_df):
    return


def boxplot(df,col_list,value,palette,tittle,y_lab,output_name,fig_size):
    # Prepare long-format data
    stab_long = df[col_list].melt(var_name="Condition", value_name=value)

    # Set style
    sns.set(style="whitegrid", context="notebook")

    # Define a color palette that you can re-use

    plt.figure(figsize=fig_size)

    # Boxplot with mean marker
    ax = sns.boxplot(
        x="Condition",
        y=value,
        data=stab_long,
        width=0.4,
        palette=palette,
        showmeans=True,
        meanprops={
            "marker": "o",
            "markerfacecolor": "black",
            "markeredgecolor": "black",
            "markersize": 6
        },
        zorder = 2
    )

    # Overlay individual points, colored by condition
    sns.stripplot(
        x="Condition",
        y=value,
        data=stab_long,
        jitter=True,
        dodge=False,
        palette=palette,
        alpha=0.7,
        zorder = 1
    )

    # Annotate the mean values as text
    mean_values = stab_long.groupby("Condition")[value].mean()
    for i, (condition, mean) in enumerate(mean_values.items()):
        ax.text(
            i, 
            mean + 0.3,  # slightly above the marker
            f"{mean:.2f}", 
            ha='center', 
            va='bottom', 
            fontsize=10, 
            color='black',
            zorder = 4
        )
    # Add breathing space above and below
    y_min = stab_long[value].min()
    y_max = stab_long[value].max()
    y_margin = (y_max - y_min) * 0.1
    plt.ylim(y_min - y_margin, y_max + y_margin)
    # Final touches
    plt.title(tittle)
    plt.ylabel(y_lab)
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(output_name)
    return

def plot_paired(df, col_wo, col_w, tittle, xlabel, output):
    # Order by col_w
    df = df.sort_values(by=col_w, ascending=False).reset_index(drop=True)

    # Generate cleaned file_IDs for y-axis labels
    file_IDs = []
    for f in df['file_ID']:
        base = f.split('.pdb')[0]
        nums = re.findall(r'\d+', f.split('.pdb', 1)[1])  # extract numbers after .pdb
        new_id = base + '_' + '_'.join(nums)
        file_IDs.append(new_id)

    # Use numeric y positions for plotting
    y_pos = range(len(df))

    # Start plotting
    plt.figure(figsize=(10, 18))
    for i, row in df.iterrows():
        plt.plot([row[col_wo], row[col_w]], [i, i], 'gray', alpha=0.4)

    plt.scatter(df[col_wo], y_pos, color='red', label='Without RF')
    plt.scatter(df[col_w], y_pos, color='blue', label='With RF')

    plt.yticks(ticks=y_pos, labels=file_IDs)
    plt.xlabel(xlabel)
    plt.ylabel("Sample")
    plt.title(tittle)
    plt.legend()
    plt.grid(True, axis='x', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output)
    return

def plot_paired2(df,col_wo,col_w,tittle,xlabel,output):
    # Order by col1
    df = df.sort_values(by=col_w,ascending=False)
    # Generate y axis
    file_IDs = list(df['file_ID'])
    y = []
    for f in file_IDs:
        base = f.split('.pdb')[0]  # before the first .pdb
        nums = re.findall(r'\d+', f.split('.pdb', 1)[1])  # numbers after first .pdb
        new_id = base + '_' + '_'.join(nums)
        y.append(new_id)
    # Start plotting
    plt.figure(figsize=(10, 18))
    for i, row in df.iterrows():
        plt.plot([row[col_wo], row[col_w]], [y[i], y[i]], 'gray', alpha=0.4)

    plt.scatter(df[col_wo], y, color='red', label='Without RF')
    plt.scatter(df[col_w], y, color='blue', label='With RF')

    plt.yticks(ticks=range(len(y)), labels=y)
    plt.xlabel(xlabel)
    plt.ylabel("Sample")
    plt.title(tittle)
    plt.legend()
    plt.grid(True, axis='x', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output)
    return

### STATISTICS ##################################################################################################################################

def statistics_columns(df, col1, col2, paired=False, alpha=0.05):
    x = df[col1].dropna()
    y = df[col2].dropna()
    
    # Ensure same length for paired tests
    if paired and len(x) != len(y):
        raise ValueError("For paired tests, both columns must have the same number of observations")

    # Normality tests
    p_norm_x = stats.shapiro(x)[1]
    p_norm_y = stats.shapiro(y)[1]
    normal = (p_norm_x > alpha) and (p_norm_y > alpha)

    # Determine which test to run
    if paired:
        if normal:
            test_stat, p = stats.ttest_rel(x, y)
            test_name = "Paired t-test"
            effect_size = (x - y).mean() / (x - y).std(ddof=1)  # Cohen's d for paired
        else:
            test_stat, p = stats.wilcoxon(x, y)
            test_name = "Wilcoxon signed-rank"
            effect_size = (x > y).mean() - (x < y).mean()  # Rank biserial
    else:
        # Check variance homogeneity
        equal_var = stats.levene(x, y).pvalue > alpha

        if normal:
            test_stat, p = stats.ttest_ind(x, y, equal_var=equal_var)
            test_name = "Unpaired t-test" if equal_var else "Welch’s t-test"
            # Pooled standard deviation
            pooled_std = np.sqrt(((len(x)-1)*np.std(x, ddof=1)**2 + (len(y)-1)*np.std(y, ddof=1)**2) / (len(x) + len(y) - 2))
            effect_size = (np.mean(x) - np.mean(y)) / pooled_std  # Cohen's d
        else:
            test_stat, p = stats.mannwhitneyu(x, y, alternative='two-sided')
            test_name = "Mann–Whitney U"
            # Rank biserial approximation
            n1, n2 = len(x), len(y)
            U = test_stat
            effect_size = 1 - (2 * U) / (n1 * n2)
    
    if p > 0.05 :
        star = "ns"
    elif p < 0.05 and p > 0.01 :
        star = "*"
    elif p < 0.01 and p > 0.001 :
        star = "**"
    else:
        star = "***"
    
    return {
        "test": test_name,
        "p_value": p,
        "star": star,
        "statistic": test_stat,
        "mean_difference": np.mean(x) - np.mean(y),
        "effect_size": effect_size,
        "normality": {"col1": p_norm_x, "col2": p_norm_y},
        "equal_variance": equal_var if not paired else None
    }

def monomer_statistics(df):
    # Set up comparison dictionary

    # Initialise output dataframe
    final_df = pd.DataFrame(index=['Comparison','stars','p-value','mean_difference','test used'])
    # Make comparisons

    return
### AUXILIARY FUNCTIONS #############################################################################################
def help():
    pass

def check_matching(folder_w, folder_wo):
    # Prepare output folders
    unmatched_w = os.path.join(folder_w, 'unmatched_w')
    unmatched_wo = os.path.join(folder_wo, 'unmatched_wo')
    os.makedirs(unmatched_w, exist_ok=True)
    os.makedirs(unmatched_wo, exist_ok=True)

    # Extract all file prefixes before '.pdb'
    w_prefixes = set()
    for f in os.listdir(folder_w):
        prefix = f.split('.pdb')[0]
        w_prefixes.add(prefix)

    wo_prefixes = set()
    for f in os.listdir(folder_wo):
        prefix = f.split('.pdb')[0]
        wo_prefixes.add(prefix)

    # Identify unique prefixes in each folder
    only_in_w = w_prefixes - wo_prefixes
    only_in_wo = wo_prefixes - w_prefixes

    # Move files with unmatched prefixes from folder_w
    for base in only_in_w:
        for fname in os.listdir(folder_w):
            if fname.startswith(base + '.pdb'):
                src = os.path.join(folder_w, fname)
                dst = os.path.join(unmatched_w, fname)
                os.rename(src, dst)
                print(f"Moved {fname} to {unmatched_w}")

    # Move files with unmatched prefixes from folder_wo
    for base in only_in_wo:
        for fname in os.listdir(folder_wo):
            if fname.startswith(base + '.pdb'):
                src = os.path.join(folder_wo, fname)
                dst = os.path.join(unmatched_wo, fname)
                os.rename(src, dst)
                print(f"Moved {fname} to {unmatched_wo}")

    return

def rename_files_to_lowercase(folder_path="."):
    for filename in os.listdir(folder_path):
        src = os.path.join(folder_path, filename)
        if os.path.isfile(src):
            new_name = filename.lower()
            dst = os.path.join(folder_path, new_name)
            if src != dst:
                os.rename(src, dst)

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


def organise_chai_folder(root_folder, destination_folder):
    """
    Uses pathlib to copy `pred.model_idx_0.cif` from subfolders ending in
    monomer, dimer, trimer, or pentamer, and renames the file to the name of its folder.

    Parameters:
    - root_folder: str or Path to the root input directory
    - destination_folder: str or Path to the output directory
    """

    root_folder = Path(root_folder)
    destination_folder = Path(destination_folder)

    # Create destination folders
    types = ["monomer", "dimer", "trimer", "pentamer"]
    for t in types:
        (destination_folder / f"{t}s").mkdir(parents=True, exist_ok=True)

    # Traverse sub-subfolders
    for subfolder in root_folder.glob("*/*"):
        if subfolder.is_dir():
            for t in types:
                if subfolder.name.endswith(t):
                    source_file = subfolder / "pred.model_idx_0.cif"
                    if source_file.exists():
                        # Rename output file to subfolder name (e.g. dimer.cif)
                        renamed_file = destination_folder / f"{t}s" / f"{subfolder.name}.cif"
                        renamed_file.write_bytes(source_file.read_bytes())
                        print(f"Copied {source_file} → {renamed_file}")
                    break

def convert_all_cif_to_pdb(input_folder, output_folder):
    command = f"conda run -n benchmark python scripts/pdb_converter.py {input_folder} {output_folder}"
    process = subprocess.run(command, shell=True)
