#
"""
improved_shapedesign.py input_folder --help[OPTIONAL]


IMPORTANT: All input files must be pdbs and must not contain points (.) other than for the extension
"""
### IMPORT MODULES ###########################################################################################################################
import functions as func
import pandas as pd
import argparse
from pathlib import Path
import sys
import logging

### INPUT CHECK ###############################################################################################################################
parser = argparse.ArgumentParser(
    description="This script runs the pipeline shapedesign with an additional RF_diffusion step"
)
parser.add_argument('--folder', help="Folder that contains all the input structures", type=str)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# If --detailed-help is provided, print custom help and exit
if args.detailed_help:
    func.help()
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
global_output_name=Path("improved_shapedesign_output_1")
count=1
while global_output_name.exists():
    global_output_name = Path(f"{global_output_name}_{count}")
    count += 1
global_output_name.mkdir(parents=True, exist_ok=True)

#Generate the rest of the outputs
Path(f"{str(global_output_name)}/rf_diffusion/scaffolds").mkdir(parents=True, exist_ok=True)
Path(f"{str(global_output_name)}/rf_diffusion/pdbs").mkdir(parents=True, exist_ok=True)
Path(f"{str(global_output_name)}/MPNN").mkdir(parents=True, exist_ok=True)
Path(f"{str(global_output_name)}/ESM").mkdir(parents=True, exist_ok=True)

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

save_to_log("Starting execution")
# 0. Format the names of all the input_files, so that they are all lowercase ###########################################
save_to_log("Processing input files")
func.rename_files_to_lowercase(input_folder)

### 1.Load the input folder as a df (["file_ID","sequence","length", "stability"]) Stability will be used later

original_structures_df = pd.DataFrame(columns=["file_ID","seq","length","sol_orig","hydro_orig"])

for file in input_folder.iterdir():
    file_name = file.name
    seq,length = func.extract_seq_pdb(file)
    sol_orig = func.predict_solubility(file_name,seq)
    hydro_orig = func.estimate_hydrophobicity(seq)
    row = pd.Series([file_name,seq,length,sol_orig,hydro_orig],index=["file_ID","seq","length","sol_orig","hydro_orig"])
    original_structures_df.loc[len(original_structures_df)] = row

original_structures_df=func.calculate_stability(input_folder,original_structures_df,f"{global_output_name}")
original_structures_df.to_csv(f"{global_output_name}/original_structures_def.csv")

save_to_log(f"{len(original_structures_df)} structures loaded")

### 2. RF diffusion step. For each row of the data frame run a partial diffussion and a fold conditioning
save_to_log("Starting RF diffussion step")
for index,row in original_structures_df.iterrows():
    func.run_rf_scaffold(row,"input_pdbs/",f"{str(global_output_name)}/rf_diffusion")
    func.run_rf_partial(row,"input_pdbs/",f"{str(global_output_name)}/rf_diffusion")


# Move all the PDBs to a sepparate folder
# Define source and destination folders
src_folder = Path(f"{str(global_output_name)}/rf_diffusion")
dst_folder = Path(f"{str(global_output_name)}/rf_diffusion/pdbs")
# Move each .pdb file
for pdb_file in src_folder.glob("*.pdb"):
    new_location = dst_folder / pdb_file.name
    pdb_file.rename(new_location)

save_to_log("Finished RF diffussion step")
### 3. MPNN step

save_to_log("Starting MPNN step")
# Generate sequences for each of the pdbs obtained from RF diffussion

func.run_MPNN(f"{str(global_output_name)}/rf_diffusion/pdbs",f"{str(global_output_name)}/MPNN")

# Save the generated sequences to a dataframe in order to filter them
MPNN_df = func.process_MPNN_folder(f"{str(global_output_name)}/MPNN/seqs")
save_to_log("Finished MPNN step")

### 4. 1D params and filter for MPNN seqs
filtered_MPNN_df = func.oned_dataframe_filter(MPNN_df)
n_filtered = len(MPNN_df) - len(filtered_MPNN_df)
save_to_log(f"Filtered {n_filtered} sequences according to 1D filters")

### 5. ESM prediction
save_to_log("Starting ESM step")

esm_output = f"{str(global_output_name)}/ESM"
for index,row in filtered_MPNN_df.iterrows():
    func.run_ESM(row,esm_output)

# Remove failed predictions (files under 1K size)
func.remove_small_files(esm_output,1024)

save_to_log("Finished ESM step")

### 6. 3D params for ESM prediction
save_to_log("Calculating 3D parameters for ESM predictions")

esm_3d_df = func.threed_params_df(esm_output,"input_pdbs",original_structures_df,pLDDT_loc=10)
esm_3d_df.to_csv(f"{str(global_output_name)}/ESM_3D.csv")

save_to_log("Finished calculating 3D parameters for ESM predictions")

### 7. 3D filtering for ESM prediction
esm_3d_df = pd.read_csv(f"{str(global_output_name)}/ESM_3D.csv")
save_to_log("Filtering ESM predictions")
filtered_esm_3d_df = func.threed_filter_df(esm_3d_df)

n_filtered = len(esm_3d_df) - len(filtered_esm_3d_df)
save_to_log(f"Filtered {n_filtered} structures according to 3D filters")

### 8. Chai prediction

save_to_log("Starting Chai step")

chai_output = f"{str(global_output_name)}/Chai"
Path(chai_output).mkdir(parents= True, exist_ok=True)
for index,row in filtered_esm_3d_df.iterrows():
    func.run_chai(row,chai_output)
save_to_log("Finishing Chai step")


# Format the chai folder adequately
chai_final_output =f"{str(global_output_name)}/Chai_cif_files"
chai_pdbs = f"{str(global_output_name)}/Chai_pdb_files"
#func.format_chai_folder(chai_output,chai_final_output)
#func.convert_all_cif_to_pdb(chai_final_output,chai_pdbs)
original_structures_df=pd.read_csv(f"{str(global_output_name)}/original_structures.csv")
### 9. 3D params for Chai
save_to_log("Calculating 3D parameters for Chai predictions")

#chai_3d_df = func.threed_params_df(chai_pdbs,"input_pdbs",original_structures_df,pLDDT_loc=10)
#chai_3d_df.to_csv(f"{str(global_output_name)}/Chai_3D.csv")

save_to_log("Finished calculating 3D parameters for Chai predictions")

### 10. 3D filtering for Chai
save_to_log("Filtering Chai predictions")
chai_3d_df = pd.read_csv(f"{str(global_output_name)}/Chai_3D.csv")
file_list,filtered_chai_3d_df = func.final_3d_filter(chai_3d_df,original_structures_df,global_output_name)
filtered_chai_3d_df.to_csv(f"{str(global_output_name)}/final_filter_small.csv")
n_filtered = len(chai_3d_df) - len(filtered_chai_3d_df)
n_dissapeared = len(original_structures_df) - len(file_list)

save_to_log(f"Filtered {n_filtered} structures according to 3D filters")
save_to_log(f"There are succesful monomers for {len(file_list)} original structures.")
save_to_log(f"{n_dissapeared} original structures have been dropped in the process.")


### 11. Minimise with pyrosetta
save_to_log("Starting structure minimisation")

minimisation_input_path = f"{str(global_output_name)}/minimisation/inputs"
minimised_folder = f"{str(global_output_name)}/minimisation/minimised_structures"
Path(minimisation_input_path).mkdir(parents=True, exist_ok=True)
Path(minimised_folder).mkdir(parents=True, exist_ok=True)

# Move the saved structures to a single folder
func.folder_from_dataframe(filtered_chai_3d_df,chai_pdbs,minimisation_input_path)
# Parse that folder minimising each structure

func.minimise_structure(minimisation_input_path,minimised_folder)

save_to_log("Finished structure minimisation")
### 12. Shapedesign iterations

save_to_log("Starting shapedesign iterations")
shapedesign_output = f"{str(global_output_name)}/shapedesign2"
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
        func.run_shapedesign(new_file_name, shapedesign_output)

        # Append metrics for each generated capsid
        indiv_capsid_metrics = func.capsid_metrics(
            indiv_capsid_metrics,
            new_file_name.name,
            f"{str(shapedesign_output)}/design.csv"
        )

    # Filter the best capsid from those generated for this file
    selected_capsid = func.filter_capsids(indiv_capsid_metrics, shapedesign_output,final_capsids)

    # Add selected capsid to the global dataframe
    global_capsid_metrics = pd.concat([global_capsid_metrics, selected_capsid], ignore_index=True)

    # Reinitialize individual metrics DataFrame
    indiv_capsid_metrics = pd.DataFrame()



global_capsid_metrics.to_csv(f"{str(global_output_name)}/capsid_metrics_2.csv")

func.move_files(f"{str(global_output_name)}/shapedesign/structures", "96_capsids_w","full.cif")


# Generate 25 capsids (5 times 5) for a random input

save_to_log("Finished shapedesign iterations")
save_to_log("FINISHED PIPELINE EXECUTION")
