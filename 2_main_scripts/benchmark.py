#
"""
benchmark.py --original_pdbs --w --wo --detailed-help

This script is designed to benchmark the improved shapedesign pipeline by comparing its output capsids with capsids generated without shapedesign. The script:

A. Checks the arguments and input folders
B. Files are renamed to lowercase and check that each sample has a corresponding with and without version
C. Extracts of monomers and oligomers from capsids with and without using pymol
D. Organises of extracted complexes in subfolders
E. Uses extracted monomers for multimer prediction with Chai1
F. Predictions are organised in subfolders
G. Different metrics are calculated for monomers, dimers, trimers and pentamers

The behaviour of each function is detailed in the functions_benchmark.py script
"""
####################################################################################################################
### IMPORT MODULES            ######################################################################################
####################################################################################################################
import functions_benchmark as func
import pandas as pd
import argparse
from pathlib import Path
import sys
import logging
import numpy as np
####################################################################################################################
### INPUT CHECK               ######################################################################################
####################################################################################################################

parser = argparse.ArgumentParser(
    description="This script runs the pipeline shapedesign with an additional RF_diffusion step"
)
parser.add_argument('--original_pdbs', help="Folder containing monomers BEFORE shapedesign", type=str)
parser.add_argument('--w', help="Folder containing capsids generated WITH an RF diffusion step", type=str)
parser.add_argument('--wo', help="Folder containing capsids generated WITHOUT an RF diffusion step", type=str)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# If --detailed-help is provided, print custom help and exit
if args.detailed_help:
    func.help()
    sys.exit()

# If --w or --wo was not provided, show error and exit
if not args.w or not args.wo :
    parser.error("--w and --wo folders are required unless --detailed-help is used")

if not args.original_pdbs :
    parser.error("--original_pdbs folder is required unless --detailed-help is used") 

# Convert to Path objects

original_pdbs = Path(args.original_pdbs)
w_folder = Path(args.w)
wo_folder = Path(args.wo)

# If original_pdbs is not an existing directory, show error and exit
if not original_pdbs.exists():
    parser.error(f"The folder '{original_pdbs}' does not exist.")
if not original_pdbs.is_dir():
    parser.error(f"The path '{original_pdbs}' is not a directory.")

# If w_folder is not an existing directory, show error and exit
if not w_folder.exists():
    parser.error(f"The folder '{w_folder}' does not exist.")
if not w_folder.is_dir():
    parser.error(f"The path '{w_folder}' is not a directory.")

# If wo_folder is not an existing directory, show error and exit
if not wo_folder.exists():
    parser.error(f"The folder '{wo_folder}' does not exist.")
if not wo_folder.is_dir():
    parser.error(f"The path '{wo_folder}' is not a directory.")



####################################################################################################################
### SET UP OUTPUT             ######################################################################################
####################################################################################################################

# Generate the global folder, making sure that it doesn't exist already by several iterations

global_output_name=Path("benchmark_output")
count=1
while global_output_name.exists():
    global_output_name = Path(f"{global_output_name}_{count}")
    count += 1

global_output_name.mkdir(parents=True, exist_ok=True)

####################################################################################################################
### SET UP LOG FILE           ######################################################################################
####################################################################################################################

logging.basicConfig(
    filename=f"{str(global_output_name)}/shapedesign_benchmark.log",       # Log file path
    filemode='a',                   # Append mode
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)

def save_to_log(message):
    logging.info(message)
####################################################################################################################
### MAIN EXECUTION            ######################################################################################
####################################################################################################################

save_to_log("INITIALISING BENCHMARKING")

# 0. Format input files to lowercase and check matching ###########################################
save_to_log("Processing input files")

func.rename_files_to_lowercase(original_pdbs)
func.rename_files_to_lowercase(w_folder)
func.rename_files_to_lowercase(wo_folder)
func.check_matching(w_folder, wo_folder) # Check that the folders have matching samples and remove those that don't match

# 1. Extraction of monomers, dimers, trimers, pentamers ############################################################
save_to_log("Extracting substructures from the input capsids")

# Create extracted folders
extracted_w_folder = f"{global_output_name}/extracted_w"
extracted_wo_folder = f"{global_output_name}/extracted_wo"
Path(extracted_w_folder).mkdir(exist_ok=True,parents=True)
Path(extracted_wo_folder).mkdir(exist_ok=True,parents=True)

func.extraction(w_folder,extracted_w_folder)
func.extraction(wo_folder,extracted_wo_folder)

save_to_log("Finished extraction")

# 2. Prediction of monomers, dimers, trimers, pentamers ############################################################
save_to_log("Starting prediction of substructures")

# Create prediction folder
prediction_w_folder = f"{global_output_name}/prediction_w"
prediction_wo_folder = f"{global_output_name}/prediction_wo"
Path(prediction_w_folder).mkdir(exist_ok=True,parents=True)
Path(prediction_wo_folder).mkdir(exist_ok=True,parents=True)

func.run_chai_folder(f"{extracted_w_folder}/monomers",prediction_w_folder) #Parse the monomer folder to extract a single sequence
func.run_chai_folder(f"{extracted_wo_folder}/monomers",prediction_wo_folder)

save_to_log("Finished prediction of substructures")

# Reorganise predictions into a clearer structure
save_to_log("Formating prediction folders")

organised_prediction_w_folder = f"{global_output_name}/organised_prediction_w"
organised_prediction_wo_folder = f"{global_output_name}/organised_prediction_wo"
Path(organised_prediction_w_folder).mkdir(exist_ok=True,parents=True)
Path(organised_prediction_wo_folder).mkdir(exist_ok=True,parents=True)

func.organise_chai_folder(prediction_w_folder,organised_prediction_w_folder)
func.organise_chai_folder(prediction_wo_folder,organised_prediction_wo_folder)

# 3. Metrics calculation ###########################################################################################
save_to_log("Starting metrics calculation")

#Generate alignment outputs for visualisation
folders = ["predicted_v_extracted_w", "predicted_v_extracted_wo","predicted_w_v_wo"]
subfolders = ["monomers","dimers","trimers","pentamers"]
for folder in folders:
    for subfolder in subfolders:
        Path(f"{global_output_name}/alignments/{folder}/{subfolder}").mkdir(exist_ok=True, parents=True)
        
# 3.1 Monomers #####################################################################################################

monomer_metrics = func.monomer_df(f"{organised_prediction_w_folder}/monomers", f"{organised_prediction_wo_folder}/monomers", global_output_name,str(original_pdbs))
monomer_metrics.to_csv(f"{global_output_name}/monomer_metrics_2.csv")
save_to_log("Finished monomers metrics calculation")

# 3.2 Dimers #######################################################################################################

dimer_metrics = func.dimer_df(f"{organised_prediction_w_folder}/dimers", f"{organised_prediction_wo_folder}/dimers", global_output_name)
dimer_metrics.to_csv(f"{global_output_name}/dimer_metrics.csv")
save_to_log("Finished dimers metrics calculation")

# 3.3 Trimers ######################################################################################################

trimer_metrics = func.mer_df(f"{organised_prediction_w_folder}/trimers", f"{organised_prediction_wo_folder}/trimers", global_output_name, "trimer")
trimer_metrics.to_csv(f"{global_output_name}/trimer_metrics.csv")
save_to_log("Finished trimers metrics calculation")
# 3.4 Pentamers ####################################################################################################

pentamer_metrics = func.mer_df(f"{organised_prediction_w_folder}/pentamers", f"{organised_prediction_wo_folder}/pentamers", global_output_name, "pentamer")
pentamer_metrics.to_csv(f"{global_output_name}/pentamer_metrics.csv")
save_to_log("Finished pentamers metrics calculation")

# 3.5 Capsids ######################################################################################################


### PLOTTING #######################################################################################################

plot_output = f"{global_output_name}/plots"
Path(plot_output).mkdir(parents=True, exist_ok=True)
# Define keys
monomer_keys = ['file_ID','stab_w','stab_wo','stab_og','RMSD_pred-ex_w', 'RMSD_pred-ex_wo','RMSD_w_wo', 
                'RMSD_pred-og_w', 'RMSD_pred-og_wo','pLDDT_w','pLDDT_wo']

dimer_keys = ['file_ID','interface_w','interface_wo','RMSD_pred-ex_w', 'RMSD_pred-ex_wo', 
              'RMSD_w_wo','pLDDT_w','pLDDT_wo']

mer_keys = ['file_ID','RMSD_pred-ex_w', 'RMSD_pred-ex_wo', 'RMSD_w_wo','pLDDT_w','pLDDT_wo']

monomer_df = pd.read_csv(f"{global_output_name}/monomer_metrics.csv")
# Start plotting
func.monomer_plots(monomer_df,plot_output)
