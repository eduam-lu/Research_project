"""
Given a structure path, it returns a minimised a minimised structure in the indicated output path
"""
### IMPORT MODULES
import pyrosetta
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.idealize import IdealizeMover
import os
import argparse

# Initialize PyRosetta
pyrosetta.init("-ignore_unrecognized_res -load_PDB_components false")  # useful for CIF support

### ARGPARSE SETUP
parser = argparse.ArgumentParser(
    description="Relax a CIF structure or all CIFs in a folder. Outputs relaxed structures to a given folder."
)
parser.add_argument('--input_structure', help="Path to input .cif file or folder with .cif files", type=str, required=True)
parser.add_argument('--output_folder', help="Folder to output the relaxed structures", type=str, required=True)

args = parser.parse_args()
input_path = args.input_structure
output_folder = args.output_folder

### RELAX FUNCTION FOR A SINGLE FILE
def relax_single(pdb_path, output_dir):
    """
    Relax a structure and write to output.
    """
    os.makedirs(output_dir, exist_ok=True)
    file_name = os.path.basename(pdb_path)
    base_name = os.path.splitext(file_name)[0]
    output_path = os.path.join(output_dir, f"{base_name}_relaxed.cif")

    try:
        pose = pyrosetta.pose_from_file(pdb_path)
        IdealizeMover().apply(pose)
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = FastRelax(scorefxn, 15)
        relax.apply(pose)  
        pose.dump_cif(output_path)
        #print(f"[✓] Relaxed saved: {output_path}")
    except Exception as e:
        print(f"[!] Failed on {file_name}: {e}")
        failed_path = output_path.replace("relaxed", "failed")
        pose.dump_pdb(failed_path)

### FUNCTION TO HANDLE DIRECTORY
def relax_folder(input_dir, output_dir):
    """
    Relax all .cif files in the given folder.
    """
    cif_files = [f for f in os.listdir(input_dir) if f.endswith(('.cif','.pdb'))]
    if not cif_files:
        print("[!] No .cif files found.")
        return

    for cif in cif_files:
        full_path = os.path.join(input_dir, cif)
        relax_single(full_path, output_dir)

### MAIN LOGIC
if os.path.isfile(input_path) and input_path.endswith('.cif'):
    relax_single(input_path, output_folder)
elif os.path.isdir(input_path):
    relax_folder(input_path, output_folder)
else:
    print(f"[!] Input must be a .cif file or a folder containing .cif files. Got: {input_path}")