import os
from glob import glob
from pymol import cmd

def convert_all_cif_to_pdb(input_folder, output_folder):
    """
    Convert all .cif files in the input folder to .pdb using PyMOL.
    """
    os.makedirs(output_folder, exist_ok=True)

    for cif_file in glob(os.path.join(input_folder, "*.cif")):
        name = os.path.splitext(os.path.basename(cif_file))[0]
        try:
            cmd.load(cif_file, name)
            output_path = os.path.join(output_folder, f"{name}.pdb")
            cmd.save(output_path, name)
            cmd.delete(name)
        except Exception as e:
            print(f"❌ Failed to convert {name}: {e}")

# Only run when executed directly, not on import
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: pymol -cq convert_cif_to_pdb.py -- /path/to/input /path/to/output")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    convert_all_cif_to_pdb(input_folder, output_folder)