import os
from glob import glob
from pymol import cmd

def convert_all_cif_to_pdb(input_folder, output_folder):
    """
    Convert all .cif files in the input folder to .pdb using PyMOL.
    Prints debug info to help trace unexpected file deletions.
    """
    os.makedirs(output_folder, exist_ok=True)
    print(f"📄 Number of files in {input_folder}: {len(os.listdir(input_folder))}")
    print(f"📄 Number of files in {output_folder}: {len(os.listdir(output_folder))}")

    cif_files = [f for f in os.listdir(input_folder) if f.endswith(".cif") or ".cif" in f]
    print(f"📂 Found {len(cif_files)} CIF files in: {input_folder}")
    print("🔍 Files before conversion:")
    for f in cif_files:
        print(f"  - {os.path.basename(f)}")

    for cif_file in cif_files:
        name = os.path.splitext(os.path.basename(cif_file))[0]
        cif_file = f"{input_folder}/{cif_file}"
        print(f"\n🔄 Attempting to convert: {name}")
        try:
            # Check file existence before load
            if not os.path.exists(cif_file):
                print(f"⚠️ File not found before loading: {cif_file}")
                continue

            cmd.load(cif_file, name)
            output_path = os.path.join(output_folder, f"{name}.pdb")
            cmd.save(output_path, name)
            cmd.delete(name)

            print(f"✅ Successfully saved: {output_path}")

        except Exception as e:
            print(f"❌ Failed to convert {name}: {e}")

        # Check if file still exists after operation
        if not os.path.exists(cif_file):
            print(f"⚠️ File {cif_file} is MISSING after conversion!")
        else:
            print(f"🟢 File {cif_file} is still present.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python pdb_converter.py /path/to/input /path/to/output")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    convert_all_cif_to_pdb(input_folder, output_folder)