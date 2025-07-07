import os
import sys

def rename_files_to_lowercase(folder_path="."):
    for filename in os.listdir(folder_path):
        src = os.path.join(folder_path, filename)
        if os.path.isfile(src):
            new_name = filename.lower()
            dst = os.path.join(folder_path, new_name)
            if src != dst:
                os.rename(src, dst)
                print(f"Renamed: {filename} → {new_name}")

if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    rename_files_to_lowercase(folder)