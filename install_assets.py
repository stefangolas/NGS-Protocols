import shutil
from pathlib import Path

# Overwrite options (global)
overwrite_all = False
skip_all = False

def prompt_overwrite(file_path):
    global overwrite_all, skip_all
    if overwrite_all:
        return True
    if skip_all:
        return False

    while True:
        choice = input(f"File exists: {file_path}. Overwrite? [y]es / [n]o / [a]ll / [s]kip all: ").lower()
        if choice == 'y':
            return True
        elif choice == 'n':
            return False
        elif choice == 'a':
            overwrite_all = True
            return True
        elif choice == 's':
            skip_all = True
            return False
        else:
            print("Invalid choice, please enter y, n, a, or s.")

def copy_item(src, dest):
    if src.is_dir():
        dest.mkdir(exist_ok=True)
        for sub_item in src.iterdir():
            copy_item(sub_item, dest / sub_item.name)
    else:
        if dest.exists():
            if not prompt_overwrite(dest):
                return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)

def main():
    # Disclaimer
    print("WARNING: This script will safely copy the contents of the 'assets' folder")
    print("into your Hamilton directory at C:\\Program Files (x86)\\Hamilton.")
    print("If this folder already exists, files may be overwritten depending on your input.")
    print("It is strongly recommended that you BACK UP this directory before continuing.")
    input("Press Enter to acknowledge and continue...")

    # Paths
    source = Path("assets")
    destination = Path(r"C:\Program Files (x86)\Hamilton")

    # Ensure destination exists
    destination.mkdir(parents=True, exist_ok=True)

    # Start copy
    for item in source.iterdir():
        copy_item(item, destination / item.name)

    print("Copy completed successfully.")

if __name__ == "__main__":
    main()
