#!/usr/bin/env python3

"""
Renames all pdbXXXX.ent files to XXXX.pdb files
(keeps only the last 4 characters before the extension).

Usage:
    python3 rename_ent_files_to_pdb.py <input_directory>

Example:
    python3 rename_ent_files_to_pdb.py pdbs
"""

import os
import sys
import glob


def rename_ent_files_to_pdb(directory):
    """
    Renames all pdbXXXX.ent files in the specified directory to XXXX.pdb files.

    Args:
        directory (str): The directory containing .ent files to be renamed.
    """
    if not os.path.isdir(directory):
        print(f"Error: The directory '{directory}' does not exist.")
        return

    for filename in glob.iglob(os.path.join(directory, '*.ent')):
        try:
            new_name = os.path.join(directory, filename[-8:-4] + '.pdb')
            os.rename(filename, new_name)
            print(f"Renamed: {filename} to {new_name}")
        except OSError as e:
            print(f"Error renaming {filename}: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 rename_ent_files_to_pdb.py <input_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]
    rename_ent_files_to_pdb(input_directory)
