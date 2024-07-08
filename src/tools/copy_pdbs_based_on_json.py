#!/usr/bin/env python

"""
Sorts files based on keys in a JSON file.
Example JSON file: http://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/P/P24941/P24941.json

Usage:
    python3 copy_pdbs_based_on_json.py -j <json_file> -i <input_directory> -o <output_directory>

Example:
    python3 copy_pdbs_based_on_json.py -j /mnt/data/P24941.json -i /path/to/source/folder -o /path/to/destination/folder
"""

import os
import shutil
import json
import argparse
import sys


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Copy files based on JSON keys.")
    parser.add_argument("-j", "--json", dest="json_file", required=True,
                        help="Path to the JSON file", metavar="FILE")
    parser.add_argument("-i", "--input", dest="input_path", required=True,
                        help="Input dataset directory path", metavar="PATH")
    parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                        help="Output directory name", metavar="STRING")
    return parser.parse_args()


def validate_directories(input_path, output_path):
    """Validate input and output directories."""
    # Convert input_path to an absolute path if it is not already
    if not os.path.isabs(input_path):
        input_path = os.path.abspath(input_path)

    if not os.path.isdir(input_path):
        print("Error: The input directory does not exist.")
        sys.exit(1)

    # Convert output_path to an absolute path if it is not already
    if not os.path.isabs(output_path):
        output_path = os.path.abspath(os.path.join(os.getcwd(), output_path))

    if not os.path.isdir(output_path):
        try:
            os.makedirs(output_path)
        except OSError as e:
            print(f"Error: Could not create output directory '{output_path}'. {e}")
            sys.exit(1)

    return input_path, output_path

def validate_json_file(json_file):
    """Validate the JSON file path."""
    if not os.path.isabs(json_file):
        json_file = os.path.abspath(json_file)

    if not os.path.isfile(json_file):
        print("Error: The JSON file does not exist.")
        sys.exit(1)

    return json_file


def copy_files_based_on_json(json_file, src_folder, dest_folder):
    """Copy files based on keys in the JSON file."""
    # Load JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Ensure the destination folder exists
    os.makedirs(dest_folder, exist_ok=True)

    # Iterate over the keys in the JSON data
    for key in data.keys():
        file_name = f"{key}.pdb"
        src_file_path = os.path.join(src_folder, file_name)
        dest_file_path = os.path.join(dest_folder, file_name)

        # Check if the file exists in the source folder
        if os.path.isfile(src_file_path):
            # Copy the file to the destination folder
            shutil.copy2(src_file_path, dest_file_path)
            print(f"Copied: {file_name}")
        else:
            print(f"File not found: {file_name}")


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)

    # Validate JSON file path
    json_file = validate_json_file(args.json_file)

    # Change working directory to output_path
    os.chdir(output_path)

    # Copy files based on JSON keys
    copy_files_based_on_json(json_file, input_path, output_path)


if __name__ == "__main__":
    main()
