#!/usr/bin/env python

"""
Sorts/copies PDB files from the input directory to specific output subdirectories that are created
based on a provided column in a provided dataframe.

Usage:
    python3 sort_pdbs_from_dataframe.py -i <input_directory> -o <output_directory> -d <dataframe.csv> -c <column_name>

Example:
    python3 sort_pdbs_from_dataframe.py -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_Bio3D/Consensus_Clusters
            -d EnsemblFlex/Analysis_Bio3D/cluster_attributions_with_consensus.csv -c Consensus_Cluster
"""

import os
import shutil
import argparse
import pandas as pd


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Subset PDB files based on a CSV dataframe.")
    parser.add_argument("-i", "--input", dest="input_directory", required=True,
                        help="Input directory path")
    parser.add_argument("-o", "--output", dest="output_directory", required=True,
                        help="Output directory path")
    parser.add_argument("-d", "--dataframe", dest="dataframe_path", required=True,
                        help="Path of the CSV dataframe containing PDB filenames and the specified column")
    parser.add_argument("-c", "--column", dest="column_name", required=True,
                        help="Name of the column in the dataframe to use for subdirectories")
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


def validate_dataframe_path(dataframe_path):
    """Validate the JSON file path."""
    if not os.path.isabs(dataframe_path):
        dataframe_path = os.path.abspath(dataframe_path)

    if not os.path.isfile(dataframe_path):
        print("Error: The JSON file does not exist.")
        sys.exit(1)

    return dataframe_path


def subset_pdbs_from_dataframe(input_dir, output_dir, dataframe_path, column_name):
    # Ensure the output directory exists, create it if not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the dataframe from the CSV file
    df = pd.read_csv(dataframe_path)

    # Iterate through the dataframe rows and copy pdb files to the corresponding subdirectories
    for index, row in df.iterrows():
        pdb_basename = row[0]+".pdb"
        cluster_value = row[column_name]
        src_path = os.path.join(input_dir, pdb_basename)
        dest_dir = os.path.join(output_dir, "cluster_"+str(cluster_value))
        dest_path = os.path.join(dest_dir, pdb_basename)

        # Ensure the subdirectory exists, create it if not
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        # Copy the pdb file
        shutil.copy(src_path, dest_path)
        print(f"Copying {pdb_basename} to {dest_dir}")


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_directory, args.output_directory)
    # Validate dataframe file path
    dataframe_path = validate_dataframe_path(args.dataframe_path)
    # Subset pdbs
    subset_pdbs_from_dataframe(input_path, output_path, dataframe_path, args.column_name)
