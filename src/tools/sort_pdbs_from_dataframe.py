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
    parser = argparse.ArgumentParser(description="Subset PDB files based on a CSV dataframe.")
    parser.add_argument("-i", "--input", dest="input_directory", required=True,
                        help="Input directory path")
    parser.add_argument("-o", "--output", dest="output_directory", required=True,
                        help="Output directory path")
    parser.add_argument("-d", "--dataframe", dest="dataframe_path", required=True,
                        help="Path of the CSV dataframe containing PDB filenames and the specified column")
    parser.add_argument("-c", "--column", dest="column_name", required=True,
                        help="Name of the column in the dataframe to use for subdirectories")

    args = parser.parse_args()

    subset_pdbs_from_dataframe(args.input_directory, args.output_directory, args.dataframe_path, args.column_name)
