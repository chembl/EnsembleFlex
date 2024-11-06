#!/usr/bin/env python

"""
Ensemble flexibility analysis through Solvent Accessibility Surface Area (SASA) differences
using the Python package Biopython.

Usage:
    python3 analyse_flex_sasa_biopython.py -i <input_directory> -o <output_directory>

Example:
    python3 analyse_flex_sasa_biopython.py -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_SASA_Biopython
"""

import os
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import Polypeptide


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform SASA calculations on protein structures.")
    parser.add_argument("-i", "--input", dest="input_path", required=True,
                        help="input dataset directory path", metavar="PATH")
    parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                        help="output directory name", metavar="STRING")
    return parser.parse_args()


def validate_directories(input_path, output_path):
    """Validate input and output directories."""
    if not os.path.isabs(input_path):
        input_path = os.path.abspath(input_path)

    if not os.path.isdir(input_path):
        print("Error: The input directory does not exist.")
        sys.exit(1)

    if not os.path.isabs(output_path):
        output_path = os.path.abspath(os.path.join(os.getcwd(), output_path))

    if not os.path.isdir(output_path):
        try:
            os.makedirs(output_path)
        except OSError as e:
            print(f"Error: Could not create output directory '{output_path}'. {e}")
            sys.exit(1)

    return input_path, output_path


def extract_residue_info(chain):
    """Extract residue numbers and names from a chain, excluding non-standard amino acids."""
    residues = [(res.id[1], res.get_resname()) for res in chain if Polypeptide.is_aa(res, standard=True)]
    return residues


def calculate_sasa(input_path, output_path):
    """Calculate SASA for all PDB files in the input directory using the first structure as a reference."""
    os.chdir(output_path)

    parser = PDBParser(PERMISSIVE=True)
    sr = ShrakeRupley()
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))

    if not pdbfiles:
        print("Error: No PDB files found in the input directory.")
        sys.exit(1)

    # Process the first structure to use as a reference for residue numbering
    reference_structure_id = os.path.basename(pdbfiles[0]).rsplit('.', 1)[0]
    reference_structure = parser.get_structure(reference_structure_id, pdbfiles[0])
    ref_chain = reference_structure[0].child_list[0]  # Assume single chain
    ref_residues = extract_residue_info(ref_chain)

    # Initialize SASA dataframe with reference residues
    sasa_df = pd.DataFrame(ref_residues, columns=['ResNum', 'ResName'])

    # Process each PDB file and calculate SASA
    for fileName in pdbfiles:
        structure_id = os.path.basename(fileName).rsplit('.', 1)[0]
        structure = parser.get_structure(structure_id, fileName)

        # Remove waters from structure to avoid SASA interference
        for chain in structure[0]:
            for residue in list(chain):
                if residue.id[0] == 'W':
                    chain.detach_child(residue.id)

        # Calculate SASA
        sr.compute(structure[0], level="R")

        # Collect SASA values for residues, mapping them to reference residue numbers
        my_list = []
        chain = structure[0].child_list[0]  # Assume single chain
        for ref_id, ref_name in ref_residues:
            matching_res = chain[ref_id] if ref_id in chain else None
            if matching_res:
                sasa_value = round(matching_res.sasa, 2)
            else:
                sasa_value = np.nan  # Fill missing residues with NaN
            my_list.append((ref_id, ref_name, sasa_value))

        # Create a DataFrame for the current structure and merge with main SASA dataframe
        my_df = pd.DataFrame(my_list, columns=['ResNum', 'ResName', structure_id])
        sasa_df = pd.merge(sasa_df, my_df, on=['ResNum', 'ResName'], how='outer')

    return sasa_df


def save_sasa_data(sasa_df):
    """Save SASA data to CSV files and generate a plot."""
    sasa_df.to_csv('SASA_per_structure.csv', index=False)
    print('SASA values per residue and structure are saved to SASA_per_structure.csv.')

    # Calculate global SASA metrics across ensemble
    calc_df = sasa_df.loc[:, ~sasa_df.columns.isin(['ResNum', 'ResName'])]
    sasa_global_df = sasa_df[['ResNum', 'ResName']]
    sasa_global_df['max'] = calc_df.max(axis=1)
    sasa_global_df['min'] = calc_df.min(axis=1)
    sasa_global_df['spread'] = calc_df.max(axis=1) - calc_df.min(axis=1)
    sasa_global_df['sd'] = calc_df.std(axis=1).round(decimals=2)
    sasa_global_df.to_csv('SASA_global.csv', index=False)
    print('SASA global metrics per residue are saved to SASA_global.csv.')

    # Plotting SASA standard deviation
    plt.plot(sasa_global_df['sd'])
    plt.title('Solvent Accessible Surface Area (SASA) differences')
    plt.xlabel('Residue index')
    plt.ylabel('SASA Standard Deviation')
    plt.savefig('SASA_sd.png')
    plt.close()
    print('SASA Standard Deviation plot saved to SASA_sd.png.')


def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)
    # Calculate SASA
    sasa_df = calculate_sasa(input_path, output_path)
    # Save SASA dataframe
    save_sasa_data(sasa_df)


if __name__ == "__main__":
    main()
