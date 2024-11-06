#!/usr/bin/env python

"""
Ensemble flexibility analysis through Solvent Accessibility Surface Area (SASA) differences
using sequence alignment for consistent residue numbering.

Usage:
    python3 analyse_flex_sasa_alignment.py -i <input_directory> -o <output_directory>

Example:
    python3 analyse_flex_sasa_alignment.py -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_SASA_Alignment
"""

import os
import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from Bio.PDB import PDBParser, Polypeptide
from Bio.PDB.SASA import ShrakeRupley
from Bio import Align, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices


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


def extract_sequence(structure):
    """Extracts the amino acid sequence from the first chain in a structure."""
    chain = structure[0].child_list[0]  # Assuming single chain per structure
    ppb = Polypeptide.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += str(pp.get_sequence())
    return SeqRecord(seq=sequence, id=structure.id)


def perform_alignment(pdbfiles):
    """Aligns sequences from all PDB files and returns a dictionary of aligned sequences with gaps."""
    parser = PDBParser(PERMISSIVE=True)
    sequences = []

    for fileName in pdbfiles:
        structure_id = os.path.basename(fileName).rsplit('.', 1)[0]
        structure = parser.get_structure(structure_id, fileName)
        seq_record = extract_sequence(structure)
        seq_record.id = structure_id
        sequences.append(seq_record)

    # Perform multiple sequence alignment
    aligner = Align.PairwiseAligner()
    matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = matrix
    alignment = MultipleSeqAlignment(sequences)

    # Convert alignment to a mapping of residue indices with gaps
    aligned_residue_map = {record.id: list(record.seq) for record in alignment}
    return aligned_residue_map


def calculate_sasa(input_path, output_path, aligned_residue_map):
    """Calculate SASA for all PDB files and align results according to sequence alignment."""
    os.chdir(output_path)

    parser = PDBParser(PERMISSIVE=True)
    sr = ShrakeRupley()
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))

    # Initialize the SASA dataframe with aligned residue numbering and gaps
    ref_seq_ids = list(range(len(aligned_residue_map[list(aligned_residue_map.keys())[0]])))
    sasa_df = pd.DataFrame({'Aligned_Index': ref_seq_ids})

    for fileName in pdbfiles:
        structure_id = os.path.basename(fileName).rsplit('.', 1)[0]
        structure = parser.get_structure(structure_id, fileName)

        # Remove waters to avoid interference in SASA calculation
        for chain in structure[0]:
            for residue in list(chain):
                if residue.id[0] == 'W':
                    chain.detach_child(residue.id)

        # Calculate SASA for residues in the structure
        sr.compute(structure[0], level="R")

        # Map calculated SASA to aligned positions
        chain = structure[0].child_list[0]  # Assume single chain
        structure_sasa = []
        for idx, aligned_res in enumerate(aligned_residue_map[structure_id]):
            if aligned_res == '-':  # Gap in alignment
                structure_sasa.append(np.nan)
            else:
                residue_num = list(chain)[idx].id[1] if idx < len(chain) else None
                if residue_num and chain[residue_num].sasa is not None:
                    structure_sasa.append(round(chain[residue_num].sasa, 2))
                else:
                    structure_sasa.append(np.nan)  # No SASA for this residue

        # Add this structure's SASA data to the main DataFrame
        sasa_df[structure_id] = structure_sasa

    return sasa_df


def save_sasa_data(sasa_df):
    """Save SASA data to CSV files and generate a plot."""
    sasa_df.to_csv('SASA_per_structure.csv', index=False)
    print('SASA values per residue and structure are saved to SASA_per_structure.csv.')

    # Calculate global SASA metrics across ensemble
    calc_df = sasa_df.loc[:, ~sasa_df.columns.isin(['Aligned_Index'])]
    sasa_global_df = sasa_df[['Aligned_Index']]
    sasa_global_df['max'] = calc_df.max(axis=1)
    sasa_global_df['min'] = calc_df.min(axis=1)
    sasa_global_df['spread'] = calc_df.max(axis=1) - calc_df.min(axis=1)
    sasa_global_df['sd'] = calc_df.std(axis=1).round(decimals=2)
    sasa_global_df.to_csv('SASA_global.csv', index=False)
    print('SASA global metrics per residue are saved to SASA_global.csv.')

    # Plotting SASA standard deviation
    plt.plot(sasa_global_df['sd'])
    plt.title('Solvent Accessible Surface Area (SASA) differences')
    plt.xlabel('Residue index (aligned)')
    plt.ylabel('SASA Standard Deviation')
    plt.savefig('SASA_sd.png')
    plt.close()
    print('SASA Standard Deviation plot saved to SASA_sd.png.')


def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))
    # Perform alignment
    aligned_residue_map = perform_alignment(pdbfiles)
    # Calculate SASA
    sasa_df = calculate_sasa(input_path, output_path, aligned_residue_map)
    # Save SASA dataframe
    save_sasa_data(sasa_df)


if __name__ == "__main__":
    main()
