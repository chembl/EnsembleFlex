#!/usr/bin/env python

"""
Sorts/copies PDB files from the input directory into two directories "structures_with_gaps" and
"structures_without_gaps" based on occurrence of missing residues in the pdb sequence, as indicated by
residue numbering.
(The sequence is filled in with ‘X’ characters to match the size of the missing region in the output txt
file "structures_with_gaps.txt", which can be used for validation.)

Usage:
    python sort_pdbs_has_gap_in_pdb_seq.py -i <input_directory> -o <output_directory>

Example:
    python sort_pdbs_has_gap_in_pdb_seq.py -i EnsemblFlex/superimposed -o EnsemblFlex/

"""

import os
import argparse
from Bio import SeqIO
from shutil import copyfile, rmtree


def subset_pdbs_on_gap_in_seq(input_directory, output_directory):
    # Remove existing subdirectories if they exist
    output_folder_with_gaps = os.path.join(output_directory, "structures_with_gaps")
    output_folder_without_gaps = os.path.join(output_directory, "structures_without_gaps")

    for subdirectory in [output_folder_with_gaps, output_folder_without_gaps]:
        if os.path.exists(subdirectory):
            rmtree(subdirectory)

    # Create new output subdirectories
    os.makedirs(output_folder_with_gaps)
    os.makedirs(output_folder_without_gaps)

    list_files_with_gaps = []
    list_files_without_gaps = []

    # Iterate through PDB files in the input directory
    for pdb_file in os.listdir(input_directory):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(input_directory, pdb_file)

            # Extract sequence using SeqIO
            # "pdb-atom" calls the Bio.SeqIO.PdbIO.PdbAtomIterator (https://biopython.org/docs/1.75/api/Bio.SeqIO.PdbIO.html):
            # Where amino acids are missing from the structure, as indicated by residue numbering,
            # the sequence is filled in with ‘X’ characters to match the size of the missing region
            seq_from_coords = next(SeqIO.parse(pdb_path, "pdb-atom"))
            seq = seq_from_coords.seq

            # Sort pdb file into subdirectory
            if 'X' in str(seq):
                list_files_with_gaps.append(">"+pdb_file)
                list_files_with_gaps.append(seq)
                # Copy PDB file to the "structures_with_gaps" folder
                copyfile(pdb_path, os.path.join(output_folder_with_gaps, pdb_file))
            else:
                list_files_without_gaps.append(">"+pdb_file)
                list_files_without_gaps.append(seq)
                # Copy PDB file to the "structures_without_gaps" folder
                copyfile(pdb_path, os.path.join(output_folder_without_gaps, pdb_file))

    output_files_with_gaps = os.path.join(output_directory, "structures_with_gaps.txt")
    with open(output_files_with_gaps, "w") as f:
        for struct in list_files_with_gaps:
            f.write(f"{struct}\n")
    print("Structures containing gaps (missing residues) are copied to "+output_files_with_gaps)

    output_files_without_gaps = os.path.join(output_directory, "structures_without_gaps.txt")
    with open(output_files_without_gaps, "w") as f:
        for struct in list_files_without_gaps:
            f.write(f"{struct}\n")
    print("Structures containing a continuous sequence (no missing residues) are copied to "+output_files_without_gaps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Subset PDB files based on gaps in sequence (taken from atom records).")
    parser.add_argument("-i", "--input", dest="input_directory", required=True,
                        help="Input directory path containing PDB files")
    parser.add_argument("-o", "--output", dest="output_directory", required=True,
                        help="Output directory path")
    args = parser.parse_args()

    subset_pdbs_on_gap_in_seq(args.input_directory, args.output_directory)
