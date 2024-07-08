#!/usr/bin/env python

"""
Ensemble structural superpositioning using the Python package ProDy.

Usage:
    python3 superimpose_prody.py -i <input_directory> -o <output_directory>

Example:
    python3 superimpose_prody.py -i pdbs -o EnsemblFlex
"""

import glob
import os
import argparse
import subprocess
from prody import *


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform global iterative superpositioning of protein structures.")
    parser.add_argument("-i", "--input", dest="input_path", required=True,
                        help="input dataset directory path", metavar="PATH")
    parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                        help="output directory name", metavar="STRING")
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


def superimpose_structures(input_path, output_path):
    """Superimpose structures and save results."""
    toolsdir = os.path.dirname(os.path.realpath(__file__)) + "/tools"

    # Change working directory to output_path
    os.chdir(output_path)

    # Create 'superimposed' directory if it does not exist
    if not os.path.exists("superimposed"):
        os.mkdir("superimposed")
    os.chdir("superimposed")

    # List of input PDB files
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))

    # Parsing structures
    structures = parsePDB(pdbfiles, compressed=False)

    # Build ensemble
    ensemble = buildPDBEnsemble(structures)
    # Method explained: buildPDBEnsemble() maps each structure against the reference structure using a function such
    # as mapOntoChain(). The reference structure is automatically the first member of list provided.

    # Perform iterative superimposition
    ensemble.iterpose()
    # Method explained: Iteratively superpose the ensemble until convergence. Initially, all conformations are aligned
    # with the reference coordinates. Then mean coordinates are calculated, and are set as the new reference
    # coordinates. This is repeated until reference coordinates do not change. This is determined by the value of RMSD
    # between the new and old reference coordinates. Note that at the end of the iterative procedure the reference
    # coordinate set will be average of conformations in the ensemble.

    # Save coordinates
    writePDB("ensemble.pdb", ensemble)
    print("Superimposed ensemble coordinates are saved to ensemble.pdb.\n")

    # Split the ensemble file into separate models
    split_ensemble_dir = os.path.join(output_path, "superimposed", "split_ensemble")
    if not os.path.exists(split_ensemble_dir):
        os.mkdir(split_ensemble_dir)
    os.chdir(split_ensemble_dir)

    # print(f"toolsdir = {toolsdir}")
    subprocess.call(['python3', os.path.join(toolsdir, 'pdb_splitmodel.py'),
                     os.path.join(output_path, 'superimposed', 'ensemble.pdb')])

    print("\nMultimodelfile split and saved to directory /superimposed/split_ensemble.\n")


def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)
    # Superimpose and save structures
    superimpose_structures(input_path, output_path)


if __name__ == "__main__":
    main()
