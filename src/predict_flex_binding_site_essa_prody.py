#!/usr/bin/env python

"""
Essential Site Scanning Analysis (ESSA)
Focused flexibility prediction with Essential Site Scanning Analysis (ESSA) using mainly the python package ProDy.

Usage:
    python3 predict_flex_nma_prody.py -i <input_directory> -o <output_directory>

Example:
    python3 predict_flex_nma_prody.py -i EnsemblFlex/superimposed -o EnsemblFlex/Prediction_NMA_ProDy
"""
# A ProDy tutorial on ESSA is provided here: http://prody.csb.pitt.edu/_static/ipynb/workshop2021/prody_essa.html

import os
import argparse
import glob
import numpy as np
from prody import *
import matplotlib.pylab as plt

plt.ion()  # The matplotlib.pyplot.ion() function turns on the interactive mode


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform Essential Site Scanning Analysis (ESSA) on a protein structure.")
    parser.add_argument("-i", "--input", dest="input_path", required=True,
                        help="input dataset directory path", metavar="PATH")
    parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                        help="output directory name", metavar="STRING")
    parser.add_argument("-b", "--bsresidues", dest="bsresidues", default=None,
                        help="binding site residues as file", metavar="STRING")
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


def validate_binding_site_file(bsresidues):
    """Validate binding site residues file."""
    if bsresidues and os.path.isfile(bsresidues):
        return os.path.abspath(bsresidues)
    else:
        print("No valid binding site residue list provided.")
        return None

def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)
    # Validate residue file
    bsresiduefile = validate_binding_site_file(args.bsresidues)

    # Change working directory to output_path
    os.chdir(output_path)

    # List of input PDB files
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))

    # Parsing only reference structure (by default, the first structure of the ensemble is taken as the reference structure)
    reference_structure = parsePDB(pdbfiles[0], compressed=False, title='reference_structure')

    # Select C-alphas of reference structure
    reference_calphas = reference_structure.select('calpha')

    if bsresiduefile:
        with open(bsresiduefile) as f:
            residues = f.read().strip()
    else:
        print("No binding site residues provided. Proceeding without highlighting binding site residues.")
        residues = ""

    ### Essential Site Scanning Analysis (ESSA)

    # Instantiate an ESSA object
    essa = ESSA()
    # Set system
    essa.setSystem(reference_structure)
    # Perform scanning
    essa.scanResidues()

    # Plot Z-Scores per residue with highlighted binding site residues
    # ProDy selection string to highlight binding residues
    if residues:
        prody_sele = 'resnum ' + residues.replace("\n", " ")
        print("\nResidue numbers used from file are: ")
        print(prody_sele)
    else:
        prody_sele = None

    # The blue dashed baseline shows the q-th quantile of the profile, which is by default q=0.75, representing the top quartile.
    with plt.style.context({'figure.figsize': (9, 7), 'figure.dpi': 100}):
        essa.showESSAProfile(rescode=True, sel=prody_sele)
    plt.savefig('ESSA_profile_of_reference_structure.png')
    plt.close()
    print('\nESSA profile plot is saved to ESSA_profile_of_reference_structure.png.')

    # Save ESSA Z-scores
    essa.saveESSAZscores()
    print("\nESSA Z-scores written to file")

    # Save PDB file with Z-scores in the B-factor column
    essa.writeESSAZscoresToPDB()
    print("\nPDB file with ESSA Z-scores in B-factor column written to file\n")


if __name__ == "__main__":
    main()
