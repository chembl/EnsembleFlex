#!/usr/bin/env python

"""
Ensemble flexibility analysis using mainly the Python package ProDy.

Usage:
    python3 analyse_flex_prody.py -i <input_directory> -o <output_directory>

Example:
    python3 analyse_flex_prody.py -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_ProDy
"""

import os
import sys
import argparse
import glob
from prody import *
import matplotlib.pylab as plt
import numpy as np
from datetime import datetime


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Perform flexibility analysis on protein structures.")
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


def setup_logging(output_path):
    """Setup logging to file and console."""
    class Logger(object):
        def __init__(self):
            self.terminal = sys.stdout
            self.log = open(os.path.join(output_path, "logfile.log"), "a")

        def write(self, message):
            self.terminal.write(message)
            self.log.write(f"{datetime.now()}: {message}")

        def flush(self):
            pass

    sys.stdout = Logger()


def perform_rmsf_analysis(ensemble, reference_structure):
    """Perform RMSF analysis and save results."""
    rmsf = ensemble.getRMSFs()
    plt.plot(rmsf)
    plt.xlabel('Residue index')
    plt.ylabel('RMSF')
    plt.savefig('RMSF.png')
    plt.close()
    print('Root mean square fluctuation (RMSF) plot saved to RMSF.png.')

    writePDB('RMSFonReference.pdb', reference_structure, beta=rmsf)
    print('RMSF values are saved in B-factor column of RMSFonReference.pdb for visualization.')


def perform_pca_analysis(ensemble, pdbfiles):
    """Perform PCA analysis and save results."""
    pca = PCA('')
    pca.buildCovariance(ensemble)
    pca.calcModes()

    if len(pdbfiles) > 3:  # if more than 3 structures are present
        # Perform Singular Value Decomposition (SVD) for faster and more memory efficient calculation of principal modes
        pca_svd = PCA('ensemblePCAsvd')
        pca_svd.performSVD(ensemble)

        abs(pca_svd.getEigvals()[:20] - pca.getEigvals()).max()
        abs(calcOverlap(pca, pca_svd).diagonal()[:20]).min()

        # Plot overlap table of PCA and PCA-SVD
        showOverlapTable(pca, pca_svd[:20])
        plt.savefig('Overlap_PCA_PCA-SVD.png')
        plt.close()
        print('Overlap plot of PCA and PCA with singular value decomposition saved to Overlap_PCA_PCA-SVD.png.')
    else:
        print("Ensembles must contain more than 3 structures to perform singular value decomposition.")

    # Print variance along PCs to console and save to file
    print('Variance along PCs:')
    for mode in pca[:3]:
        var = calcFractVariance(mode) * 100
        print(f'{mode}  % variance = {var:.2f}')
    with open('PCA_variance.txt', 'w') as f:
        print('Variance along PCs:', file=f)
        for mode in pca[:3]:
            var = calcFractVariance(mode) * 100
            print(f'{mode}  % variance = {var:.2f}', file=f)

    # Print PCA mode collectivity to console and save to file
    print('PCA mode collectivity:')
    for mode in pca[:3]:
        coll = calcCollectivity(mode)
        print(f'{mode}  collectivity = {coll:.2f}')
    with open('PCA_mode_collectivity.txt', 'w') as f:
        print('PCA mode collectivity:', file=f)
        for mode in pca[:3]:
            coll = calcCollectivity(mode)
            print(f'{mode}  collectivity = {coll:.2f}', file=f)

    # Save PCA eigenvectors to file
    writeArray('PCA_eigvecs.txt', pca.getEigvecs())
    print('PCA eigenvectors are written to PCA_eigvecs.txt.')

    # Plot PCA square fluctuations and save to file
    showSqFlucts(pca[:3])
    plt.savefig('Fluctuations_PCA.png')
    plt.close()
    print('PCA square fluctuations plot is saved to Fluctuations_PCA.png.')

    # Plot cross correlations of PC1 per residue and save to file
    showCrossCorr(pca[0])
    plt.savefig('CrossCorrelations_PCA1perResidue.png')
    plt.close()
    print('Cross Correlations of PC1 per residue is saved to CrossCorrelations_PCA1perResidue.png.')


def perform_eda_analysis(ensemble):
    """Perform EDA analysis and save results."""
    eda = EDA('')
    eda.buildCovariance(ensemble)
    eda.calcModes()

    # Print variance along modes to console and save to file
    print('Variance along modes:')
    for mode in eda[:3]:
        var = calcFractVariance(mode) * 100
        print(f'{mode}  % variance = {var:.2f}')
    with open('EDA_variance.txt', 'w') as f:
        print('Variance along modes:', file=f)
        for mode in eda[:3]:
            var = calcFractVariance(mode) * 100
            print(f'{mode}  % variance = {var:.2f}', file=f)

    # Plot EDA square fluctuations and save to file
    showSqFlucts(eda[:3])
    plt.savefig('Fluctuations_EDA.png')
    plt.close()
    print('EDA square fluctuations plot is saved to Fluctuations_EDA.png.')


def main():
    # Parse command-line arguments
    args = parse_arguments()
    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)
    # Set up logging to file and console
    setup_logging(output_path)

    # Change working directory to output_path
    os.chdir(output_path)

    # List of input PDB files
    pdbfiles = glob.glob(os.path.join(input_path, "*.pdb"))

    # Parse structures and select C-alphas (ca) of the first structure as reference
    structures = parsePDB(pdbfiles, subset='ca', compressed=False)
    reference_structure = structures[0] # by default first structure of ensemble is taken as reference structure and
    # For unresolved atoms, the coordinates of the reference structure is assumed in RMSD calculations and superpositions.
    # select only C-alphas of reference structure (not explicitly needed here)
    #reference_calphas = reference_structure.select('calpha')

    # Build PDB ensemble without superposition
    # buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
    # The reference structure is automatically the first member of list provided.
    # Here superpositioning is omitted by setting it to 'False'.
    ensemble = buildPDBEnsemble(structures, superpose=False)
    # set c-alpha reference (not explicitly needed here)
    # ensemble.setAtoms(reference_structure.calpha)

    ### 1. RMSF Analysis
    # Perform RMSF analysis
    perform_rmsf_analysis(ensemble, reference_structure)

    ### 2. Principal Component Analysis (PCA)
    # PCA is a method that identifies the components which account for the greatest amount of variability in the ensemble.
    # Perform PCA analysis
    perform_pca_analysis(ensemble, pdbfiles)

    ### 3. Essential Dynamics Analysis (EDA)
    # Perform EDA analysis
    perform_eda_analysis(ensemble)

    # ### Generate Report # Not in use
    # import analyse_flex_prody_reporting as reporting
    # reporting.generate_report(str_input_path=args.input_path, output_path=output_path)


if __name__ == "__main__":
    main()
