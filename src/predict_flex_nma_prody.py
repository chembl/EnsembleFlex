#!/usr/bin/env python

"""
Flexibility prediction with Normal Mode Analysis (NMA) of elastic network models
using mainly the python package ProDy.

Usage:
    python3 predict_flex_nma_prody.py -i <input_directory> -o <output_directory>

Example:
    python3 predict_flex_nma_prody.py -i EnsemblFlex/superimposed -o EnsemblFlex/Prediction_NMA_ProDy
"""

import os
import argparse
import glob
import numpy as np
from prody import *
import matplotlib.pylab as plt

plt.ion()  # The matplotlib.pyplot.ion() function turns on the interactive mode


def parse_arguments():
    parser = argparse.ArgumentParser(description="Perform NMA on a protein structure.")
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


def perform_anm(reference_calphas):
    """Perform ANM analysis on the C-alpha atoms of a protein structure."""
    anm = ANM('')  # Instantiate a ANM instance
    anm.buildHessian(reference_calphas)  # Build Hessian for the reference chain
    anm.calcModes()  # Calculate slowest non-trivial 20 modes
    return anm


def save_anm_results(anm, reference_calphas, reference_structure, output_path):
    """Save results of ANM analysis."""
    # Write modes to file
    writeModes(os.path.join(output_path, 'ANM_modes.txt'), anm)
    print('\nAnisotropic Network Model (ANM) modes are written to ANM_modes.txt.\n')

    # Variance along ANM modes
    with open(os.path.join(output_path, 'ANM_mode_variance.txt'), 'w') as f:
        print('Variance along ANM modes:', file=f)
        for mode in anm[:3]:
            var = calcFractVariance(mode) * 100
            print(f"{mode}  % variance = {var:.2f}", file=f)

    # Collectivity of ANM modes
    with open(os.path.join(output_path, 'ANM_mode_collectivity.txt'), 'w') as f:
        print('ANM mode collectivity:', file=f)
        for mode in anm[:3]:
            coll = calcCollectivity(mode)
            print(f"{mode}  collectivity = {coll:.2f}", file=f)

    # Contact Map
    showContactMap(anm)
    plt.savefig(os.path.join(output_path, 'ContactMap_ANM.png'))
    plt.close()
    print('\nContactMap plot is saved to ContactMap_ANM.png.\n')

    # Cross Correlations
    showCrossCorr(anm[:3])
    plt.savefig(os.path.join(output_path, 'CrossCorrelations_ANM123.png'))
    plt.close()
    print('\nCross Correlation plot is saved to CrossCorrelations_ANM123.png.\n')

    showCrossCorr(anm)
    plt.savefig(os.path.join(output_path, 'CrossCorrelations_ANM.png'))
    plt.close()
    print('\nCross Correlation plot is saved to CrossCorrelations_ANM.png.\n')

    # Square fluctuations
    for i in range(3):
        showSqFlucts(anm[i], hinges=True)
        plt.savefig(os.path.join(output_path, f'Fluctuations_perResidue_ANM{i + 1}.png'))
        plt.close()
        print(
            f'\nGNM square fluctuations per residue plot of ANM{i + 1} is saved to Fluctuations_perResidue_ANM{i + 1}.png.\n')

    # Comparison of MSFs with experimental B-factors (Debye-Waller factor)
    bfactors = reference_calphas.getBetas()
    anm_msfs = calcSqFlucts(anm)
    anm_msfs_rescaled = anm_msfs / np.mean(anm_msfs) * np.mean(bfactors)

    plt.figure(figsize=(9, 5), dpi=300)
    plt.plot(bfactors, 'orange', label='Experimental')
    plt.plot(anm_msfs_rescaled, 'g', lw=1., label='ANM')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'B-factors_vs_ANM-msfs.png'))
    plt.close()

    ## Writing protein structures
    # Extend the ANM model based on Cα-atoms to all heavy atoms by the method extendModel
    anm_aa, atoms_all = extendModel(anm, reference_calphas, reference_structure)
    # write NMD (can be visualized with the VMD plugin NMWizard: VMD Main menu --> Analysis --> Normal Mode Wizard)
    writeNMD(os.path.join(output_path, 'ANM_model_aa.nmd'), anm_aa, atoms_all)
    print("ANM model written to NMD file ANM_model_aa.nmd")

    # Generate trajectories by the method traverseMode
    # This takes steps in both directions starting from the provided structure to generate conformers along the chosen mode.
    # Generate and save trajectories
    for i in range(2):
        anm_traj_aa = traverseMode(anm_aa[i], atoms_all)
        anm_traj_aa.setAtoms(atoms_all)
        writePDB(os.path.join(output_path, f'ANM{i + 1}_traj_aa.pdb'), anm_traj_aa)
        print(f"ANM{i + 1}-allatom trajectory saved as multimodel to ANM{i + 1}_traj_aa.pdb")

    # Save mean square fluctuations computed with sets of ANM modes stored in the B-factor column
    # mean square fluctuations based on all 20 ANM modes, based on the first 10 ANM modes, and based on the first 3 modes
    for num_modes in [20, 10, 3]:
        msf = calcSqFlucts(anm[:num_modes])
        writePDB(os.path.join(output_path, f'ANM_msf_first{num_modes}_on_bfactor.pdb'), reference_calphas, beta=msf)
    print("Mean square fluctuations computed with sets of ANM modes stored in the B-factor column:\n"
          "ANM_msf_first20_on_bfactor.pdb, ANM_msf_first10_on_bfactor.pdb, ANM_msf_first3_on_bfactor.pdb")


def perform_gnm(reference_calphas):
    """Perform GNM analysis on the C-alpha atoms of a protein structure."""
    gnm = GNM('')  # Instantiate a GNM instance
    gnm.buildKirchhoff(reference_calphas)  # Build Kirchhoff matrix for the reference chain
    gnm.calcModes()  # Calculate slowest non-trivial 20 modes
    return gnm


def save_gnm_results(gnm, reference_calphas, reference_structure, output_path):
    """Save results of GNM analysis."""

    # Write modes to file
    writeModes(os.path.join(output_path, 'GNM_modes.txt'), gnm)
    print('\nGaussian Network Model (GNM) modes are written to GNM_modes.txt.\n')

    # Variance along GNM modes
    with open(os.path.join(output_path, 'GNM_mode_variance.txt'), 'w') as f:
        print('Variance along GNM modes:', file=f)
        for mode in gnm[:3]:
            var = calcFractVariance(mode) * 100
            print(f"{mode}  % variance = {var:.2f}", file=f)

    # Collectivity of GNM modes
    with open(os.path.join(output_path, 'GNM_mode_collectivity.txt'), 'w') as f:
        print('GNM mode collectivity:', file=f)
        for mode in gnm[:3]:
            coll = calcCollectivity(mode)
            print(f"{mode}  collectivity = {coll:.2f}", file=f)

    # Contact Map
    showContactMap(gnm)
    plt.savefig(os.path.join(output_path, 'ContactMap_GNM.png'))
    plt.close()
    print('\nContactMap plot is saved to ContactMap_GNM.png.\n')

    # Cross-correlations
    showCrossCorr(gnm)
    plt.savefig(os.path.join(output_path, 'CrossCorrelations_GNM.png'))
    plt.close()
    print('\nCross Correlation plot is saved to CrossCorrelations_GNM.png.\n')

    # Square fluctuations
    for i in range(3):
        showSqFlucts(gnm[i], hinges=True)  # [0], [1], [2]
        plt.savefig(os.path.join(output_path, f'Fluctuations_perResidue_GNM{i + 1}.png'))
        plt.close()
        print(
            f'\nGNM square fluctuations per residue plot of GNM{i + 1} is saved to Fluctuations_perResidue_GNM{i + 1}.png.\n')

    # Comparison of MSFs with experimental B-factors (Debye-Waller factor)
    bfactors = reference_calphas.getBetas()
    gnm_msfs = calcSqFlucts(gnm)
    gnm_msfs_rescaled = gnm_msfs / np.mean(gnm_msfs) * np.mean(bfactors)

    plt.figure(figsize=(9, 5), dpi=300)
    plt.plot(bfactors, 'orange', label='Experimental')
    plt.plot(gnm_msfs_rescaled, 'g', lw=1., label='GNM')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, 'B-factors_vs_GNM-msfs.png'))
    plt.close()

    ## Access and write hinge sites
    # hinges = calcHinges(gnm)
    # node indices in the GNM object
    # calcHinges(gnm[0])  # Hinge sites from the slowest mode
    # calcHinges(gnm[:2])  # Hinge sites identified from multiple modes (e.g. 2 modes)
    # translate to residue numbers corresponding to the hinges
    resnums = reference_calphas.getResnums()
    # mode1_hinges = calcHinges(gnm[0])
    # resnums[mode1_hinges]
    # mode2_hinges = calcHinges(gnm[1])
    # resnums[mode2_hinges]
    # save to file
    with open('GNM_hinge_residues.txt', 'w') as f:
        print('GNM hinge residues', file=f)
        hinge_res = resnums[calcHinges(gnm[:3])]
        print(' - identified from modes 1-3:\n{}'.format(hinge_res), file=f)
        hinge_res = resnums[calcHinges(gnm[0])]
        print(' - identified from mode 1:\n{}'.format(hinge_res), file=f)
        hinge_res = resnums[calcHinges(gnm[1])]
        print(' - identified from mode 2:\n{}'.format(hinge_res), file=f)
        hinge_res = resnums[calcHinges(gnm[2])]
        print(' - identified from mode 3:\n{}'.format(hinge_res), file=f)

    ## Writing protein structures
    # extend the GNM model based on Cα-atoms to all heavy atoms by the method extendModel
    gnm_aa, atoms_all = extendModel(gnm, reference_calphas, reference_structure)
    # write NMD (can be visualized with the VMD plugin NMWizard: VMD Main menu --> Analysis --> Normal Mode Wizard)
    writeNMD('GNM_model_aa', gnm_aa, atoms_all)
    print("GNM model written to NMD file GNM_model_aa.nmd")

    # Save mean square fluctuations computed with sets of GNM modes stored in the B-factor column
    for num_modes in [20, 10, 3]:
        msf = calcSqFlucts(gnm[:num_modes])
        writePDB(os.path.join(output_path, f'GNM_msf_first{num_modes}_on_bfactor.pdb'), reference_calphas, beta=msf)
    print("Mean square fluctuations computed with sets of GNM modes stored in the B-factor column:\n"
          "GNM_msf_first20_on_bfactor.pdb, GNM_msf_first10_on_bfactor.pdb, GNM_msf_first3_on_bfactor.pdb")


def perform_ddd(reference_structure):
    """Perform Dynamical Domain Decomposition of reference_structure (using GNM)."""
    # Dynamical domain decomposition using GNM modes
    ddd_gnm, _ = calcGNM(reference_structure, selstr='all')
    # The number of dynamical domains obtained is approximately equal to the number of modes used (here 8).
    domains = calcGNMDomains(ddd_gnm[:3])
    # assign this data to the AtomGroup:
    reference_structure.setData('domain', domains)
    # using the B-factor field for writing the domains
    writePDB('dynamic_domains_GNM_allatom.pdb', reference_structure, beta=domains)
    print('\nDynamical domains of reference structure are saved in B-factor column of dynamic_domains_GNM_allatom.pdb '
          'for visualisation.\n')


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)

    # List of input PDB files
    pdb_files = glob.glob(os.path.join(input_path, "*.pdb"))

    if not pdb_files:
        raise ValueError("No PDB files found in the input directory.")

    # by default first structure of ensemble is taken as reference structure
    reference_structure = parsePDB(pdb_files[0], compressed=False)  # , subset='ca'
    pdb_id = os.path.splitext(os.path.basename(pdb_files[0]))[0]
    print(f"Reference structure set to {pdb_id}.")

    # select C-alphas of reference structure
    reference_calphas = reference_structure.select('calpha')

    # Perform ANM analysis
    my_anm = perform_anm(reference_calphas)
    save_anm_results(my_anm, reference_calphas, reference_structure, output_path)

    # Perform GNM analysis
    my_gnm = perform_gnm(reference_calphas)
    save_gnm_results(my_gnm, reference_calphas, reference_structure, output_path)

    # Perform Dynamical Domain Decomposition analysis
    perform_ddd(reference_structure)


if __name__ == "__main__":
    main()
