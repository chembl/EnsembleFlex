#!/usr/bin/env python

import os
import argparse
import glob
import numpy as np
from prody import *
import matplotlib.pylab as plt

plt.ion()  # The matplotlib.pyplot.ion() function turns on the interactive mode


#------- File Argument/Option Parser -------#
parser = argparse.ArgumentParser(description="Perform NMA on a protein structure.")
parser.add_argument("-i", "--input", dest="input_path", required=True,
                    help="input dataset directory path", metavar="PATH")
parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                    help="output directory name", metavar="STRING")

args = parser.parse_args()
if os.path.isdir(args.input_path) and os.path.isdir(args.output_dir):
    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_dir)
elif os.path.isdir(args.input_path) and not os.path.isdir(args.output_dir):
    input_path = os.path.abspath(args.input_path)
    if os.path.isdir(os.path.dirname(args.output_dir)):
        output_path = args.output_dir
        os.mkdir(output_path)
    else:
        output_path = args.input_path+"/"+args.output_dir
    if not os.path.exists(output_path):
        os.mkdir(output_path)
else:
    print("Please specify a valid input directory path to your structure files.\n")
#-------------------------------------------#

# change working directory to output_path
os.chdir(output_path)


# list of input pdb files
pdbfiles = glob.glob(input_path+"/*.pdb")
# parsing only reference structure
# by default first structure of ensemble is taken as reference structure
reference_structure = parsePDB(pdbfiles[0], compressed=False)  #, subset='ca'
# select C-alphas of reference structure
reference_calphas = reference_structure.select('calpha')




### Anisotropic Network Model (ANM) Normal Mode Analysis (NMA)
# The ANM allows for the identification of the most impactful (slowest) modes in dynamics of a single protein,
# which we can compare to the principal components from PCA.
# ANM calculations

anm = ANM('')                 # Instantiate a ANM instance
anm.buildHessian(reference_calphas)   # Build Hessian for the reference chain
anm.calcModes()                   # Calculate slowest non-trivial 20 modes

# write modes to file
writeModes('ANM_modes.txt', anm) # This function is based on writeArray
print('\nAnisotropic Network Model (ANM) modes are written to ANM_modes.txt.\n')


## Variance along ANM modes
# print to sys.stdout
print('Variance along ANM modes:')
for mode in anm[:3]:
    var = calcFractVariance(mode)*100
    print(str(mode) + '  % variance = {:.2f}'.format(var))
    #print('{0:s}  % variance = {1:.2f}'.format(mode, var))
# save to file
with open('ANM_mode_variance.txt', 'w') as f:
    print('Variance along ANM modes:', file=f)
    for mode in anm[:3]:
        var = calcFractVariance(mode) * 100
        print(str(mode) + '  % variance = {:.2f}'.format(var), file=f)

## Collectivity of ANM modes
# print to sys.stdout
print('ANM mode collectivity:')
for mode in anm[:3]:    # Print ANM mode collectivity
    coll = calcCollectivity(mode)
    print(str(mode) + '  collectivity = {:.2f}'.format(coll))
    #print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))
# save to file
with open('ANM_mode_collectivity.txt', 'w') as f:
    print('ANM mode collectivity:', file=f)
    for mode in anm[:3]:  # Print ANM mode collectivity
        coll = calcCollectivity(mode)
        print(str(mode) + '  collectivity = {:.2f}'.format(coll), file=f)


## Contact Map
showContactMap(anm)
plt.savefig('ContactMap_ANM.png')
plt.close()
print('\nContactMap plot is saved to ContactMap_ANM.png.\n')


## Cross Correlations
# We can also see how correlated the motions for each residue are with each other residue.
# We see similar patterns for the two methods (PCA and ANM), especially when using a large number of modes.
showCrossCorr(anm[:3])#;
plt.savefig('CrossCorrelations_ANM123.png')
plt.close()
print('\nCross Correlation plot is saved to CrossCorrelations_ANM123.png.\n')
showCrossCorr(anm)#;
plt.savefig('CrossCorrelations_ANM.png')
plt.close()
print('\nCross Correlation plot is saved to CrossCorrelations_ANM.png.\n')


## Square fluctuations
showSqFlucts(anm[0], hinges=True)#;
plt.savefig('Fluctuations_perResidue_ANM1.png')
plt.close()
print('\nGNM square fluctuations per residue plot of ANM1 is saved to Fluctuations_perResidue_ANM1.png.\n')

showSqFlucts(anm[1], hinges=True)#;
plt.savefig('Fluctuations_perResidue_ANM2.png')
plt.close()
print('\nGNM square fluctuations per residue plot of ANM2 is saved to Fluctuations_perResidue_ANM2.png.\n')

showSqFlucts(anm[2], hinges=True)#;
plt.savefig('Fluctuations_perResidue_ANM3.png')
plt.close()
print('\nGNM square fluctuations per residue plot of ANM3 is saved to Fluctuations_perResidue_ANM3.png.\n')

## Comparison of MSFs with experimental B-factors (Debye-Waller factor)
# rescale MSFs
bfactors = reference_calphas.getBetas()
anm_msfs = calcSqFlucts(anm)
anm_msfs_rescaled = anm_msfs / np.mean(anm_msfs) * np.mean(bfactors)

plt.figure(figsize=(9, 5), dpi=300)
plt.plot(bfactors, 'orange', label='Experimental')
plt.plot(anm_msfs_rescaled, 'g', lw=1., label='ANM')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('B-factors_vs_ANM-msfs.png')
plt.close()

# ## Protein structure bipartition # not working for ANM modes!
# showProtein(reference_calphas, mode=anm[0])
# plt.savefig('ProteinStructureBipartition_ANM1.png')
# plt.close()
# print('\nProtein structure bipartition plot based on ANM1 is saved to ProteinStructureBipartition_ANM1.png.\n')


## Writing protein structures
# extend the ANM model based on Cα-atoms to all heavy atoms by the method extendModel
anm_aa, atoms_all = extendModel(anm, reference_calphas, reference_structure)
# write NMD (can be visualized with the VMD plugin NMWizard: VMD Main menu --> Analysis --> Normal Mode Wizard)
writeNMD('ANM_model_aa', anm_aa, atoms_all)
print("ANM model written to NMD file ANM_model_aa.nmd")
# generate a trajectory by the method traverseMode
# This takes steps in both directions starting from the provided structure to generate conformers along the chosen mode.
anm1_traj_ca = traverseMode(anm[0], reference_calphas)  #, rmsd=1.5
anm1_traj_ca.setAtoms(reference_calphas)
# write PDB
writePDB('ANM1_traj_ca.pdb', anm1_traj_ca)
print("ANM1-Calpha trajectory saved as multimodel to ANM1_traj_ca.pdb")
# generate trajectory for all-atom ANM1
anm1_traj_aa = traverseMode(anm_aa[0], atoms_all)
anm1_traj_aa.setAtoms(atoms_all)
# write PDB
writePDB('ANM1_traj_aa.pdb', anm1_traj_aa)
print("ANM1-allatom trajectory saved as multimodel to ANM1_traj_aa.pdb")
# generate trajectory for all-atom ANM2
anm2_traj_aa = traverseMode(anm_aa[1], atoms_all)
anm2_traj_aa.setAtoms(atoms_all)
# write PDB
writePDB('ANM2_traj_aa.pdb', anm2_traj_aa)
print("ANM2-allatom trajectory saved as multimodel to ANM2_traj_aa.pdb")

# Mean square fluctuations computed with sets of ANM modes stored in the B-factor column
msf_all20 = calcSqFlucts(anm)  # mean square fluctuations based on all 20 ANM modes
writePDB('ANM_msf_all20_on_bfactor.pdb', reference_calphas, beta=msf_all20)
msf_first10 = calcSqFlucts(anm[:10])  # mean square fluctuations based on the first 10 ANM modes
writePDB('ANM_msf_first10_on_bfactor.pdb', reference_calphas, beta=msf_first10)
msf_first3 = calcSqFlucts(anm[:3])  # mean square fluctuations based on the first 3 ANM modes
writePDB('ANM_msf_first3_on_bfactor.pdb', reference_calphas, beta=msf_first3)
print("Mean square fluctuations computed with sets of ANM modes stored in the B-factor column:\n"
      "ANM_msf_all20_on_bfactor.pdb, ANM_msf_first10_on_bfactor.pdb, ANM_msf_first3_on_bfactor.pdb")



### Gaussian Network Model (GNM) Normal Mode Analysis (NMA)

# GNM calculations
gnm = GNM('')  # Instantiate a GNM instance
gnm.buildKirchhoff(reference_calphas)  # default parameters: buildKirchhoff(calphas, cutoff=10.0, gamma=1.0)
#gnm.getKirchhoff()  # get a copy of the Kirchhoff matrix
# calculate normal modes from the Kirchhoff matrix
gnm.calcModes() # by default 20 non-zero (or non-trivial) modes are calculated
# Normal mode indices start from 0, so slowest mode has index 0.

# write modes to file
writeModes('GNM_modes.txt', gnm) # This function is based on writeArray
print('\nGaussian Network Model (GNM) modes are written to GNM_modes.txt.\n')


## Variance along GNM modes
# save to file
with open('GNM_mode_variance.txt', 'w') as f:
    print('Variance along GNM modes:', file=f)
    for mode in gnm[:3]:
        var = calcFractVariance(mode) * 100
        print(str(mode) + '  % variance = {:.2f}'.format(var), file=f)

## Collectivity of GNM modes
# save to file
with open('GNM_mode_collectivity.txt', 'w') as f:
    print('GNM mode collectivity:', file=f)
    for mode in gnm[:3]:  # Print GNM mode collectivity
        coll = calcCollectivity(mode)
        print(str(mode) + '  collectivity = {:.2f}'.format(coll), file=f)


# Access hinge sites
#hinges = calcHinges(gnm)
# node indices in the GNM object
#calcHinges(gnm[0])  # Hinge sites from the slowest mode
#calcHinges(gnm[:2])  # Hinge sites identified from multiple modes (e.g. 2 modes)
# translate to residue numbers corresponding to the hinges
resnums = reference_calphas.getResnums()
#mode1_hinges = calcHinges(gnm[0])
#resnums[mode1_hinges]
#mode2_hinges = calcHinges(gnm[1])
#resnums[mode2_hinges]
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


## Plotting results

## Contact Map
showContactMap(gnm)
plt.savefig('ContactMap_GNM.png')
plt.close()
print('\nContactMap plot is saved to ContactMap_GNM.png.\n')


## Cross-correlations
showCrossCorr(gnm)
plt.savefig('CrossCorrelations_GNM.png')
plt.close()
print('\nCross Correlation plot is saved to CrossCorrelations_GNM.png.\n')


## Square fluctuations per residue (per mode)
showSqFlucts(gnm[0], hinges=True)
plt.savefig('Fluctuations_perResidue_GNM1.png')
plt.close()
print('\nGNM square fluctuations per residue plot of GNM1 is saved to Fluctuations_perResidue_GNM1.png.\n')

showSqFlucts(gnm[1], hinges=True)
plt.savefig('Fluctuations_perResidue_GNM2.png')
plt.close()
print('\nGNM square fluctuations per residue plot of GNM2 is saved to Fluctuations_perResidue_GNM2.png.\n')

showSqFlucts(gnm[2], hinges=True)
plt.savefig('Fluctuations_perResidue_GNM3.png')
plt.close()
print('\nGNM square fluctuations per residue plot of GNM3 is saved to Fluctuations_perResidue_GNM3.png.\n')


## Comparison of MSFs with experimental B-factors (Debye-Waller factor)
# rescale MSFs
bfactors = reference_calphas.getBetas()
gnm_msfs = calcSqFlucts(gnm)
gnm_msfs_rescaled = gnm_msfs / np.mean(gnm_msfs) * np.mean(bfactors)

plt.figure(figsize=(9, 5), dpi=300)
plt.plot(bfactors, 'orange', label='Experimental')
plt.plot(gnm_msfs_rescaled, 'g', lw=1., label='GNM')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('B-factors_vs_GNM-msfs.png')
plt.close()


## Protein structure bipartition
showProtein(reference_calphas, mode=gnm[0])
plt.savefig('ProteinStructureBipartition_GNM1.png')
plt.close()
print('\nProtein structure bipartition plot based on GNM1 is saved to ProteinStructureBipartition_GNM1.png.\n')


## Writing protein structures
# extend the GNM model based on Cα-atoms to all heavy atoms by the method extendModel
gnm_aa, atoms_all = extendModel(gnm, reference_calphas, reference_structure)
# write NMD (can be visualized with the VMD plugin NMWizard: VMD Main menu --> Analysis --> Normal Mode Wizard)
writeNMD('GNM_model_aa', gnm_aa, atoms_all)
print("GNM model written to NMD file GNM_model_aa.nmd")

# Mean square fluctuations computed with sets of GNM modes stored in the B-factor column
msf_all20 = calcSqFlucts(gnm)  # mean square fluctuations based on all 20 GNM modes
writePDB('GNM_msf_all20_on_bfactor.pdb', reference_calphas, beta=msf_all20)
msf_first10 = calcSqFlucts(gnm[:10])  # mean square fluctuations based on the first 10 GNM modes
writePDB('GNM_msf_first10_on_bfactor.pdb', reference_calphas, beta=msf_first10)
msf_first3 = calcSqFlucts(gnm[:3])  # mean square fluctuations based on the first 3 GNM modes
writePDB('GNM_msf_first3_on_bfactor.pdb', reference_calphas, beta=msf_first3)
print("Mean square fluctuations computed with sets of GNM modes stored in the B-factor column:\n"
      "GNM_msf_all20_on_bfactor.pdb, GNM_msf_first10_on_bfactor.pdb, GNM_msf_first3_on_bfactor.pdb")

# generate a trajectory by the method traverseMode
# This takes steps in both directions starting from the provided structure to generate conformers along the chosen mode.
# generate a trajectory for all-atom GNM1 (does not work as mode is not 3D)
#gnm1_traj_aa = traverseMode(gnm_aa[0], atoms_all)
#gnm1_traj_aa.setAtoms(atoms_all)
# write PDB
#writePDB('GNM1_traj_aa.pdb', gnm1_traj_aa)



### Dynamical Domain Decomposition of reference_structure (using GNM)

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
