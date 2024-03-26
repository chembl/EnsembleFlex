#!/usr/bin/env python

import os
import argparse
import glob
from prody import *
import matplotlib.pylab as plt
#from pylab import *  # Matplotlib package
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


# ## Protein structure bipartition # not working for ANM modes!
# showProtein(reference_calphas, mode=anm[0])
# plt.savefig('ProteinStructureBipartition_ANM1.png')
# plt.close()
# print('\nProtein structure bipartition plot based on ANM1 is saved to ProteinStructureBipartition_ANM1.png.\n')




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


## Protein structure bipartition
showProtein(reference_calphas, mode=gnm[0])
plt.savefig('ProteinStructureBipartition_GNM1.png')
plt.close()
print('\nProtein structure bipartition plot based on GNM1 is saved to ProteinStructureBipartition_GNM1.png.\n')
# writePDB('ProteinStructureBipartition_GNM1.pdb', reference_calphas, beta=gnm[0])
# print('\nProtein structure bipartition based on GNM1 is saved in B-factor column of ProteinStructureBipartition_GNM1.pdb.\n')




### Dynamical Domain Decomposition of reference_structure (using GNM)

# Dynamical domain decomposition using GNM modes
ddd_gnm, _ = calcGNM(reference_structure, selstr='all')
# The number of dynamical domains obtained is approximately equal to the number of modes used (here 8).
domains = calcGNMDomains(ddd_gnm[:3])
# assign this data to the AtomGroup:
reference_structure.setData('domain', domains)
# using the B-factor field for writing the domains
writePDB('reference_structure_dynamic_domains_GNM.pdb', reference_structure, beta=domains)
print('\nDynamical domains of reference structure are saved in B-factor column of reference_structure_dynamic_domains_GNM.pdb '
      'for visualisation.\n')
