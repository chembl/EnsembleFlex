#!/usr/bin/env python

import os
import argparse
import glob
from prody import *
import matplotlib.pylab as plt


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



### Anisotropic Network Model (ANM) Normal Mode Analysis (NMA)
# The ANM allows for the identification of the most impactful (slowest) modes in dynamics of a single protein,
# which we can compare to the principal components from PCA.
# ANM calculations

anm = ANM('')                 # Instantiate a ANM instance
anm.buildHessian(reference_structure)   # Build Hessian for the reference chain
anm.calcModes()                   # Calculate slowest non-trivial 20 modes


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

writeModes('ANM_modes.txt', anm) # This function is based on writeArray
print('\nAnisotropic Network Model (ANM) modes are written to ANM_modes.txt.\n')

# Square fluctuations - ANM
showSqFlucts(anm[:3])#;
plt.savefig('Fluctuations_ANM.png')
plt.close()
print('\nANM square fluctuations plot is saved to Fluctuations_ANM.png.\n')


## Cross Correlations
# We can also see how correlated the motions for each residue are with each other residue.
# We see similar patterns for the two methods (PCA and ANM), especially when using a large number of modes.

showCrossCorr(anm[0])#;
plt.savefig('CrossCorrelations_ANM1perResidue.png')
plt.close()
print('\nCross Correlations of ANM1 per residue is saved to CrossCorrelations_ANM1perResidue.png.\n')


### 4. Dynamical Domain Decomposition of reference_structure

# Dynamical domain decomposition using GNM modes
# The number of dynamical domains obtained is approximately equal to the number of modes used (here 8).
gnm, _ = calcGNM(reference_structure, selstr='all')
domains = calcGNMDomains(gnm[:8])
# assign this data to the AtomGroup:
reference_structure.setData('domain', domains)
# using the B-factor field for writing the domains
writePDB('reference_structure_dynamic_domains.pdb', reference_structure, beta=domains)
print('\nDynamical domains of reference structure are saved in B-factor column of reference_structure_dynamic_domains.pdb '
      'for visualisation.\n')

