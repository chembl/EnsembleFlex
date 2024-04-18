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
#from matplotlib.pylab import *


#------- File Argument/Option Parser -------#
parser = argparse.ArgumentParser(description="Perform flexibility analysis on protein structures.")
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

#------- Logging to file and console -------#
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("logfile.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass
sys.stdout = Logger()
#-------------------------------------------#


# list of input pdb files
pdbfiles = glob.glob(input_path+"/*.pdb")
# parsing structures
structures = parsePDB(pdbfiles, subset='ca', compressed=False)  #, subset='ca'
reference_structure = structures[0] # by default first structure of ensemble is taken as reference structure and
# For unresolved atoms, the coordinates of the reference structure is assumed in RMSD calculations and superpositions.
# select C-alphas of reference structure
#reference_calphas = reference_structure.select('calpha')

# buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
# The reference structure is automatically the first member of list provided.
# Here superpositioning is omitted by setting it to 'False'.
ensemble = buildPDBEnsemble(structures, superpose=False)
# set c-alpha reference
#ensemble.setAtoms(reference_structure.calpha)


### 1. RMSF
rmsf = ensemble.getRMSFs()
plt.plot(rmsf);
plt.xlabel('Residue index');
plt.ylabel('RMSF');
plt.savefig('RMSF.png')
plt.close()
print('\nRoot mean square fluctuation (RMSF) plot saved to RMSF.png.\n')

# Save reference structure with RMSF in B-factor column
writePDB('RMSFonReference.pdb', reference_structure, beta=rmsf)
print('\nRMSF values are saved in B-factor column of RMSFonReference.pdb for visualisation.\n')



### 2. Principal Component Analysis (PCA)
# PCA is a method that identifies the components which account for the greatest amount of variability in the ensemble.
# PCA calculations

pca = PCA('')           # Instantiate a PCA instance
pca.buildCovariance(ensemble)   # Build covariance for the ensemble
pca.calcModes()                 # Calculate modes (20 of the by default)

# if more than 3 structures are present:
if len(pdbfiles) > 3:
    # using singular value decomposition for faster and more memory efficient calculation of principal modes
    pca_svd = PCA('ensemblePCAsvd')
    pca_svd.performSVD(ensemble)

    abs(pca_svd.getEigvals()[:20] - pca.getEigvals()).max()
    abs(calcOverlap(pca, pca_svd).diagonal()[:20]).min()

    # plot
    showOverlapTable(pca, pca_svd[:20])#;
    plt.savefig('Overlap_PCA_PCA-SVD.png')
    plt.close()
    print('\nOverlap plot of PCA and PCA with singular value decomposition saved to Overlap_PCA_PCA-SVD.png.\n')

else:
    print("Ensembles must contain more than 3 structures to perform singular value decomposition. \n")

## Variance along PCs
# print to sys.stdout
print('Variance along PCs:')
for mode in pca[:3]:
    var = calcFractVariance(mode)*100
    print(str(mode) + '  % variance = {:.2f}'.format(var))
    #print('{0:s}  % variance = {1:.2f}'.format(mode, var))
# save to file
with open('PCA_variance.txt', 'w') as f:
    print('Variance along PCs:', file=f)
    for mode in pca[:3]:
        var = calcFractVariance(mode) * 100
        print(str(mode) + '  % variance = {:.2f}'.format(var), file=f)

## Collectivity of modes
# print to sys.stdout
print('PCA mode collectivity:')
for mode in pca[:3]:    # Print PCA mode collectivity
    coll = calcCollectivity(mode)
    print(str(mode) + '  collectivity = {:.2f}'.format(coll))
    #print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))
# save to file
with open('PCA_mode_collectivity.txt', 'w') as f:
    print('PCA mode collectivity:', file=f)
    for mode in pca[:3]:  # Print PCA mode collectivity
        coll = calcCollectivity(mode)
        print(str(mode) + '  collectivity = {:.2f}'.format(coll), file=f)

writeArray('PCA_eigvecs.txt', pca.getEigvecs() ) # PCA eigenvectors
print('\nPCA eigenvectors are written to PCA_eigvecs.txt.\n')

# Square fluctuations - PCA
showSqFlucts(pca[:3])#;
plt.savefig('Fluctuations_PCA.png')
plt.close()
print('\nPCA square fluctuations plot is saved to Fluctuations_PCA.png.\n')


## Cross Correlations
# We can also see how correlated the motions for each residue are with each other residue. We see similar patterns
# for the two methods, especially when using a large number of modes.

showCrossCorr(pca[0])#;
plt.savefig('CrossCorrelations_PCA1perResidue.png')
plt.close()
print('\nCross Correlations of PC1 per residue is saved to CrossCorrelations_PCA1perResidue.png.\n')



### 3. Essential Dynamics Analysis (EDA)

eda = EDA('')
eda.buildCovariance( ensemble )
eda.calcModes()
#saveModel(eda)

## Variance along modes
# print to sys.stdout
print('Variance along modes:')
for mode in eda[:3]:
    var = calcFractVariance(mode)*100
    print(str(mode) + '  % variance = {:.2f}'.format(var))
    #print('{0:s}  % variance = {1:.2f}'.format(mode, var))
# save to file
with open('EDA_variance.txt', 'w') as f:
    print('Variance along modes:', file=f)
    for mode in eda[:3]:
        var = calcFractVariance(mode) * 100
        print(str(mode) + '  % variance = {:.2f}'.format(var), file=f)

# Square fluctuations - PCA
showSqFlucts(eda[:3])#;
plt.savefig('Fluctuations_EDA.png')
plt.close()
print('\nEDA square fluctuations plot is saved to Fluctuations_EDA.png.\n')



### Generate Report
import analyse_flex_prody_reporting as reporting
reporting.generate_report(str_input_path=args.input_path, output_path=output_path)
