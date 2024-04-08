#!/usr/bin/env python

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
structures = parsePDB(pdbfiles, compressed=False)  #, subset='ca'
reference_structure = structures[0] # by default first structure of ensemble is taken as reference structure and
# For unresolved atoms, the coordinates of the reference structure is assumed in RMSD calculations and superpositions.

# buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
# The reference structure is automatically the first member of list provided.
# Here superpositioning is omitted by setting it to 'False'.
ensemble = buildPDBEnsemble(structures, superpose=False)

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



### 3. Anisotropic Network Model (ANM) Normal Mode Analysis (NMA)
# The ANM allows for the identification of the most impactful (slowest) modes in dynamics of a single protein,
# which we can compare to the principal components from PCA.
# ANM calculations

anm = ANM('')                 # Instantiate a ANM instance
anm.buildHessian(reference_structure)   # Build Hessian for the reference chain
anm.calcModes()                   # Calculate slowest non-trivial 20 modes

## Variance along PCs
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

## Collectivity of modes
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


## PCA - ANM overlap
# We can also look at how well the modes produced from each method correlate with each other using the
# overlap (correlation cosine).
printOverlapTable(pca[:3], anm[:3]) # Top 3 PCs vs slowest 3 ANM modes

# Plotting overlap of PCA and NMA

showOverlapTable(pca[:6], anm[:6]); # Top 6 PCs vs slowest 6 ANM modes
#plt.title('PCA - ANM Overlap Table');
plt.savefig('OverlapTable_PCA_ANM.png')
plt.close()
print('\nOverlap table of PCA and ANM is saved to OverlapTable_PCA_ANM.png.\n')

# The cumulative overlap is the square root of the sum of squared overlaps.
showOverlap(pca[0], anm)#; # PC1 with all ANMs
showCumulOverlap(pca[0], anm, color='r')#;
plt.savefig('OverlapCumulative_PC1_ANM.png')
plt.close()
print('\nThe cumulative overlap table of PC1 with all ANMs is saved to OverlapCumulative_PC1_ANM.png.\n')

# Square fluctuations along PCA and ANM modes in the same plot
# ANM modes are scaled to have the same mean as PCA modes with showScaledSqFlucts.
# Alternatively, we could plot normalized square fluctuations with showNormedSqFlucts
showScaledSqFlucts(pca[0], anm[0])#;
plt.legend()#;
plt.savefig('SquareFluctuations_PC1_ANM1.png')
plt.close()
print('\nThe combined square fluctuations plot along PC1 and ANM1 is saved to SquareFluctuations_PC1_ANM1.png.\n')

showScaledSqFlucts(pca[1], anm[1])#;
plt.legend()#;
plt.savefig('SquareFluctuations_PC2_ANM2.png')
plt.close()
print('\nThe combined square fluctuations plot along PC2 and ANM2 is saved to SquareFluctuations_PC2_ANM2.png.\n')

if len(pca) >= 3:
    showScaledSqFlucts(pca[2], anm[2])#;
    plt.legend()#;
    plt.savefig('SquareFluctuations_PC3_ANM3.png')
    plt.close()
    print('\nThe combined square fluctuations plot along PC3 and ANM3 is saved to SquareFluctuations_PC3_ANM3.png.\n')

# Project the ensemble onto PC 1 and 2 using showProjection()
showProjection(ensemble, pca[:2])#;
#plt.axis([-0.8, 0.8, -0.8, 0.8])#;
plt.savefig('Projection_PC1_2.png')
plt.close()
print('\nProjection of the ensemble onto PC 1 and 2 is saved to Projection_PC1_2.png.\n')

# cross-projection plot comparing PCA modes and ANM modes using showCrossProjection()
showCrossProjection(ensemble, pca[0], anm[2], scale="y");
plt.legend(loc='upper left');
#plt.plot([-0.8, 0.8], [-0.8, 0.8], 'k');
#plt.axis([-0.8, 0.8, -0.8, 0.8]);
plt.savefig('CrossProjection_PCA_ANM.png')
plt.close()
print('\nCross-projection plot comparing PCA modes and ANM modes is saved to CrossProjection_PCA_ANM.png.\n')

# correlation between these projections
pca_coords, anm_coords = calcCrossProjection(ensemble, pca[0], anm[2])
print('\nCorrelation between these projections:')
print(np.corrcoef(pca_coords, anm_coords))

## Cross Correlations
# We can also see how correlated the motions for each residue are with each other residue. We see similar patterns
# for the two methods, especially when using a large number of modes.

showCrossCorr(pca[0])#;
plt.savefig('CrossCorrelations_PCA1perResidue.png')
plt.close()
print('\nCross Correlations of PC1 per residue is saved to CrossCorrelations_PCA1perResidue.png.\n')

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


### Generate Report
import analyse_flex_prody_reporting as reporting
reporting.generate_report(str_input_path=args.input_path, output_path=output_path)
