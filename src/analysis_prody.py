
import glob
import os
from optparse import OptionParser
from prody import *

parser = OptionParser()
parser.add_option("-i", "--in", dest="input_path",
                  help="input dataset directory path", metavar="PATH")
parser.add_option("-o", "--out", dest="output_dir", default="outdir",
                  help="output directory name", metavar="STRING")

(opts, args) = parser.parse_args()
if os.path.isdir(opts.input_path) and os.path.isdir(opts.output_dir):
    output_path = os.path.abspath(opts.output_dir)
elif os.path.isdir(opts.input_path) and not os.path.isdir(opts.output_dir):
    if os.path.isdir(os.path.dirname(opts.output_dir)):
        output_path = opts.output_dir
        os.mkdir(output_path)
    else:
        output_path = opts.input_path+"/"+opts.output_dir
    if not os.path.exists(output_path):
        os.mkdir(output_path)
else:
    print("Please specify a valid input directory path to your structure files.\n")

os.chdir(output_path)

# list of pdb files
pdbfiles = glob.glob(opts.input_path+"/*.pdb")
# parsing structures
structures = parsePDB(pdbfiles, compressed=False)  #, subset='ca'

# buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
# The reference structure is automatically the first member of list provided.
ensemble = buildPDBEnsemble(structures)


### PCA calculations

pca = PCA('ensemblePCA')           # Instantiate a PCA instance
pca.buildCovariance(ensemble)   # Build covariance for the ensemble
pca.calcModes()                 # Calculate modes (20 of the by default)

# using singular value decomposition for faster and more memory efficient calculation of principal modes
pca_svd = PCA('ensemblePCAsvd')
pca_svd.performSVD(ensemble)

abs(pca_svd.getEigvals()[:20] - pca.getEigvals()).max()
abs(calcOverlap(pca, pca_svd).diagonal()[:20]).min()

# plot
showOverlapTable(pca, pca_svd[:20]).savefig('PCA_PCA-SVD_overlap.png')#;

## Variance along PCs
for mode in pca[:3]:
    var = calcFractVariance(mode)*100
    print('{0:s}  % variance = {1:.2f}'.format(mode, var))

## Collectivity of modes
for mode in pca[:3]:    # Print PCA mode collectivity
    coll = calcCollectivity(mode)
    print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))


writeArray('PCA_eigvecs.txt', pca.getEigvecs() ) # PCA eigenvectors



### ANM calculations

anm = ANM('ANM')                 # Instantiate a ANM instance
anm.buildHessian(structures[1])   # Build Hessian for the reference chain
anm.calcModes()                   # Calculate slowest non-trivial 20 modes

## Collectivity of modes
for mode in anm[:3]:    # Print PCA mode collectivity
    coll = calcCollectivity(mode)
    print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))

writeModes('ANM_modes.txt', anm) # This function is based on writeArray


## Plotting overlap of PCA and NMA

showOverlapTable(pca[:6], anm[:6]);
title('PCA - ANM Overlap Table');

# Square fluctuations
showSqFlucts(pca[:3]).savefig('PCA_fluctuations.png')#;
showSqFlucts(anm[:3]).savefig('ANM_fluctuations.png')#;

# square fluctuations along PCA and ANM modes in the same plot
# ANM modes are scaled to have the same mean as PCA modes with showScaledSqFlucts.
# Alternatively, we could plot normalized square fluctuations with showNormedSqFlucts
showScaledSqFlucts(pca[0], anm[2]);
legend();

showScaledSqFlucts(pca[1], anm[0]);
legend();

showScaledSqFlucts(pca[0], anm[2]);
legend();


# project the ensemble onto PC 1 and 2 using showProjection()
showProjection(ensemble, pca[:2]);
#axis([-0.8, 0.8, -0.8, 0.8]);

# cross-projection plot comparing PCA modes and ANM modes using showCrossProjection()
showCrossProjection(ensemble, pca[0], anm[2], scale="y");
#legend(loc='upper left');
#plot([-0.8, 0.8], [-0.8, 0.8], 'k');
#axis([-0.8, 0.8, -0.8, 0.8]);

# correlation between these projections
pca_coords, anm_coords = calcCrossProjection(ensemble, pca[0], anm[2])
print(np.corrcoef(pca_coords, anm_coords))

