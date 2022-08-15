
import glob
import os
from optparse import OptionParser
# import sys
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
    output_path = opts.input_path+"/"+opts.output_dir
    if not os.path.exists(output_path):
        os.mkdir(output_path)
else:
    print("Please specify a valid input directory path to your structure files.\n")

os.chdir(output_path)

pdbfiles = glob.glob(opts.input_path+"/*.pdb")

structures = parsePDB(pdbfiles, compressed=False) #, subset='ca'

# buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
# The reference structure is automatically the first member of list provided.
ensemble = buildPDBEnsemble(structures)
# Perform an iterative superimposition
ensemble.iterpose()
# Save coordinates
writePDB('ensemble.pdb', ensemble)
print("Superimposed ensemble coordinates are saved to ensemble.pdb\n")

print(output_path)

if not os.path.exists(output_path+"/split_ensemble"):
    os.mkdir(output_path+"/split_ensemble")
os.chdir(output_path +"/split_ensemble")

from tools import pdb_splitmodel
pdb_splitmodel.run(output_path+"/ensemble.pdb")

print("Multimodelfile split and saved to directory split_ensemble\n")
