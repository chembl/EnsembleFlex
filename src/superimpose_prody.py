#!/usr/bin/env python

"""
Ensemble structural superpositioning using the Python package ProDy.

Usage:
    python3 superimpose_prody.py -i <input_directory> -o <output_directory>

Example:
    python3 superimpose_prody.py -i pdbs -o EnsemblFlex
"""

import glob
import os
import argparse
import subprocess
from prody import *

#------- File Argument/Option Parser -------#
parser = argparse.ArgumentParser(description="Perform global iterative superpositioning of protein structures.")
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

toolsdir = os.path.dirname(os.path.realpath(__file__))+"/tools"

# change working directory to output_path
os.chdir(output_path)

if not os.path.exists(output_path+"/superimposed"):
    os.mkdir('superimposed')
os.chdir('superimposed')

# list of input pdb files
pdbfiles = glob.glob(input_path+"/*.pdb")

# parsing structures
structures = parsePDB(pdbfiles, compressed=False) #, subset='ca'

# buildPDBEnsemble() maps each structure against the reference structure using a function such as mapOntoChain().
# The reference structure is automatically the first member of list provided.
ensemble = buildPDBEnsemble(structures)
# Perform an iterative superimposition:
# Iteratively superpose the ensemble until convergence. Initially, all conformations are aligned with the reference
# coordinates. Then mean coordinates are calculated, and are set as the new reference coordinates. This is repeated
# until reference coordinates do not change. This is determined by the value of RMSD between the new and old reference
# coordinates. Note that at the end of the iterative procedure the reference coordinate set will be average of
# conformations in the ensemble.
ensemble.iterpose()
# Save coordinates
writePDB('ensemble.pdb', ensemble)
print("Superimposed ensemble coordinates are saved to ensemble.pdb.\n ")

# Using pdb_splitmodel.py from tools directory to split the ensemble file into directory "split_ensemble"

if not os.path.exists(output_path+"/superimposed/split_ensemble"):
    os.mkdir(output_path+"/superimposed/split_ensemble")
os.chdir(output_path +"/superimposed/split_ensemble")

print("toolsdir = "+toolsdir)
subprocess.call(['python3', '{}/pdb_splitmodel.py'.format(toolsdir), '{}/superimposed/ensemble.pdb'.format(output_path)])

print("\nMultimodelfile split and saved to directory /superimposed/split_ensemble.\n ")
