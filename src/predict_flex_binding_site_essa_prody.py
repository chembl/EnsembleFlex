#!/usr/bin/env python
'''
Essential Site Scanning Analysis (ESSA) - using ProDy
'''
# A ProDy tutorial on ESSA is provided here: http://prody.csb.pitt.edu/_static/ipynb/workshop2021/prody_essa.html

import os
import argparse
import glob
import numpy as np
from prody import *
import matplotlib.pylab as plt

plt.ion()  # The matplotlib.pyplot.ion() function turns on the interactive mode


#------- File Argument/Option Parser -------#
parser = argparse.ArgumentParser(description="Perform Essential Site Scanning Analysis (ESSA) on a protein structure.")
parser.add_argument("-i", "--input", dest="input_path", required=True,
                    help="input dataset directory path", metavar="PATH")
parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                    help="output directory name", metavar="STRING")
parser.add_argument("-b", "--bsresidues", dest="bsresidues", default=None,
                    help="binding site residues as file", metavar="STRING")

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

if os.path.isfile(args.bsresidues):
    bsresiduefile = os.path.abspath(args.bsresidues)
    print("Residue file in use: "+bsresiduefile)
else:
    print("No binding site residue list provided.\n")

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

# read binding site residues
with open(bsresiduefile) as f:
    residues = f.read()
#print(residues)


### Essential Site Scanning Analysis (ESSA)

# instantiate an ESSA object
essa = ESSA()
# set system
essa.setSystem(reference_structure)
# perform scanning
essa.scanResidues()

# Plot Z-Scores per residue with highlighted binding site residues
# ProDy selection string to highlight binding residues
prody_sele = 'resnum '+str(residues).replace("\n", " ")
print("\nResidue numbers used from file are: ")
print(prody_sele)
# The blue dashed baseline shows the q-th quantile of the profile, which is by default q=0.75, representing the top quartile.
with plt.style.context({'figure.figsize': (9, 7), 'figure.dpi': 100}):
    essa.showESSAProfile(rescode=True, sel=prody_sele)
plt.savefig('ESSA_profile_of_reference_structure.png')
plt.close()
print('\nESSA profile plot is saved to ESSA_profile_of_reference_structure.png.')

# Save ESSA z-scores
essa.saveESSAZscores()
print("\nESSA Z-scores written to file")

# Save PDB file with Z-scores in the B-factor column
essa.writeESSAZscoresToPDB()
print("\nPDB file with ESSA Z-scores in B-factor column written to file\n")
