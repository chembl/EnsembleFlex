#!/usr/bin/env python

import os
import argparse
import glob
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd

#------- File Argument/Option Parser -------#
parser = argparse.ArgumentParser(description="Perform SASA calculations on protein structures.")
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

parser = PDBParser(PERMISSIVE=True)
sr = ShrakeRupley()

pdbfiles = glob.glob(input_path+"/*.pdb")

sasa_df = pd.DataFrame(columns=['ResNum', 'ResName'])

#metals = ['ZN', 'FE', 'CU', 'MG', 'CA', 'MN']
residue_names = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
                 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

for fileName in pdbfiles:
    # parse a structure
    structure_id = fileName.rsplit('/', 1)[1][:-4]
    structure = parser.get_structure(structure_id, fileName)

    # remove water molecules (all atoms in a structure may affect SASA calculation)
    for chain in structure[0]:
        for residue in list(chain):
            if residue.id[0] == 'W':
                chain.detach_child(residue.id)

    # Calculate surface accessibility surface area for an entity
    # Level options: “A” (Atom), “R” (Residue), “C” (Chain), “M” (Model), or “S” (Structure)
    # if level=”R”, all residues will have a .sasa attribute.
    sr.compute(structure[0], level="R")

    # get SASA data from each structure into list
    my_list = []
    for chain in structure[0]:
        for res in chain:
            my_list.append((res.id[1], res.get_resname(), round(res.sasa, 2)))
    # combine lists into df
    my_df = pd.DataFrame(my_list, columns=['ResNum', 'ResName', structure_id])
    # keep only calculations for natural residues (filter other molecules)
    my_df = my_df[my_df['ResName'].isin(residue_names)]
    # merge data frames
    sasa_df = sasa_df.merge(my_df, on=['ResNum', 'ResName'], how='outer')

#print(sasa_df) # for debugging

# save the dataframe containing per residue SASA values for all structures as CSV file
sasa_df.to_csv('SASA_per_structure.csv', index=False)
print('\nSASA values per residue and structure are saved to SASA_per_structure.csv.\n')

# get numeric data to perform calculations
calc_df = sasa_df.loc[:, ~sasa_df.columns.isin(['ResNum', 'ResName'])]
# do calculations: max, min, spread, sd
sasa_global_df = sasa_df[['ResNum', 'ResName']]
sasa_global_df['max'] = calc_df.max(axis=1)
sasa_global_df['min'] = calc_df.min(axis=1)
sasa_global_df['spread'] = calc_df.max(axis=1)-calc_df.min(axis=1)
sasa_global_df['sd'] = calc_df.std(axis=1).round(decimals=2)
# save the dataframe containing global per residue calculations as CSV file
sasa_global_df.to_csv('SASA_global.csv', index=False)
print('\nSASA global metrics per residue are saved to SASA_global.csv.\n')
# save only sd to file
#sasa_global_df['sd'].to_csv('SASA_sd_global.csv', index=False, header=False)
