#!/bin/bash

# Runs the specified pdb_delXXX.py script on all pdb files in a specified directory
# Usage: run_pdb_del_on_directory.sh <path_to_pdb_del_script> <path_to_directory> <option> <outpath(optional)>

# to be run with scripts pdb_delelem.py, pdb_delresname.py, pdb_delhetatm.py, pdb_delhetatm_ions.py, and pdb_delres.py

# example1: run_pdb_del_on_directory.sh ~/ARISE/ARISE-Project/EnsembleFlex/src/tools/pdb_delelem.py \
#EnsemblFlex/superimposed MN EnsemblFlex/superimposed_no_MN

# example2: run_pdb_del_on_directory.sh ~/ARISE/ARISE-Project/EnsembleFlex/src/tools/pdb_delhetatm.py \
#EnsemblFlex/superimposed hetatm EnsemblFlex/superimposed_no_hetatm

# example3: run_pdb_del_on_directory.sh ~/ARISE/ARISE-Project/EnsembleFlex/src/tools/pdb_delhetatm_ions.py \
#EnsemblFlex/superimposed ions EnsemblFlex/superimposed_no_ions

# example4: run_pdb_del_on_directory.sh ~/ARISE/ARISE-Project/EnsembleFlex/src/tools/pdb_delres.py \
#EnsemblFlex/split_pdbs_filtered negative EnsemblFlex/superimposed_no_negative


# path_to_script=$1
# path_to_directory=$2
# option=$3
# outpath=$4 (optional)


if [ "$#" -eq 4 ]; then # if number of arguments equals 4
  outpath=$4
else
	outpath=$2_no_$3
fi

mkdir $outpath

if [ "$3" = "ions" ] || [ "$3" = "hetatm" ]; then # if <option> argument is "ions" or "hetatm"
  for f in $2/*.pdb; do
    python $1 "$f" > "$outpath/$(basename "$f")" # <option> argument is not used in execution
  done
elif [ "$3" = "negative" ]; then # if <option> argument is "negative"
  for f in $2/*.pdb; do
    # <option> argument "negative" is not used in execution, instead -:0 (which deletes residues from START to 0)
    python $1 -:0 "$f" > "$outpath/$(basename "$f")"
  done
else
  for f in $2/*.pdb; do
    python $1 -$3 "$f" > "$outpath/$(basename "$f")"
  done
fi