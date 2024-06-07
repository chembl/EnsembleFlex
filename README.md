
# EnsembleFlex - Flexibility Analysis of Structure Ensembles

## Overview

This program provides flexibility analysis tools for protein structure ensembles via a graphical user interface (GUI).
Nevertheless, all tools are provided as separately executable Python or R scripts and can therefore be integrated in 
custom in-house pipelines, without using the GUI.
The program and tools are particularly designed to work with heterogeneous pdb ensembles (non-identical number of atoms 
or residues), but can also be used for homogenous ensembles (originating e.g. from molecular dynamics simulations).
The analysis will be performed for a single chain (monomeric) protein. If your protein of interest contains several 
chains, they need to be split first (e.g. using the provided tool `pdb_splitchain` from the 
[pdb-tools](http://www.bonvinlab.org/pdb-tools/) utility suit) and the equivalent chains of the ensemble need to be placed together into one folder for 
analysis. 

### Version
v1.0.0 2024.0

### Installation

Using ANACONDA:
can be installed by cloning this repository and setting up an environment using your favourite environment manager 
(I recommend [mamba](https://github.com/conda-forge/miniforge#mambaforge)).

* Installation:

      git clone https://github.com/EnsembleFlex.git
      cd EnsembleFlex
      mamba env create -f environment.yml

  OR with conda (instead of mamba):

      conda env create -f environment.yml

* Usage: With conda/mamba installation EnsemblFlex can be used with the 
[Graphical User Interface (GUI) documentation](https://ensembleflex.readthedocs.io/en/latest/modules.html) and the 
[Command Line documentation](https://ensembleflex.readthedocs.io/en/latest/command_line.html)
      
      mamba activate ensembleflex

  For the webapp GUI:

      streamlit run path/to/EnsembleFlex/src/streamlit_app/streamlit_app.py



Using DOCKER:

* Installation:

        docker pull quay.io/...

* Usage:

        docker run quay.io/... <command>
...

### Operating instructions  
#### Input files
After launching the GUI the user needs to provide an input directory containing all (monomeric) ensemble structure 
files that shall be analysed. Example input data is provided in the folder `example_input_data`.

#### Output files
The execution of provided tools will produce the following types of output files:
- analysis plots and graphs in png format, 
- pdb structure files with modified coordinates or different information in the b-factor column,
- structure ensemble pdb files
- PyMol scripts (to be executed for visualization in PyMol)  

The user is free to use any convenient structure visualisation tool to examine the structural output (except for the 
PyMol scripts, which are only executable with PyMol).
Example output data is provided in the folder `example_results`.

#### Analysis
Please look into `user_guide.md` for a detailed interpretation help and description of methods.  

- Workflow description:  

The analysis result is in many cases (when coordinates are used) dependent on the relative superposition of the 
structures contained in the ensemble. Therefore, structural superposing (also known as structural alignment) needs to 
be performed prior to analysis. 
The user can directly supply a superposed ensemble, or choose one of the two provided superposing methods - 
from Bio3D (recommended) or ProDy. 
Upon superposing the ensemble can be analysed using Bio3D tools. 
(Additionally, equivalent ProDy tools are provided via a command line script.).  

When using the user interface analysis results are presented in a scrollable container which will appear when analysis 
is finished.  

For binding site analysis the user has the option to perform a superposing on binding site residues only (recommended 
when there is domain movement shifting the binding site). 
Again this is possible with Bio3D or ProDy. The user can also opt for using the globally superimposed structures 
(recommended in most minor-movement cases).

Additionally, there is the option to perform flexibility prediction based on elastic network model Normal Mode Analysis
(NMA), as well as Essential Site Scanning Analysis (ESSA) on the reference structure.

### File manifest  
```
EnsembleFlex  
├── LICENSE  
├── README.md  
├── documentation  
│   └── OutputDirectoryStructure.rtf
├── environment.yml
├── example_input_data
├── example_results
├── requirements.r
├── requirements.txt
├── src  
│   ├── analyse_flex_binding_site_bio3d.R
│   ├── analyse_flex_bio3d.R
│   ├── analyse_flex_bio3d_reporting.py
│   ├── analyse_flex_prody.py
│   ├── analyse_flex_prody_reporting.py
│   ├── analyse_flex_sasa_biopython.py
│   ├── identify_binding_site_bio3d.R
│   ├── predict_flex_binding_site_essa_prody.py
│   ├── predict_flex_nma_bio3d.R
│   ├── predict_flex_nma_prody.py
│   ├── superimpose_binding_site_bio3d.R
│   ├── superimpose_bio3d.R
│   ├── superimpose_prody.py
│   ├── streamlit_app
│   │   ├── directorypicker.py
│   │   ├── requirements.txt
│   │   ├── stoc.py
│   │   ├── streamlit_app.py
│   └── tools
│       ├── data_on_structure.R
│       ├── pdb_delhetatm.py
│       ├── pdb_delhetatm_ions.py
│       ├── pdb_delresname.py
│       ├── run_pdb_del_on_directory.sh
│       ├── sort_pdbs_from_dataframe.py
│       ├── sort_pdbs_from_file.py
│       ├── sort_pdbs_has_gap_in_seq.py
│       ├── sort_pdbs_has_ligand.R
│       ├── split_pdbs_bio3d.R
│       ├── subset_alignment_gap_fraction_drop.py
│       └── subset_alignment_has_gap_in_seq.sh
└── tests
```


### Copyright and licensing information  
This program is distributed under a permissive open source licence 
(GNU General Public License, version 3 (GPL-3.0)). A copy is provided in file `LICENSE_GPL-3.txt`.

### Contact information  
The author Melanie Schneider can be contacted at melanie@ebi.ac.uk.