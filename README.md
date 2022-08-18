
# EnsembleFlex - Flexibility Analysis of Structure Ensembles

This program provides flexibility analysis tools for protein structure ensembles via a graphical user interface (GUI).
Nevertheless, all tools are provided as separately executable Python or R scripts and can therefore be integrated in 
custom in-house pipelines, without using the GUI.
The program and tools are particularly designed to work with heterogeneous ensembles (non-identical number of atoms or residues), 
but can also be used for homogenous ensembles (originating e.g. from molecular dynamics simulations).
The analysis will be performed for a single chain (monomeric) protein. If your protein of interest contains several 
chains, they need to be split first (e.g. using the provided tool `pdb_splitchain` from the 
[pdb-tools](http://www.bonvinlab.org/pdb-tools/) utility suit) and the equivalent chains of the ensemble need to be 
placed together into one folder for analysis. 

### Installation instructions  
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
The analysis result is in many cases (when coordinates are used) dependent on the relative superposition of the 
structures contained in the ensemble. Therefore, structural superposing (also known as structural alignment) needs to 
be performed prior to further analysis. 
The user can directly supply a superposed ensemble, or choose one of the two provided superposing methods - 
from Bio3D (recommended) or ProDy. 
Upon superposing the ensemble can be analysed using Bio3D and/or ProDy tools.
For binding site analysis a dedicated superposing on binding site residues is recommended (again possible with Bio3D 
or ProDy), but the user can also opt for using the globally superimposed structures.


### File manifest  
```
.  
├── LICENSE_GPL-3.txt  
├── README.md  
├── documentation  
├── example_input_data  
├── example_results  
├── src  
│   ├── analysis_bio3d.R  
│   ├── analysis_prody.py  
│   ├── main.py  
│   ├── styles.qss  
│   ├── superimpose_bio3d.R  
│   ├── superimpose_prody.py  
│   └── tools  
│       └── pdb_splitmodel.py  
└── tests  
```


### Copyright and licensing information  
This program is distributed under a permissive open source licence 
(GNU General Public License, version 3 (GPL-3.0)). A copy is provided in file `LICENSE_GPL-3.txt`.

### Contact information  
The author Melanie Schneider can be contacted at melanie@ebi.ac.uk.