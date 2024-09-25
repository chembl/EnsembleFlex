
# EnsembleFlex - Flexibility Analysis of Structure Ensembles

## Overview

This program provides flexibility analysis tools for protein structure ensembles via a graphical user interface (GUI).
Nevertheless, all tools are provided as separately executable Python or R scripts and can therefore be integrated in 
custom in-house pipelines, without using the GUI. 
For more details, please see section 'Operating instructions' below, as well as the additional `user_guide.md` file.  

**Input**: The program and tools are particularly designed to work with heterogeneous pdb ensembles (non-identical 
number of atoms or residues), but can also be used for homogenous ensembles (originating e.g. from molecular dynamics 
simulations). The analysis will be performed for a single chain (monomeric) protein. If your protein of interest 
contains several chains, they need to be split first (e.g. using the provided tool `split_pdbs_bio3d.R`) and the 
equivalent chains of the ensemble need to be placed together into one folder for analysis. 

### Version
v1.0.0 (2024.07)

## General installation notes
EnsembleFlex is basically composed of Python and R scripts that require several Python and R packages. 
You have different options for installing EnsembleFlex.
Choose one of the following options depending on your expertise and required environment restrictions. 
Here is the list of options ordered from highest environment reproducibility to lowest:
- Using **Docker** - This is the most reproducible environment setup and comes with all benefits but also some 
inconveniences of containers.
- Using the `conda-lock.yml` file with **Conda-lock** or **Micromamba** - This provides also a highly reproducible 
environment, but can only be used on Mac or Linux systems (not for Windows).
- Using the `environment_versioned.yml` file with **Conda** or **Mamba** to create the environment - This ensures 
package compatibility, but may give you slightly older package versions compared to the next option below.
- Using the `environment.yml` file with **Conda** or **Mamba** to create the environment - This installs more recent 
versions of required packages.
- [NOT recommended] Installing all packages listed in the `environment.yml` **manually** on your system / in your custom environment - 
This is the highest customizable installation option, but you need to manage the package dependencies.


## Installation \& launch with DOCKER

**Installation**:

        git clone https://gitlab.ebi.ac.uk/melanie/ensembleflex.git
        cd EnsembleFlex
        docker build -t ensembleflex-image .

**Launch**:  

        docker run -d -it --name ensembleflex -p 8501:8501 -p 80:80 --volume ./docker-data:/app/docker-data ensembleflex-image

Then open your browser and navigate to `http://localhost:8501/`. There you should see the graphical user interface.
Note that your input pdb files need to be put into `EnsembleFlex/docker-data/` (or subdirectories), 
as this is a mounted directory in the Docker image and will therefore be accessible from within the Docker container 
and the user interface.


## Installation \& launch with CONDA/MAMBA

EnsembleFlex can be installed by cloning this repository and setting up an environment using your favourite environment 
manager. (I recommend [`mamba` through miniforge](https://github.com/conda-forge/miniforge#mambaforge)).

[//]: # (* Installation:)
Clone the repository and go to the EnsembleFlex folder:

      git clone https://gitlab.ebi.ac.uk/melanie/ensembleflex.git
      cd EnsembleFlex
  
* install with `mamba` (if you prefer `conda`, just replace 'mamba' with 'conda'):
  
      mamba env create -f environment.yml
  
  **Troubleshooting**: If you encounter package incompatibilities try the next installation suggestion, which is the recommended way. 
  If you are on Windows, replace `environment.yml` with `environment_versioned.yml` in the previous command.

* OR (recommended): On **Linux** or **macOS**: For an exact environment reproduction, to make sure that there are no 
package version inconsistencies, you can use the provided 'conda-lock' file (only on Linux or macOS). Be aware that 
this needs to be installed with the `conda-lock install` or `micromamba install` commands, as `mamba/conda create` will 
ignore the pip package dependencies.
  - using ['conda-lock'](https://conda.github.io/conda-lock/) (will be installed):

        mamba install -c conda-forge conda-lock
        conda-lock install --name ensembleflex --file conda-lock.yml
  
  - OR using ['micromamba'](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) (to be installed):
  
        micromamba install --name ensembleflex --file conda-lock.yml

Finally:
Install the R package(s) not available from conda: This is particularly `bio3d.eddm`. In case installation fails you can 
still run EnsembleFlex, only ensemble Difference Distance Matrix (eDDM) calculations won't be available within ensemble 
analysis.

    mamba activate ensembleflex
    Rscript requirements.R


**Launch**:  
With conda/mamba installation EnsemblFlex can be used with the 
browser based Graphical User Interface (GUI) and the 
Command Line (see detailed documentation in `user_guide.md`).  

First, activate the environment:
      
      mamba activate ensembleflex

  For the browser-based graphical user interface run (replace `path/to/` with yor local path):

      streamlit run path/to/EnsembleFlex/src/streamlit_app/streamlit_app.py

  The user interface should automatically open in a new tab of your default browser.



## Operating instructions  

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

When using the graphical user interface analysis results are presented in a scrollable container which will appear when 
analysis is finished.  

### Analysis
Please look into `user_guide.md` for a detailed interpretation help and description of methods.  

Primary tools and packages used in the workflow are:
- [Bio3D](http://thegrantlab.org/bio3d/) (R package) - main flexibility analysis & prediction
- [ProDy](http://www.bahargroup.org/prody/) (Python package) - alternative main flexibility analysis & prediction
- [Biopython](https://biopython.org) (Python package) - SASA & radius of gyration analysis, data handling tools
- [umap](https://github.com/tkonopka/umap) (R package) - UMAP projection
- [cluster](https://cran.r-project.org/web/packages/cluster/index.html) and [clvalid](https://cran.r-project.org/web/packages/clValid/index.html) - clustering and cluster validation
- [vanddraabe](external_packages) (R package - modified package provided in `external_packages`) - conserved water analysis

for the user interface:
- [Streamlit](https://streamlit.io) - main interface
- [py3Dmol](https://pypi.org/project/py3Dmol/) - molecular visualisations


#### Workflow description:  

1. **Settings**: The user needs to provide input and output directories.

2. **Structural superposing**: The analysis result is in many cases (when coordinates are used) dependent on the relative superposition of the 
structures contained in the ensemble. Therefore, structural superposing (also known as structural alignment) needs to 
be performed prior to analysis. 
The user can directly supply a superposed ensemble, or choose one of the two provided superposing methods - 
from Bio3D (recommended) or ProDy. 

3. **Flexibility Analysis**: Upon superposing the ensemble can be analysed using Bio3D tools (RMSD, RMSF, dimension reduction and projection with 
PCA and UMAP, ...) and Biopython tools (SASA, radius of gyration). 
(Additionally, equivalent ProDy tools are provided via a command line script.).  

4. **Binding Site Analysis** [optional]: For binding site analysis the user has the option to perform a superposing on binding site residues only (recommended 
only when there is domain movement shifting the binding site). 
Again this is possible with Bio3D or ProDy. The user can also opt for using the globally superimposed structures 
(recommended in most minor-movement cases).
The binding site analysis includes an automatised identification of the binding site (based on user-defined distance to ligands), 
a statistical analysis of the residues identified as binding site (occurrence), as well as an all-atom flexibility analysis.
Additionally a conserved water analysis can be performed, which takes into account all waters from all provided structures (not only the binding site), 
and outputs an additional PyMol script dedicated to the binding site of the reference structure.

5. **Flexibility Prediction** [optional]: Additionally, there is the option to perform flexibility prediction based on elastic network model Normal Mode Analysis
(NMA), as well as Essential Site Scanning Analysis (ESSA) on the reference structure.


## File manifest  
```
EnsembleFlex  
├── LICENSE  
├── README.md  
├── docker-data
├── documentation  
│   └── OutputDirectoryStructure.rtf
├── environment.yml
├── example_input_data
├── example_results
├── external_packages
│   └── vanddraabe.zip
├── requirements.r
├── requirements.txt
├── src  
│   ├── analyse_conserved_waters.R
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

## Compliance with FAIR principles for research software
- **F**indable: The EnsembleFlex code is available on GitHub, the Docker image is also available on Docker Hub 
and the tool is listed at Elixir [bio.tools](https://bio.tools/).
- **A**ccessible: All code is available on GitHub and can be installed via Conda environment manager or simply 
executed as Docker image. 
Documentation is provided via this `README.md` and an additional `user_guide.md` file.
- **I**nteroperable: The code is constructed in a building block manner, where all steps of the workflow can be 
executed independently via command line execution of scripts having a consistent execution syntax (see `user_guide.md`), 
ensuring for customizability and interoperability. 
- **R**eusable: All provided scripts are well documented with execution information in the file headers. 
The provided scripts are wrapped into a complete workflow, executable through a browser-based graphical user interface, 
that does not require any coding experience. 
Environment reproducibility is achieved through installation with a conda-lock file and full isolation through Docker.

## Copyright and licensing information  
This program is distributed under a permissive open source licence
(MIT). A copy is provided in file `LICENSE`.

## Contact information  
The author Melanie Schneider can be contacted at melanie@ebi.ac.uk.

## Funding
Funding was provided to Melanie Schneider as ARISE Fellowship from the European Union’s Horizon 2020 research and 
innovation programme under the Marie Skłodowska-Curie grant agreement No 945405.