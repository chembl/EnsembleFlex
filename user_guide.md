
# User Guide
# EnsembleFlex - Flexibility Analysis of Structure Ensembles

## Introduction

Welcome to an exploration of the dynamics within protein structure ensembles!

EnsembleFlex serves as a exploration gateway into the field of protein structural flexibility 
in an ***automated*** manner, where tools are ***streamlined*** and analysis is easily ***reproducible*** and 
***fast***.

**Why investigating protein flexibility?**

    Proteins, the fundamental molecular machines exhibit a remarkable diversity in structure and function. 
    While classical structural biology has provided invaluable insights into the static three-dimensional 
    architectures of proteins, the dynamic nature of these macromolecules is increasingly recognized as a 
    critical determinant of their biological activity.  
    
    The inherent flexibility of proteins manifests at various levels, spanning from local fluctuations 
    of individual amino acid side chains to large-scale domain movements. Understanding the dynamic 
    behavior of proteins is essential for deciphering the mechanisms underlying their diverse functions, 
    including enzymatic catalysis, molecular recognition, and signal transduction. Moreover, protein 
    flexibility plays a pivotal role in modulating interactions with other biomolecules, enabling 
    proteins to adapt to different cellular environments and respond to external stimuli.
    
    Understanding protein structural flexibility is not just about unraveling molecular intricacies; 
    it holds the key to deciphering diseases, designing targeted therapies, and engineering proteins 
    with tailored functionalities.

## Input

- `<Input path>`: The primary input is the path to a folder that contains all PDB structure files to be analysed.
- Calculations can only be performed on monomers! 
- In case your PDB files contain several chains or are multi-model files they need to be split. 
(Tools to do this are provided - see below.)


## Command-line usage of analysis scripts explained

The following analysis scripts can be used independently from the browser-based interface.  
Note that script usage follows a general common rule (`[Rscript/python3] path/to/script -i <input> -o <output> <more options>`)
and is specified in the header of each script.
```
EnsembleFlex  
└── src  
   ├── analyse_flex_binding_site_bio3d.R
   ├── analyse_flex_bio3d.R
   ├── analyse_flex_bio3d_reporting.py
   ├── analyse_flex_prody.py
   ├── analyse_flex_prody_reporting.py
   ├── analyse_flex_sasa_biopython.py
   ├── identify_binding_site_bio3d.R
   ├── predict_flex_binding_site_essa_prody.py
   ├── predict_flex_nma_bio3d.R
   ├── predict_flex_nma_prody.py
   ├── superimpose_binding_site_bio3d.R
   ├── superimpose_bio3d.R
   ├── superimpose_prody.py
   └── tools
       ├── data_on_structure.R
       ├── pdb_delhetatm.py
       ├── pdb_delhetatm_ions.py
       ├── pdb_delresname.py
       ├── run_pdb_del_on_directory.sh
       ├── sort_pdbs_from_dataframe.py
       ├── sort_pdbs_from_file.py
       ├── sort_pdbs_has_gap_in_seq.py
       ├── sort_pdbs_has_ligand.R
       ├── split_pdbs_bio3d.R
       ├── subset_alignment_gap_fraction_drop.py
       └── subset_alignment_has_gap_in_seq.sh
```

Here is an example of how to execute the whole pipeline from the command line.

Go to your desired output directory:

    cd <your/output/directory>

Note that execution from inside your output directory is not required, but otherwise you need to specify explicitly the 
path to your input and output directories in each command.

Don't forget to activate the conda environment:

    conda activate EnsembleFlex

Now you can run the scripts as in the example below!  
Just make sure that you use the correct paths in your system  
a) to the scripts: You can run the scripts from your output directory by specifying the correct path to the EnsembleFlex scripts 
(replace `~/path/to` with the correct path on your system).  
b) to your input and output directory:
The following commands assume that your input structure files are located in a folder called `pdbs` and you want to 
create an output directory called `EnsembleFlex` besides your input folder where all outputs will be located. 

- [optional] Tool: PDB splitting on whole directory (if multiple chains are present in your PDB files)

      Rscript ~/path/to/EnsembleFlex/src/tools/split_pdbs_bio3d.R -i pdbs -o split_pdbs

- [optional] Tool: subset PDBs based on gap occurrence (if you want to analyse "structures_with_gaps" and 
"structures_without_gaps" separately)

      python3 ~/path/to/EnsembleFlex/src/tools/sort_pdbs_has_gap_in_pdb_seq.py -i pdbs/ -o .

- Superimpositioning  

      Rscript ~/path/to/EnsembleFlex/src/superimpose_bio3d.R -i pdbs -o EnsembleFlex

      [OR] python3 ~/path/to/EnsembleFlex/src/superimpose_prody.py -i pdbs/ -o EnsembleFlex/

- Analysis
    
      Rscript ~/path/to/EnsembleFlex/src/analyse_flex_bio3d.R -i EnsembleFlex/superimposed -o EnsembleFlex/Analysis_Bio3D -n 3
        
      [optional] python3 ~/path/to/EnsembleFlex/src/analyse_flex_prody.py -i EnsembleFlex/superimposed -o EnsemblFlex/Analysis_ProDy
      
      python3 ~/path/to/EnsembleFlex/src/analyse_flex_sasa_biopython.py -i EnsembleFlex/superimposed/ -o EnsemblFlex/Analysis_SASA_Biopython

      Rscript ~/path/to/EnsembleFlex/src/tools/data_on_structure.R -i EnsembleFlex/superimposed -o EnsemblFlex/Analysis_SASA_Biopython -d EnsembleFlex/Analysis_SASA_Biopython/SASA_global.csv -c sd
    
- Subset clusters based on dataframe (where option `-c ` specifies the column name of the data frame)

      python3 ~/path/to/EnsembleFlex/src/tools/pdb_sorter_from_dataframe.py -i EnsembleFlex/superimposed -o EnsembleFlex/Analysis_Bio3D/Consensus_Clusters -d EnsembleFlex/Analysis_Bio3D/cluster_attributions_with_consensus.csv -c Consensus_Cluster

**Binding site**
- Checking for ligands

      Rscript ~/path/to/EnsembleFlex/src/tools/pdb_sorter_has_ligand.R -i EnsembleFlex/superimposed -o EnsemblFlex

- [optional] Removing crystallization factor DMSO (using resname DMS)

      ~/path/to/EnsembleFlex/src/tools/run_pdb_del_on_directory.sh ~/path/to/EnsembleFlex/src/tools/pdb_delresname.py EnsembleFlex/superimposed DMS [EnsembleFlex/superimposed_no_DMS]

- [optional] Removing Ions

      ~/path/to/EnsembleFlex/src/tools/run_pdb_del_on_directory.sh ~/path/to/EnsembleFlex/src/tools/pdb_delhetatm_ions.py EnsembleFlex/superimposed[_no_DMS] ions [EnsembleFlex/superimposed[_no_DMS]_no_ions]

- [if one of the optional "removing" steps has been performed] Checking again for ligands

      Rscript ~/path/to/EnsembleFlex/src/tools/pdb_sorter_has_ligand.R -i EnsembleFlex/superimposed_no_DMS_no_ions -o EnsembleFlex

- Binding site Identification and Analysis

      Rscript ~/path/to/EnsembleFlex/src/identify_binding_site_bio3d.R -i EnsembleFlex/structures_with_ligand -o EnsembleFlex/BindingSite_ident_Bio3D -d 4.0[3.5]

- [optional] Superimpose only on binding site

      Rscript ~/path/to/EnsembleFlex/src/superimpose_binding_site_bio3d.R -i EnsembleFlex/superimposed -o EnsembleFlex/superimposed_only_on_BindingSite -b EnsembleFlex/BindingSite_ident_Bio3D/binding_site_residue_numbers.txt

- Binding site Flexibility Analysis

      Rscript ~/path/to/EnsembleFlex/src/analyse_flex_binding_site_bio3d.R -i EnsembleFlex/structures_with_ligand -o EnsembleFlex/BindingSite_analysis_Bio3D -b EnsembleFlex/BindingSite_ident_Bio3D/binding_site_residue_numbers.txt

      python3 ~/path/to/EnsembleFlex/src/predict_flex_binding_site_essa_prody.py -i EnsembleFlex/superimposed/ -o EnsembleFlex/Prediction_BindingSite_ESSA_ProDy -b EnsembleFlex/BindingSite_ident_Bio3D/binding_site_residue_numbers.txt


**Flexibility Prediction**

- NMA

      Rscript ~/path/to/EnsembleFlex/src/predict_flex_nma_bio3d.R -i EnsemblFlex/superimposed -o EnsemblFlex/Prediction_NMA_Bio3D -e

      python3 ~/path/to/EnsembleFlex/src/predict_flex_nma_prody.py -i EnsemblFlex/superimposed -o EnsemblFlex/Prediction_NMA_ProDy

- Essential Site Scanning Analysis (ESSA)

      python3 ~/path/to/EnsembleFlex/src/predict_flex_binding_site_essa_prody.py -i EnsemblFlex/superimposed -o EnsemblFlex/Prediction_BindingSite_ESSA_ProDy


### Troubleshooting

##### Problems with finding the right path to Rscript  
If you have R already installed on your system the execution of Rscript may default back to your previous R 
installation. In such case just replace `Rscript` by `${CONDA_PREFIX}/bin/Rscript` in all commands. 
This will ensure usage of the R installation inside the conda environment.

##### Problems with package versions and dependencies (Python an R)
In case you have errors during installation indicating incompatibilities of Python or R packages, try first using the 
provided `conda-lock.yml` file (as indicated in installation instructions, for macOS or Linux), and if this does not 
work on your system (e.g. for Windows) try using the `environment_versioned.yml` file instead of `environment.yml`. 
This will potentially result in slightly older package versions, but they are compatible.


## Output structure explained
*Several types of output files will be generated by the analysis tools. This includes image files (png), 
text files (txt), data tables (csv/tsv), structure files (pdb), pymol scripts (pml).  
Most output will be visualised directly in the browser-based user interface when analysis is performed. 
All generated output will be saved in the respective output directory with potential further information. 
The output directory will contain the following structure of subdirectories, where 
uppercase directory names contain analysis output files of diverse format and 
lowercase directory names only contain .pdb structure files:*

```
EnsembleFlex_output_folder (defined by you)  
├── superimposed  
├── Analysis_Bio3D  
│   └── pymol_pdbs  
├── Analysis_SASA_Biopython  
├── superimposed_no_ions (optional)  
├── superimposed_no_ions_<MOL-ID>  (optional)  
├── structures_with_ligand  
├── structures_without_ligand  
├── BindingSite_ident_Bio3D  
│   └── structures_labeled_binding_site  
├── BindingSite_analysis_Bio3D  
├── Prediction_NMA_Bio3D  
├── Prediction_NMA_ProDy  
└── Prediction_ESSA_ProDy  
```


## Output interpretation help

**1. Flexibility analysis: How to interpret RMSD, RMSF, PCA & UMAP analysis?**

### RMSD & RMSF

    ### RMSD and RMSF analysis for flexibility analysis

    RMSD (Root Mean Square Deviation) and RMSF (Root Mean Square Fluctuation) are widely used metrics for analyzing 
    the flexibility and dynamics of protein structures, particularly in the context of structural ensembles 
    obtained from techniques like NMR spectroscopy, molecular dynamics (MD) simulations, or other computational 
    methods.

    #### RMSD

    RMSD measures the overall structural deviation between two conformations of a protein. It quantifies the 
    average distance between the atoms of a target structure and a reference structure after optimal 
    superimposition. A higher RMSD value indicates a greater structural difference between the two conformations.

    In the context of protein ensembles, RMSD is often used to:

    1. Assess the structural heterogeneity within the ensemble by calculating pairwise RMSDs between all 
    conformations[4].
    2. Monitor the structural deviation of a protein over time in MD simulations by calculating the RMSD from the 
    starting structure at each time step[2].
    3. Evaluate the convergence of NMR structure calculations by examining the RMSD between the calculated 
    models[1].

    #### RMSF

    RMSF, on the other hand, measures the fluctuation of individual atoms or residues around their average 
    positions within an ensemble or trajectory. It provides a residue-level measure of flexibility, highlighting 
    the most mobile regions of the protein.

    RMSF is commonly used to:

    1. Identify flexible loops, hinges, or other dynamic regions in proteins by analysing the RMSF profile along 
    the protein sequence[1][3].
    2. Compare the flexibility patterns between different ensembles, such as NMR models, MD trajectories, or 
    crystal structures[3].
    3. Relate the observed flexibility to functional mechanisms, such as conformational changes or binding 
    events[1][3].

    It's important to note that RMSD and RMSF provide complementary information. While RMSD captures global 
    structural changes, RMSF highlights local flexibility. Additionally, the interpretation of these metrics can 
    be influenced by factors like the alignment strategy, reference structure, and ensemble size[2][4].

    In summary, RMSD and RMSF are invaluable tools for analyzing protein flexibility and dynamics from structural 
    ensembles, enabling researchers to gain insights into conformational heterogeneity, identify mobile regions, 
    and relate structural dynamics to protein function[1][2][3][4].

    Citations:  
    [1] https://europepmc.org/article/MED/24735558  
    [2] https://europepmc.org/article/MED/25816325  
    [3] https://europepmc.org/article/MED/33803249  
    [4] https://europepmc.org/article/MED/20197040  

### PCA & UMAP

    ### Dimension reduction techniques (PCA and UMAP) for flexibility analysis
    
    Dimensionality reduction techniques like Principal Component Analysis (PCA) and Uniform Manifold 
    Approximation and Projection (UMAP) are valuable tools for visualizing conformational changes in protein 
    structure ensembles.
    
    #### PCA for Visualizing Protein Conformations
    
    Principal Component Analysis (PCA) is a linear technique that projects the high-dimensional conformational 
    data (e.g., Cartesian coordinates of atoms) onto a lower-dimensional subspace while preserving the maximum 
    variance in the original data. The key steps are:  
    1.	Calculate the covariance matrix of the conformational ensemble.
    2.	Compute the eigenvectors and eigenvalues of the covariance matrix.
    3.	Project the conformations onto the eigenvectors corresponding to the largest eigenvalues (principal 
    components).[1][2]
    
    In the context of protein conformations, PCA can be applied to the Cartesian coordinates of atoms, or other 
    higher dimensional conformational data, such as distance matrices, or torsion angles, to obtain a reduced 
    representation that captures the dominant modes of motion by projecting the data onto the principal 
    components that capture the largest variance. PCA facilitates the comparison of conformational ensembles 
    from different conditions or mutations by projecting them onto a common subspace defined by the principal 
    components, enabling the identification of differences in conformational sampling.[1][2]  
    The advantage of PCA is that it provides a straightforward way to visualize the conformational space 
    sampled by the ensemble by projecting the high-dimensional conformations onto the first two or three 
    principal components. However, PCA has limitations in accurately representing the local structure and 
    nonlinear relationships in the data, as it assumes linearity.
    
    #### UMAP for Visualizing Protein Conformations
    
    Uniform Manifold Approximation and Projection (UMAP) is a nonlinear dimensionality reduction technique that 
    aims to preserve both local and global structure of the high-dimensional data in the lower-dimensional 
    representation. UMAP works by:
    1.	Constructing a fuzzy topological representation of the high-dimensional data.
    2.	Optimizing a low-dimensional representation that preserves the fuzzy topological structure.  
    
    UMAP can capture complex, nonlinear relationships in the data, making it suited for visualizing the 
    conformational space of proteins, which often involves nonlinear motions like bond rotations and global 
    structural changes. Studies have shown that UMAP can outperform PCA and other linear techniques in 
    preserving the structure of protein conformations in the reduced space, as evidenced by higher Pearson 
    correlation coefficients between high- and low-dimensional distances.
    UMAP's ability to accurately represent the conformational space can aid in understanding free energy 
    barriers, conformational transitions, and other biologically relevant information.
    
    In summary, dimensionality reduction techniques like PCA and UMAP are powerful tools for visualizing and 
    understanding the conformational changes and dynamics of protein structure ensembles. While PCA offers a 
    linear projection and therefore enables a very straightforward molecular visualisation including atom 
    contributions to the principal components (when performed on coordinates), UMAP's nonlinear approach may 
    better capture the complex motions involved in protein conformational changes, making it a more accurate 
    representation of the conformational space.[3][4]
    
    Citations:  
    [1] https://europepmc.org/article/MED/24061923  
    [2] https://doi.org/10.3390/j5020021  
    [3] https://europepmc.org/article/MED/33973773  
    [4] https://europepmc.org/article/MED/27179343  
    [5] https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/Dimension%20Reduction.pdf  


**2. How to interpret analysis performed on two different scales: Backbone vs. All-atom**

    ### Analysis performed on two different scales: Backbone vs. All-atom
    
    The values of RMSD (Root Mean Square Deviation) and RMSF (Root Mean Square Fluctuation) can vary significantly 
    depending on which atoms are included in the computation, such as considering only the backbone atoms versus all 
    heavy atoms. The same holds true for principal component analysis (PCA). The choice of atoms included in the PCA 
    of protein structures can significantly impact the resulting principal components and the interpretation of 
    conformational changes. Here's how these metrics are affected:
    
    RMSD:
    - When calculating RMSD between two protein structures, including only the backbone atoms (e.g., Cα, N, C, O) 
    generally results in lower RMSD values compared to including all heavy atoms[1].
    - Backbone atoms tend to be more structurally conserved, while side chains can exhibit greater conformational 
    variability, leading to higher RMSD values when included[1].
    - The difference in RMSD values between backbone-only and all-atom calculations can be substantial for proteins 
    with flexible loops, termini, or side chains undergoing conformational changes[1].
    - For highly rigid proteins or small conformational changes, the difference may be less pronounced[1].
    
    RMSF:
    - When calculating RMSF, including only backbone atoms will capture the fluctuations of the protein's main chain, 
    but may underestimate the mobility of side chains.
    - Including all heavy atoms in the RMSF calculation provides a more comprehensive picture of residue-level 
    flexibility, as it accounts for side chain motions[3].
    - The RMSF profile along the protein sequence can differ significantly between backbone-only and all-atom 
    calculations, especially for regions with flexible side chains or residues involved in conformational changes[3].
    - Backbone-only RMSF may be more suitable for identifying flexible loops or hinges, while all-atom RMSF can reveal 
    side chain dynamics relevant for binding or functional mechanisms.
    
    PCA:
    - PCA performed on just the backbone atoms (N, Cα, C, O) captures the large-scale, global motions of the protein 
    backbone.
    - It highlights the dominant conformational changes in the overall protein fold and secondary structure elements.
    - Backbone-only PCA may miss local conformational changes involving side-chain motions or loop rearrangements.
    - Including all non-hydrogen atoms (backbone and side-chains) in the PCA captures both global and local 
    conformational changes.
    - It can reveal side-chain motions, loop rearrangements, and subtle local structural variations that may be 
    functionally relevant.
    - The principal components will reflect a combination of backbone and side-chain motions.
    - However, the increased number of atoms can make the PCA computationally more expensive and may introduce noise 
    from flexible side-chains.
    
    In general, including only backbone atoms in RMSD and RMSF calculations can provide a more conservative estimate of 
    structural differences and flexibility, focusing on the core protein fold[5]. Including all heavy atoms captures 
    additional conformational variability from side chains, which can be important for understanding functional 
    dynamics or comparing structures with significant side chain rearrangements[5][6]. 
    Backbone-only PCA captures the dominant global motions but may miss local conformational changes. Including all 
    heavy atoms in the PCA provides a more comprehensive representation of the conformational variability, but it may 
    also introduce noise from highly flexible regions. 
    The choice depends on the specific research question and the level of detail required in the analysis. If the focus 
    is on large-scale conformational changes or global protein motions, backbone-only RMSD/RMSF/PCA may be sufficient. 
    However, if local structural variations or side-chain motions are of interest, including all heavy atoms in the 
    analysis is recommended.
    It's also possible to combine both approaches, using backbone-only PCA to identify the dominant global motions and 
    then analysing the side-chain conformations within the subspace defined by the principal components.[4]
    
    Citations:  
    [1] https://europepmc.org/article/MED/22323224  
    [2] https://europepmc.org/article/MED/28267328  
    [3] https://europepmc.org/article/MED/25816325  
    [4] https://en.wikipedia.org/wiki/Root_mean_square_deviation_of_atomic_positions  
    [5] https://europepmc.org/article/MED/11420442  
    [6] https://europepmc.org/article/MED/11430756  
    [7] http://thegrantlab.org/bio3d_v2/tutorials/principal-component-analysis  


**3. Flexibility prediction: How to interpret Normal Mode Analysis (NMA) results?**

    Normal Mode Analysis (NMA) based on elastic network models like Anisotropic Network Model (ANM) and Gaussian 
    Network Model (GNM) is a widely used computational technique for predicting the intrinsic flexibility of protein 
    structures. These coarse-grained models represent proteins as a network of nodes (atoms) connected by springs 
    (interactions), allowing for efficient calculations of low-frequency normal modes that describe the large-scale 
    collective motions.[1][2][3]
    
    ## Elastic Network Model NMA
    
    - ANM and GNM treat proteins as elastic networks, capturing the overall shape and topology rather than atomic details.
    - They are computationally inexpensive and provide an approximation of the slowest modes that dominate functional motions.
    - The modes describe the directions and relative amplitudes of residue fluctuations around an equilibrium state.
    - ANM is anisotropic, considering directional preferences, while GNM is isotropic, treating all directions equivalently.
    - These modes can reveal mechanisms of conformational changes, functional motions, and allosteric communication pathways.[1][2][3]
    
    ### Differences between ANM and GNM 

    The Anisotropic Network Model (ANM) and the Gaussian Network Model (GNM) differ in their treatment of protein 
    flexibility prediction as follows:
    
    1. Directionality of motions:
     - ANM is anisotropic, meaning it can predict the directionality of residue fluctuations and provide information on 
      the individual x, y, z components of the deformation vectors associated with each mode.[3][4]
     - GNM is isotropic, predicting only the magnitudes of the fluctuation vectors but not their directions. It cannot 
      resolve the individual x, y, z components.[3][5]

    2. Accuracy of fluctuation magnitudes:
     - GNM generally provides a more accurate prediction of the magnitudes of residue fluctuations compared to 
      experimental B-factors or mean-square fluctuations.[3][4]
     - The fluctuation magnitudes predicted by GNM tend to correlate better with experimental data than those from ANM.[3][4]

    3. Computational cost:
     - GNM operates in an N-dimensional configurational space, where N is the number of residues, making it 
      computationally more efficient than ANM.[4]
     - ANM operates in a 3N-dimensional space, considering all 3N degrees of freedom, resulting in a higher 
      computational cost.[4]
    
    In summary, GNM is preferred when accurately predicting the distribution and magnitudes of residue fluctuations is 
    the primary goal, while ANM is better suited for assessing the directionality and mechanisms of motions, albeit 
    with a higher computational expense.[3][4] The choice depends on whether the focus is on the amplitudes or the 
    directions of the motions.

    ## All-Atom NMA
    
    In contrast, all-atom NMA considers the detailed atomic interactions within the protein structure, providing a more 
    accurate but computationally expensive description of the potential energy surface and corresponding normal modes.[1]
    
    - All-atom force fields explicitly model bonded and non-bonded atomic interactions.
    - This level of detail captures finer motions and higher-frequency modes inaccessible to coarse-grained models.
    - All-atom NMA is better suited for studying local flexibility, side-chain motions, and subtle conformational changes.
    - However, the high computational cost limits its application to smaller systems or lower modes.[1][2]
    
    Both coarse-grained and all-atom NMA offer valuable insights into protein dynamics, with the choice depending on 
    the system size, desired accuracy, and the nature of the motions under investigation.[1][2]
    
    Citations:  
    [1] https://europepmc.org/article/MED/33349977  
    [2] https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0427-6  
    [3] https://europepmc.org/article/MED/31510014  
    [4] https://www.taylorfrancis.com/books/edit/10.1201/9781420035070/normal-mode-analysis-ivet-bahar-qiang-cui?refId=7b1a66b2-71fe-4a55-b39a-1212bd8da7b5&context=ubx  
    [4] https://www.ccbb.pitt.edu/archive/Faculty/bahar/b14.pdf  
    [5] https://link.springer.com/article/10.1007/s42493-024-00097-8  
    [6] https://www.sciencedirect.com/science/article/pii/S0166128008005435  


