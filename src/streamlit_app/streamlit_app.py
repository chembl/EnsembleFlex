
import sys
import os
import subprocess
import json

import pandas as pd
import streamlit as st
from stoc import stoc
#from directorypicker import st_directory_picker
from pathlib import Path
import py3Dmol
from stmol import showmol
import Bio.PDB # to get b-factor values for display
import tkinter as tk
from tkinter import filedialog


# make app accessible only through localhost - Not working!
# st.set_option('browser.serverAddress', 'localhost')

filepath = Path(__file__).parent.resolve()
parentfilepath = Path(__file__).parent.parent.resolve()
basefilepath = Path(__file__).parent.parent.parent.resolve()

## Personalisation

st.set_page_config(
    page_title="EnsembleFlex", # the title of the page
    page_icon='➰', # the page’s icon
    # page_icon=str(filepath)+"/EnsembleFlex_icon.png", # the page’s icon
    layout="centered", # ["centered", "wide"] # set Streamlit app to use wide mode
    initial_sidebar_state="auto", # whether sidebar will be initially loaded
)

# Table of contents in sidebar
toc = stoc()

# use custom css styles
# with open(str(filepath) + '/wave.css') as f:
#     st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
# use css file for styling
# with open('./bootstrap_minty.css') as f:
#     css = f.read()
# st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)
#st.markdown('<link rel="stylesheet" type="text/css" href="https://www.example.com/style.css">', unsafe_allow_html=True)

# Function to detect if running in Docker
def is_running_in_docker():
    """Check if the script is running inside a Docker container."""
    return os.path.exists("/.dockerenv")

# Function to list directories for browser-based selection
def list_directories(base_path, show_hidden=False):
    """
    List directories with enhanced filtering and error handling.
    Args:
        base_path (str): Base directory to list subdirectories from
        show_hidden (bool): Whether to show hidden directories
    Returns:
        list: Filtered list of directories
    """
    try:
        # Use pathlib.Path for robust path handling
        path = Path(base_path)
        # Filter directories with additional checks
        directories = [
            d.name for d in path.iterdir()
            if d.is_dir() and
               (show_hidden or not d.name.startswith('.'))
        ]
        return sorted(directories)
    except PermissionError:
        st.warning(f"Permission denied: Cannot access {base_path}")
        return []
    except FileNotFoundError:
        st.warning(f"Directory not found: {base_path}")
        return []
    except Exception as e:
        st.error(f"Unexpected error listing directories: {e}")
        return []


# Function to select a folder using the appropriate method
# Use system-based folder selection for non-Docker environments
def select_folder_non_docker(folder_key):
    """Select a folder using a system-native method (Windows, macOS, Linux)."""
    folder_path = None

    if sys.platform == "darwin":  # macOS (AppleScript)
        script = """
        tell application "System Events"
            activate
            set folderPath to choose folder with prompt "Please select a folder"
            set folderPath to POSIX path of folderPath
        end tell
        return folderPath
        """
        try:
            folder_path = subprocess.check_output(["osascript", "-e", script]).decode("utf-8").strip()
        except subprocess.CalledProcessError:
            folder_path = None

    elif sys.platform.startswith("linux"):  # Linux (Tkinter)
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes('-topmost', 1)
        folder_path = filedialog.askdirectory()
        root.destroy()

    elif sys.platform.startswith("win"):  # Windows (Win32 API)
        try:
            import ctypes
            from ctypes import wintypes
            from win32com.shell import shell, shellcon

            def win32_folder_picker():
                pidl, display_name, image_list = shell.SHBrowseForFolder(
                    None, None, "Select a folder", shellcon.BIF_RETURNONLYFSDIRS
                )
                if pidl:
                    return shell.SHGetPathFromIDList(pidl)
                return None

            folder_path = win32_folder_picker()
        except ImportError:  # Fallback to Tkinter if PyWin32 is unavailable
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            root.wm_attributes('-topmost', 1)
            folder_path = filedialog.askdirectory()
            root.destroy()

    if folder_path:
        st.session_state[folder_key] = folder_path
        st.success(f"📁 {folder_key.replace('_', ' ').title()} set to: `{folder_path}`")

# Function to select a folder using the appropriate method
# Use browser-based folder selection for Docker environments
# **Function for Docker Folder Selection (Now Uses Callbacks for Partial Rerun)**
def update_navigation(folder_key):
    """Callback to update navigation without full script rerun."""
    selected_dir = st.session_state[f"dir_select_{folder_key}"]
    current_path = st.session_state[f"current_path_{folder_key}"]

    new_path = os.path.dirname(current_path) if selected_dir == ".." else os.path.join(current_path, selected_dir)
    if os.path.isdir(new_path):
        st.session_state[f"current_path_{folder_key}"] = new_path

# Function to check for .pdb files in the input directory
def contains_pdb_files(directory):
    """Check if a directory contains any .pdb files."""
    if not os.path.isdir(directory):
        return False
    return any(f.endswith(".pdb") for f in os.listdir(directory))

def create_new_directory(folder_key):
    """Callback to create new directory without full script rerun."""
    current_path = st.session_state[f"current_path_{folder_key}"]
    new_dir_name = st.session_state[f"new_dir_{folder_key}"]

    if new_dir_name:
        new_dir_path = os.path.join(current_path, new_dir_name)
        try:
            os.makedirs(new_dir_path, exist_ok=True)
            st.session_state[f"current_path_{folder_key}"] = new_dir_path
            st.success(f"📁 Created and selected: `{new_dir_path}`")
        except OSError as e:
            st.error(f"❌ Failed to create directory: {e}")

# **Docker Folder Selection (Now Uses Callbacks)**
def select_folder_in_docker(folder_key, title, allow_creation=False):
    """Streamlit-based folder selection for Docker environments with navigation and optional directory creation."""
    current_path = st.session_state[f"current_path_{folder_key}"]
    directories = list_directories(current_path)

    # st.markdown(f"### 📁 {title}")
    st.write(f"**Current Path:** `{current_path}`")

    # **Navigation Section (Uses Callbacks)**
    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        st.write("📂 Navigate Directories:")
    with col2:
        st.selectbox("", [".."] + directories, key=f"dir_select_{folder_key}", label_visibility="collapsed",
                     on_change=update_navigation, args=(folder_key,))
    with col3:
        st.button("🔄 Navigate", key=f"navigate_{folder_key}", on_click=update_navigation, args=(folder_key,))

    # **Folder Creation (Optional, Uses Callback)**
    if allow_creation:
        col1, col2, col3 = st.columns([1, 2, 1])
        with col1:
            st.write("🆕 New Directory Name:")
        with col2:
            st.text_input("", placeholder="Enter directory name", key=f"new_dir_{folder_key}", label_visibility="collapsed")
        with col3:
            st.button("📂 Create Directory", key=f"create_{folder_key}", on_click=create_new_directory, args=(folder_key,))

    # **Final Selection (Only Saves Confirmed Selection)**
    col1, col2 = st.columns([1, 3])
    with col1:
        st.write("☑️ Confirm Selection:")
    with col2:
        if st.button("Set Directory", key=f"set_dir_{folder_key}"):
            selected_path = st.session_state[f"current_path_{folder_key}"]  # Uses latest path
            if os.access(selected_path, os.W_OK):
                st.session_state[folder_key] = selected_path  # Saves only when confirmed
                st.success(f"✅ Selected directory: `{selected_path}`")
            else:
                st.warning("⚠️ Selected directory is not writable.")


def b_fac_on_structure_vis(pdbfilepath):
    '''
    Visualize PDB structure with b-factor coloring (with py3Dmol and stmol)
    :param pdbfilepath: path to PDB file
    :return: visualisation
    '''
    # get b-factor min and max values using biopython
    p = Bio.PDB.PDBParser()
    structure = p.get_structure('', pdbfilepath)
    bfacs = [a.get_bfactor() for a in structure.get_atoms()] #a.get_name()=='CA'
    minbfacs = min(bfacs)
    maxbfacs = max(bfacs)
    print(minbfacs, maxbfacs)
    # visualisation with py3Dmol
    view = py3Dmol.view()
    view.addModel(open(pdbfilepath, 'r').read(), 'pdb')
    view.setBackgroundColor('white')
    # view.setStyle({'cartoon': {'color': 'lightgrey'}})
    view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': maxbfacs, 'max': minbfacs}}})
    #, 'min': min({'prop': 'b'}), 'max': max({'prop': 'b'})  # 'rwb' 'roygb' 'sinebow'
    view.addSurface(py3Dmol.VDW, {'opacity': 0.6, 'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': maxbfacs, 'max': minbfacs}}, #, 'min': 0
                    {"hetflag": False})
    view.addStyle({"elem": "C", "hetflag": True},
                  {"stick": {"color": "lightgrey", "radius": 0.2}})
    view.addStyle({"hetflag": True},
                  {"stick": {"radius": 0.2}})
    view.setHoverable({}, True, '''function(atom,viewer,event,container) {
                       if(!atom.label) {
                        atom.label = viewer.addLabel(atom.resn+atom.resi+":"+atom.b,
                        {position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                       }}''',
                   '''function(atom,viewer) { 
                       if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                       }
                    }''') #+":"+atom.atom
    view.zoomTo()
    # final display with stmol
    showmol(view)

def multimodel_animation(pdbfilepath):
    '''
    Visualize animated multimodel PDB structure (with py3Dmol and stmol)
    :param pdbfilepath: path to multimodel PDB file
    :return: visualisation
    '''
    view = py3Dmol.view()
    # view.setBackgroundColor('white')
    with open(pdbfilepath) as ifile:
        pdbsystem = "".join([x for x in ifile])
    view.addModelsAsFrames(pdbsystem)
    # view.setStyle({'model': -1},{'cartoon': {'color': 'lightgrey'}})
    view.setStyle({'sphere':{'color': 'lightgrey', "radius": 0.5}})
    view.addSurface(py3Dmol.VDW, {'opacity': 0.6, 'color': 'lightgrey'},
                    {"hetflag": False})
    view.addStyle({"elem": "C", "hetflag": True},
                  {"stick": {"color": "lightgrey", "radius": 0.2}})
    view.addStyle({"hetflag": True},
                  {"stick": {"radius": 0.2}})
    view.zoomTo()
    view.animate({'loop': "backAndForth"}) #forward
    showmol(view)


# Session states initialisation
if "input_directory" not in st.session_state:
    st.session_state.input_directory = "/app/docker-data" if is_running_in_docker() else ""
    input_directory = st.session_state.input_directory
if "output_directory" not in st.session_state:
    st.session_state.output_directory = "/app/docker-data" if is_running_in_docker() else ""
    output_directory = st.session_state.output_directory
if 'input_directory_liganded' not in st.session_state:
    st.session_state.input_directory_liganded = st.session_state.output_directory
    input_directory_liganded = st.session_state.input_directory_liganded
# Initialize navigation state
if "current_path_input_directory" not in st.session_state:
    st.session_state.current_path_input_directory = st.session_state.input_directory
if "current_path_output_directory" not in st.session_state:
    st.session_state.current_path_output_directory = st.session_state.output_directory
if "current_path_input_directory_liganded" not in st.session_state:
    st.session_state.current_path_input_directory_liganded = st.session_state.output_directory

if 'superimposed' not in st.session_state:
    st.session_state.superimposed = ""
if 'superimposeisdone' not in st.session_state:
    st.session_state.superimposeisdone = False
if 'superimposeBio3Disdone' not in st.session_state:
    st.session_state.superimposeBio3Disdone = False
if 'Bio3Danalysisdone' not in st.session_state:
    st.session_state.Bio3Danalysisdone = False
if 'SASAanalysisdone' not in st.session_state:
    st.session_state.SASAanalysisdone = False
if 'PDBhasLigandSortisdone' not in st.session_state:
    st.session_state.PDBhasLigandSortisdone = False
if 'BSidentifyisdone' not in st.session_state:
    st.session_state.BSidentifyisdone = False
if 'loadinitialresidues' not in st.session_state:
    st.session_state.loadinitialresidues = False
if 'superimposed_bs' not in st.session_state:
    st.session_state.superimposed_bs = ""
if 'BSanalysisdone' not in st.session_state:
    st.session_state.BSanalysisdone = False
if 'residue_numbers_file' not in st.session_state:
    st.session_state.residue_numbers_file = ""
if 'wateranalysisdone' not in st.session_state:
    st.session_state.wateranalysisdone = False
if 'NMABio3Disdone' not in st.session_state:
    st.session_state.NMABio3Disdone = False
if 'NMAProDyisdone' not in st.session_state:
    st.session_state.NMAProDyisdone = False
if 'aaNMABio3Disdone' not in st.session_state:
    st.session_state.aaNMABio3Disdone = False
if 'ESSAisdone' not in st.session_state:
    st.session_state.ESSAisdone = False

# def main():

#st.title("EnsembleFlex - Flexibility Analysis of Structure Ensembles")
toc.title("EnsembleFlex - Flexibility Analysis of Structure Ensembles")
st.markdown('''
_An interactive tool for analysing structure ensembles._
''')

toc.header("Introduction")
st.markdown(":rainbow[Welcome to an exploration of the dynamics within protein structure ensembles!]")
st.markdown('''This interface serves as exploration gateway into the fascinating field of protein structural flexibility 
    in an ***automated*** manner, where tools are ***streamlined*** and analysis is easily ***reproducible*** and 
    ***fast***.''')
with st.expander("**Why investigating protein flexibility?**"):
    st.markdown('''
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
    ''')

toc.subheader("Methods explained")
st.markdown("Please read the following dropdown explanations when using for the first time:")
with st.expander("**1. Flexibility analysis: How to interpret RMSD, RMSF, PCA & UMAP analysis?**"):
    tab1, tab2 = st.tabs(["RMSD & RMSF", "PCA & UMAP"])
    with tab1:
        st.markdown('''
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
        ''')
# with st.expander("Methods explained: PCA, UMAP"):
    with tab2:
        st.markdown('''
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
        
        ''')

with st.expander("**2. How to interpret analysis performed on two different scales: Backbone vs. All-atom**"):
    st.markdown('''
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
    - For highly rigid proteins or small conformational changes, the difference may be less pronounced[1][2].
    
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
    then analysing the side-chain conformations within the subspace defined by the principal components.[7]
    
    Citations:  
    [1] https://europepmc.org/article/MED/22323224  
    [2] https://europepmc.org/article/MED/28267328  
    [3] https://europepmc.org/article/MED/25816325  
    [4] https://en.wikipedia.org/wiki/Root_mean_square_deviation_of_atomic_positions  
    [5] https://europepmc.org/article/MED/11420442  
    [6] https://europepmc.org/article/MED/11430756  
    [7] http://thegrantlab.org/bio3d_v2/tutorials/principal-component-analysis  

    ''')

with st.expander("**3. Flexibility prediction: How to interpret Normal Mode Analysis (NMA) results?**"):
    st.markdown('''

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

    ''')

toc.subheader("Workflow overview")
st.image(str(basefilepath) + '/documentation/EnsembleFlex_Workflow.png') #, caption='EnsembleFlex_Workflow'
st.markdown('''
    The workflow steps are ideally performed in sequential order, as some steps are dependant on previous output, 
    but they can be rerun if upstream parameters have been changed:  
    :one: *Settings* - [mandatory] The place to provide your input and output directory  
        - *Display previous analysis results* - [optional] Use this to display analysis you have already run  
    :two: *Superimpose* - [mandatory confirmation] You can choose to run superpositioning with provided tools or not  
    :three: ***Flexibility Analysis*** - The main analysis computations are performed here  
    :four: *Binding Site Analysis* - [optional] Recommended analysis if you have ligands in your structures  
    :five: *Conserved Water Analysis* - [optional]  
    :six: *Compare with predicted flexibility* - [optional]  
    ''')

toc.subheader('Output structure explained')
st.markdown('*Several types of output files will be generated by the analysis tools. '
            'This includes image files (png), text files (txt), data tables (csv/tsv), structure files (pdb), '
            'pymol scripts (pml).  \n'
            'Most output will be visualised directly in the interface when analysis is performed. '
            'All generated output will be saved in the respective output directory with potential further information. '
            'The output directory will contain the following structure of subdirectories, where '
            'uppercase directory names contain analysis output files of diverse format and '
            'lowercase directory names only contain .pdb structure files:*')
st.code('''
EnsembleFlex_output_folder (defined by you)  
├── superimposed  
├── Analysis_Bio3D  
│   └── pymol_pdbs  
├── Analysis_SASA_Biopython  
├── Analysis_Waters  
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
'''
)

st.divider()  # Draws a horizontal rule
#######################################################################################################################

toc.header(":one: Settings")
st.markdown('''
_Here is the place for data input and output and additional variable settings._
''')
# st.markdown('''### Data loading''')
# structureFiles = st.file_uploader("Select your pdb files", accept_multiple_files=True, type=['pdb','ent','cif','cif.txt'])
# st.write('You selected `%s` files.' % len(structureFiles))
# st.write(structureFiles)

st.markdown('''#### 📁 Input Directory''')
st.write('Please select your input directory (where your structure files are located):')
# input_directory = st_directory_picker_input(key="input_directory")

# Input Folder Picker
if not is_running_in_docker():
    if st.button("Select Input Folder"):
        select_folder_non_docker("input_directory")
else:
    select_folder_in_docker("input_directory", "Select Input Directory")

if "input_directory" in st.session_state:
    st.write(f"📂 **Selected Input Folder:** `{st.session_state.input_directory}`")
    input_directory = st.session_state.input_directory
    # ✅ **Check if .pdb files exist in the input directory**
    if contains_pdb_files(input_directory):
        st.success("✅ PDB files found in the input directory.")
    else:
        st.warning("⚠️ No PDB files found in the selected input directory. Please ensure the correct folder is selected.")


st.markdown('''#### 📁 Output Directory''')
st.write('Please select your output directory (where all calculation results and files will be saved)\n '
         'It is recommended to select an empty folder:')

# Output Folder Picker
if not is_running_in_docker():
    if st.button("Select Output Folder"):
        select_folder_non_docker("output_directory")
else:
    select_folder_in_docker("output_directory", "Select Output Directory", allow_creation=True)

if "output_directory" in st.session_state:
    st.write(f"📂 **Selected Output Folder:** `{st.session_state.output_directory}`")
    output_directory = st.session_state.output_directory
    # **Check if output directory is empty**
    if not len(os.listdir(output_directory)) == 0:
        st.warning("⚠️ The selected output directory is not empty. Existing results may be overwritten.")

# output_directory = st_directory_picker_output(key="output_directory")

if st.session_state.input_directory and st.session_state.output_directory:
    st.success("Input and output directories are set. You're ready to start the analysis.", icon="✅")
else:
    st.error('Please select both input and output directories to proceed.', icon="🚨")
st.divider()
#######################################################################################################################

# Define the path where the session state will be saved
if st.session_state.output_directory:
    session_state_file = str(output_directory)+"/session_state.json"

# List of keys that should be treated as paths (and excluded from saving)
EXCLUDE_PATH_KEYS = [
    "outpath", "inpath_liganded", "output_directory",
    "input_directory_liganded", "input_directory", "inpath",
    "new_dir", "display_previous_btn"
]
PATH_KEYS = [
    "superimposed", "superimposed_bs", "residue_numbers_file"
]

# List of suffixes for keys to be excluded
EXCLUDE_SUFFIXES = ["subdirectory", "left", "right", "btn"]

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Path):
            return str(obj)
        return super().default(obj)

def save_session_state():
    # Convert session state to a dictionary, excluding streamlit's internal keys,
    # keys ending with specified suffixes, and path keys
    session_state_dict = {
        key: value for key, value in st.session_state.items()
        if not key.startswith('_')
        and not any(key.endswith(suffix) for suffix in EXCLUDE_SUFFIXES)
        and key not in EXCLUDE_PATH_KEYS
    }
    # Convert the dictionary to a JSON string using the custom encoder
    json_str = json.dumps(session_state_dict, indent=2, cls=CustomEncoder)
    # Save the JSON string to a file
    with open(session_state_file, 'w') as f:
        f.write(json_str)


# Save button in sidebar
if st.sidebar.button("Save Session State"):
    save_session_state()


def load_session_state():
    if os.path.exists(session_state_file):
        try:
            # Read the JSON content
            with open(session_state_file, 'r') as f:
                content = f.read()

            session_state_dict = json.loads(content)

            # Update the session state
            for key, value in session_state_dict.items():
                # Convert string paths back to Path objects if necessary
                if key in PATH_KEYS and value:
                    st.session_state[key] = Path(value)
                else:
                    st.session_state[key] = value

            st.success("Session state loaded successfully!")
        except json.JSONDecodeError:
            st.error("Invalid JSON file. Please check the session state file.")
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
    else:
        st.warning("No saved session state file found.")


st.markdown('''#### Display previous analysis results''')
st.write('In case you have already performed the analysis using this interface and you just want to re-display the '
         'results, use one of the buttons below (after having set the correct `Output Directory` above!).')
st.write('Be aware that pushing any further downstream button may overwrite your previous results.')
st.write('''There are two options:\n  
A) You have run parts of the pipeline and want to continue from there AND you have previously saved the state of your 
analysis through the `Save Session State` button in the side bar. (The file must be located in your Output Directory.)
''')

if st.button('A) Display previous analysis results - use saved session state', key="display_previous_ss_btn"):
    # Load session state from file
    load_session_state()
    # # Display current session state
    # st.write("Current Session State:", st.session_state)

st.write('''  
B) You have run the whole pipeline and want to visualise all results. (No need for having saved the session state.)  
*Note that if the output files are not found in their destined place (e.g. you have moved them or you haven't run parts 
of the pipeline) then the interface will likely throw an Error at that point, but you can still continue from there and 
(re-)run the missing parts.*
''')

if st.button('B) Display previous analysis results - set all run', key="display_previous_all_btn"):
    st.session_state.superimposed = str(output_directory)+'/superimposed'
    st.session_state.superimposeisdone = True
    st.session_state.superimposeBio3Disdone = True
    st.session_state.Bio3Danalysisdone = True
    st.session_state.SASAanalysisdone = True
    st.session_state.PDBhasLigandSortisdone = True
    st.session_state.BSidentifyisdone = True
    st.session_state.superimposed_bs = ""
    st.session_state.BSanalysisdone = True
    st.session_state.wateranalysisdone = True
    st.session_state.loadinitialresidues = True
    st.session_state.NMABio3Disdone = True
    st.session_state.NMAProDyisdone = True
    st.session_state.aaNMABio3Disdone = True
    st.session_state.ESSAisdone = True
    st.success("Whole pipeline set to run.")

st.divider()
#######################################################################################################################
st.markdown('''#### Optional Preparation Tools''')
st.write('*After running any of the preparation tools, don\'t forget to reset your Input (and Output) '
         'Directory to the correct one and verify that your new Input Directory contains only the structures you want '
         'to include in the analysis.*')

st.write('##### ✂️ Split PDBs (in case they contain multiple chains or are multi-model files)')
st.write('This is an additional optional tool that might help you with preparing you pdb ensemble.  '
         'In case your PDB structure files contain multiple chains or are multi-model files you need to split them '
         'before running the analysis.')
if st.button('Split PDBs', key="split_pdbs_btn"):
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/tools/split_pdbs_bio3d.R",
         '-i', str(input_directory),
         '-o', str(output_directory) + "/split_pdbs"])
    st.write('PDB files are split and saved in: ', str(output_directory) + '/split_pdbs')

st.write('##### ↩️↪️ Sort/subset your pdb files based on gap occurrence')
st.write('This is an additional optional tool that might help you with sorting/subsetting you pdb ensemble.  '
         'In case some of your structures have a missing loop region, you may want to consider subsetting your '
         'ensemble and analyse the two subsets "structures_with_gaps" and "structures_without_gaps" separately '
         '(in a second run), as the missing residues may hide some movement that is present in the structures with '
         'full sequences.')
if st.button('Subset gap structures', key="subset_gap_pdbs_btn"):
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/tools/sort_pdbs_has_gap_in_pdb_seq.py",
         '-i', str(input_directory),
         '-o', str(output_directory)])
    st.write('PDB files are sorted into subfolders "structures_with_gaps" and "structures_without_gaps" in: ', output_directory)

st.divider()
#######################################################################################################################

toc.header(":two: Superimpose")
st.markdown('''
_Superimpose your structures globally with one of the provided methods._
''')
superimp_method = st.radio(
    "Select the method you would like to use",
    ["Bio3D", "ProDy", "None"],
    captions = ["This method is recommended due to its speed and dedicated 'rigid core' detection method.",
                "This is an alternative method (not recommended, slow).",
                "My structures are already superimposed the way I want them."],
    key="superimp_method")
st.write("Your choice is: ", superimp_method)


def superimpose(superimp_method):
    if superimp_method == "Bio3D":
        result = subprocess.run(
            ["Rscript", str(parentfilepath)+"/superimpose_bio3d.R", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("Bio3D calculations are running...")
        # # Show stdout for external command
        # for line in iter(lambda: result.stdout.readline(), b""):
        #     st.text(line.decode("utf-8"))
        st.session_state.superimposed = str(output_directory)+'/superimposed'
        st.session_state.superimposeisdone = True
        st.session_state.superimposeBio3Disdone = True
    if superimp_method == "ProDy":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath)+"/superimpose_prody.py", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("ProDy calculations are running...")
        st.session_state.superimposed = str(output_directory)+'/superimposed/split_ensemble'
        st.session_state.superimposeisdone = True
    if superimp_method == "None":
        st.session_state.superimposed = str(input_directory)
        st.write("Superimposed structures from input directory will be used.")

st.button('Confirm/Go!', key="superimp_btn", on_click=superimpose, args=[superimp_method])

superimposed = st.session_state.superimposed

if st.session_state.superimposeisdone == True:
    with st.container(border=True, height=80):
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed)

if st.session_state.superimposeBio3Disdone == True:
    with st.container(border=True, height=500):
        st.markdown('##### Identified "rigid core" residues that were used for superpositioning \n '
                    '(highlighted on reference structure)')
        b_fac_on_structure_vis(str(output_directory) + '/superimp_core_labelled_on_ref.pdb')

st.divider()
#######################################################################################################################

toc.header(":three: Flexibility Analysis")
st.markdown('''
_Investigate structural flexibility globally using the provided methods._
''')

st.markdown('''##### Variables''')
number_of_groups = st.number_input("Desired number of clusters\n\n"
                                   "(You may want to run it first with the default and adjust the value in subsequent "
                                   "runs based on the output - see 'Overall clustering results'.)", value=3)
st.write("The selection is "+str(number_of_groups)+" clusters.")

# st.subheader("A) Structural variability analysis")
toc.subheader("A) Structural variability analysis")
st.markdown("[Used package/tool: Bio3D, umap (R)]")
st.markdown("Analysis performed on the protein backbone and additionally on 'all-atom' including side chains:\n"
            " - Pairwise Root Mean Square Deviations (RMSD) between structures\n"
            " - Root Mean Square Fluctuations (RMSF) per residue\n"
            " - Difference Torsion/Dihedral analysis between structures\n"
            " - Difference Distance Matrix (DDM) analysis between structures (only all-atom)\n"
            " - Principal Component Analysis (PCA) on coordinates, on torsion data, and on distance matrices\n"
            " - Uniform Manifold Approximation and Projection (UMAP) analysis on coordinates\n"
            " - Clustering based on RMSD values, PCA, and UMAP analysis\n"
            )
            # " - Ensemble Normal Mode Analysis (eNMA)")

outputdirBio3D = str(output_directory) + '/Analysis_Bio3D'

def run_analysis_Bio3D():
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/analyse_flex_bio3d.R", '-i', superimposed, '-o', outputdirBio3D, '-n', str(number_of_groups)])
    st.write("Superimposed structures are taken from ", superimposed)
    st.write("Bio3D calculations are running...")
    st.write("Files are saved in: ", outputdirBio3D)
    st.session_state.Bio3Danalysisdone = True


st.button('Run', key="analyse_Bio3D_btn", on_click=run_analysis_Bio3D)
st.markdown("Results in form of images will be displayed below.  "
            "Please note that more information is available in form of data connected to structures (.pdb .pml) "
            "in the output folder.")

if st.session_state.Bio3Danalysisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Flexibility Analysis Results")
        tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
            "Base info", "RMSF & RMSD analysis", "PCA", "UMAP analysis",
            "All-atom RMSF & RMSD analysis", "All-atom PCA & DDM analysis",
            "Clustering"])
        with tab1:
            st.markdown("#### Base info: Alignment overview and B-factors of reference")
            st.image(outputdirBio3D + '/alignment_overview.png', caption='alignment overview')
            st.image(outputdirBio3D + '/B-factors.png', caption='B-factors of reference structure')
            #st.image("https://static.streamlit.io/examples/cat.jpg", width=200)
        with tab2:
            st.markdown("#### RMSF and RMSD analysis")
            st.markdown("#### - on backbone coordinates")
            st.image(outputdirBio3D + '/RMSF.png', caption='RMSF')
            st.write("Reference structure residues colored by backbone RMSF")
            b_fac_on_structure_vis(outputdirBio3D + '/RMSF_data_on_structure.pdb')
            st.image(outputdirBio3D + '/RMSD_hist.png', caption='RMSD histogram')
            st.image(outputdirBio3D + '/RMSD_heatmap.png', caption='RMSD heatmap')
        with tab3:
            st.markdown("#### Principal Component Analysis (PCA) on coordinates and on torsion data")
            st.markdown("#### - on coordinates (backbone)")
            st.image(outputdirBio3D + '/PCA.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_residue_contribution.png', caption='PCA residue contribution')
            st.image(outputdirBio3D + '/PCA_dendrogram.png', caption='PCA dendrogram')
            st.write("Sampling along PC1")
            multimodel_animation(outputdirBio3D + '/PC1.pdb')
            st.write("Sampling along PC2")
            multimodel_animation(outputdirBio3D + '/PC2.pdb')
            st.write("Sampling along PC3")
            multimodel_animation(outputdirBio3D + '/PC3.pdb')
            st.markdown("#### - on torsion data (backbone)")
            st.image(outputdirBio3D + '/PCA_on_Torsion.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_on_Torsion_loadings.png', caption='PCA loadings')
            st.image(outputdirBio3D + '/PCA_on_Torsion_dendrogram.png', caption='PCA dendrogram')
        with tab4:
            st.markdown("#### Uniform Manifold Approximation and Projection (UMAP) analysis")
            st.markdown("#### - on coordinates (backbone)")
            st.image(outputdirBio3D + '/UMAP.png', caption='2D UMAP plot')
            st.image(outputdirBio3D + '/UMAP_dendrogram.png', caption='2D UMAP dendrogram')
        with tab5:
            st.markdown("#### All-atom RMSF and RMSD analysis")
            st.markdown("#### - on all-atom coordinates (including side chains)")
            st.image(outputdirBio3D + '/RMSF_allatom.png', caption='All-atom RMSF')
            st.write("Reference structure residues colored by 'all-atom' RMSF")
            b_fac_on_structure_vis(outputdirBio3D + '/RMSF_allatom_data_on_structure.pdb')
            st.image(outputdirBio3D + '/RMSD_hist_allatom.png', caption='All-atom RMSD histogram')
            st.image(outputdirBio3D + '/RMSD_heatmap_allatom.png', caption='All-atom RMSD heatmap')
        with tab6:
            st.markdown("#### All-atom PCA")
            st.markdown("#### - on coordinates (all-atom)")
            st.image(outputdirBio3D + '/PCA_allatom.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_loadings_allatom.png', caption='PCA loadings')
            st.image(outputdirBio3D + '/PCA_dendrogram_allatom.png', caption='PCA dendrogram')
            st.write("PC1")
            multimodel_animation(outputdirBio3D + '/PC1_allatom.pdb')
            st.write("PC2")
            multimodel_animation(outputdirBio3D + '/PC2_allatom.pdb')
            st.write("PC3")
            multimodel_animation(outputdirBio3D + '/PC3_allatom.pdb')
            st.markdown("#### - on difference distance matrices (all-atom)")
            try:
                st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix.png', caption='PCA overall info')
                st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix_loadings.png', caption='PCA loadings')
                st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix_dendrogram.png', caption='PCA dendrogram')
                # st.markdown("#### All-atom Uniform Manifold Approximation and Projection (UMAP) analysis")
                # st.markdown("#### - on coordinates (all-atom)")
                # st.image(outputdirBio3D + '/UMAP_allatom.png', caption='2D UMAP plot')
                # st.image(outputdirBio3D + '/UMAP_dendrogram_allatom.png', caption='2D UMAP dendrogram')
            except:
                st.write("ERROR: No output available.")
            st.markdown("#### Difference Distance Matrix (DDM) analysis between structures (only all-atom)")
            try:
                st.image(outputdirBio3D + '/eDDM_complete.png', caption='Ensemble Difference Distance Matrix')
            except:
                st.write("ERROR: No output available.")
        with tab7:
            st.markdown("#### Overall clustering results")
            st.image(outputdirBio3D + '/cluster_attributions_heatmap.png', caption='Cluster attributions')
            st.markdown("#### Coassignment of clustering results")
            st.markdown("With the coassignment of clustering results you can assess the consensus clustering.")
            try:
                st.image(outputdirBio3D + '/coassignment_heatmap.png', caption='Cluster coassignment heatmap')
            except:
                st.write("---")
            st.markdown("#### Clustering validation metrics")
            st.write("""
            INTERNAL VALIDATION MEASURES rely on information in the data only, that is the characteristics of the 
            clusters themselves, such as compactness and separation. Ideally one would want the clusters to be as 
            compact and separated as possible. This can be measured using the following metrics:  
            - CONNECTIVITY:  
            This measure reflects the extent to which items that are placed in the same cluster are also considered 
            their nearest neighbors in the data space - or, in other words, the degree of connectedness of the clusters.
            [It should be minimised.]  
            - DUNN INDEX:  
            Dunn Index represents the ratio of the smallest distance between observations not in the same cluster to 
            the largest intra-cluster distance. [The nominator should be maximised and the denomitor minimised, 
            therefore the index should be maximized.]  
            - SILHOUETTE INDEX:  
            This index defines compactness based on the pairwise distances between all elements in the cluster, and 
            separation based on pairwise distances between all points in the cluster and all points in the closest other 
            cluster (Van Craenendonck & Blockeel 2015) [The values as close to (+) 1 as possible are most desirable.]  
            """)
            st.markdown("##### - on backbone RMSD data")
            with open(outputdirBio3D + '/cluster_validation_RMSD.txt') as f:
                st.write(f.readlines())
            st.markdown("##### - on backbone coordinate PCA results")
            with open(outputdirBio3D + '/cluster_validation_PCA.txt') as f:
                st.write(f.readlines())
            st.markdown("##### - on backbone coordinate UMAP results")
            with open(outputdirBio3D + '/cluster_validation_UMAP.txt') as f:
                st.write(f.readlines())



st.markdown("##### Helper Tool: Sort structures based on clustering")
st.markdown("[Used package/tool:  (Python)]")

cluster_option = st.selectbox(
    'Which clustering result would you like to use to sort structures into subdirectories?',
    ("backbone_RMSD", "backbone_PCA_onCoords", "backbone_UMAP_onCoords", "backbone_PCA_onTorsion",
     "allatom_RMSD", "allatom_PCA_onCoords", "allatom_UMAP_onCoords", "allatom_PCA_onDist")) #"Consensus_Cluster",
st.write('Your choice is:', cluster_option)

if st.button('Run', key="subset_clusters_btn"):
    outputdirBio3D_clusters = str(output_directory) + '/Analysis_Bio3D/' + cluster_option
    cluster_dataframe = str(output_directory) + '/Analysis_Bio3D/cluster_attributions.csv'
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/tools/sort_pdbs_from_dataframe.py", '-i', superimposed, '-o',
         outputdirBio3D_clusters, '-d', cluster_dataframe, '-c', cluster_option])
    st.write("PDB files are copied to: ", outputdirBio3D_clusters)


# toc.subheader("B) using ProDy")
# st.markdown(" - Root Mean Square Fluctuations (RMSF)\n"
#             " - Principal Component Analysis (PCA)\n"
#             " - Anisotropic Network Model (ANM) Normal Mode Analysis (NMA)\n"
#             " - Dynamical Domain Decomposition of reference structure")
#
# def run_analysis_ProDy():
#     outputdirProDy = str(output_directory) + '/Analysis_ProDy'
#     result = subprocess.run(
#         [f"{sys.executable}", str(parentfilepath) + "/analyse_flex_prody.py", '-i', superimposed, '-o', outputdirProDy])
#     st.write("Superimposed structures are taken from ", superimposed)
#     st.write("ProDy calculations are running...")
#     st.write("Files are saved in: ", outputdirProDy)
#     show_report(parentpath=outputdirProDy, filename="/analysis_prody_html_report.html")
#
# st.button('Run', key="analyse_ProDy_btn", on_click=run_analysis_ProDy)


toc.subheader("B) SASA Difference Analysis")
st.markdown("[Used package/tool: Biopython (Python)]")
st.markdown('''
_Investigate differences in Solvent Accessible Surface Area (SASA) of protein residues across the whole ensemble._  
_This highlights residues that experience differences in solvent accessibility across the ensemble._  
''')

outputdirSASA = str(output_directory) + '/Analysis_SASA_Biopython'

sasa_choice = st.radio(
    "Select the method you would like to use",
    ["Fast SASA", "SASA using reference", "SASA using alignment"],
    captions = ["This is recommended due to its speed, but it expects coherent residue numbering.",
                "This uses the reference structure as residue indices (for incoherent residue numbering).",
                "This uses a sequence alignment as residue indices (for incoherent residue numbering) - may be slow."],
    key="sasa_choice")
st.write("Your choice is: ", sasa_choice)

# if st.button('Run fast SASA (expecting coherent residue numbering)', key="analyse_SASA_btn"):
def run_analysis_SASA(sasa_choice):
    if sasa_choice == "Fast SASA":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath) + "/analyse_flex_sasa_biopython.py", '-i', superimposed, '-o', outputdirSASA])
    if sasa_choice == "SASA using reference":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath) + "/analyse_flex_sasa_biopython_using_reference.py", '-i', superimposed, '-o', outputdirSASA])
    if sasa_choice == "SASA using alignment":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath) + "/analyse_flex_sasa_biopython_using_alignment.py", '-i', superimposed, '-o', outputdirSASA])
    st.write("Superimposed structures are taken from ", superimposed)
    st.write("SASA calculations are done.")
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/tools/data_on_structure.R", '-i', superimposed, '-o', outputdirSASA,
         '-d', outputdirSASA+"/SASA_global.csv", '-c', "sd"])
    st.write("SASA standard deviation per residue saved in b-factor on reference structure.")
    st.write("Files are saved in: ", outputdirSASA)
    st.session_state.SASAanalysisdone = True

st.button('Run SASA', key="analyse_SASA_btn", on_click=run_analysis_SASA, args=[sasa_choice])


if st.session_state.SASAanalysisdone == True:
    with st.container(border=True, height=600):
        st.markdown("#### SASA difference analysis result")
        st.write("More data is provided in structural and table format in the output directory: ", outputdirSASA)
        st.image(outputdirSASA + '/SASA_sd.png', caption='SASA standard deviation per residue')
        b_fac_on_structure_vis(outputdirSASA + '/sd_data_on_structure.pdb')

st.divider()
#######################################################################################################################


toc.header(":four: Binding Site Analysis")
#st.header(":four: Binding Site Analysis")
st.markdown('''Important note: Make sure that your structures have consistent residue numbering! (The following tools
are based on residue numbers only, not on alignment position.)
''')

toc.subheader("1 - Identify the Binding Site")
st.markdown("[Used package/tool: Bio3D (R)]")
st.markdown('''
_Please provide a directory of structures with ligands (no apo-structures)._  
_Be aware that any molecule that is not protein, nucleic acid or water will be considered as "ligand"._
''')


st.markdown('''##### Helper Tool: Sort liganded structures and remove unwanted molecules''')
st.markdown("[Used package/tool: Bio3D (R)]")
st.write('This tool automatically sorts copies of your superimposed structures into the subfolders "structures_with_ligand" and "structures_without_ligand".')
st.write('Additionally, you have the option to remove *ions* and other *non-relevant/unwanted* small molecules (such as crystallisation additives or co-factors) before sorting in order to exclude them from the analysis.')

removeIons = st.toggle('Remove all *ions* from all structures.', value=True)
removeLigs = st.toggle('Remove *non-relevant/unwanted* molecules (such as crystallisation additives or co-factors) from all structures.')
if removeLigs:
    ligs_to_remove = st.text_input('Ligand IDs to remove', key='rem_ligs_input', placeholder='DMS,PEG',
                                   help='Please specify the ligand ID codes as used in the pdb files.')

def sortPDBs():
    if not removeIons and not removeLigs:
        result = subprocess.run(
            ['Rscript', str(parentfilepath)+'/tools/sort_pdbs_has_ligand.R', '-i', superimposed,
             '-o', output_directory])
    elif removeIons and not removeLigs:
        result = subprocess.run(
            ['sh', str(parentfilepath)+'/tools/run_pdb_del_on_directory.sh',
             str(parentfilepath)+'/tools/pdb_delhetatm_ions.py',
             superimposed, 'ions', str(output_directory)+'/superimposed_no_ions'])
        result = subprocess.run(
            ['Rscript', str(parentfilepath)+'/tools/sort_pdbs_has_ligand.R',
             '-i', str(output_directory)+'/superimposed_no_ions',
             '-o', output_directory])
        st.write("Ions are removed and structures are saved in: ", str(output_directory)+'/superimposed_no_ions')
    elif removeIons and removeLigs:
        result = subprocess.run(
            ['sh', str(parentfilepath)+'/tools/run_pdb_del_on_directory.sh',
             str(parentfilepath)+'/tools/pdb_delhetatm_ions.py',
             superimposed, 'ions', str(output_directory)+'/superimposed_no_ions'])
        result = subprocess.run(
            ['sh', str(parentfilepath)+'/tools/run_pdb_del_on_directory.sh',
             str(parentfilepath)+'/tools/pdb_delresname.py',
             str(output_directory)+'/superimposed_no_ions', str(ligs_to_remove), str(output_directory)+'/superimposed_no_ions_'+str(ligs_to_remove)])
        result = subprocess.run(
            ['Rscript', str(parentfilepath) + '/tools/sort_pdbs_has_ligand.R',
             '-i', str(output_directory)+'/superimposed_no_ions_'+str(ligs_to_remove),
             '-o', output_directory])
        st.write("Ions and specified molecules are removed and structures are saved in: ", str(output_directory)+'/superimposed_no_ions_'+str(ligs_to_remove))
    liganded = str(output_directory) + '/structures_with_ligand'
    not_liganded = str(output_directory) + '/structures_without_ligand'
    st.write("The subfolders are located in ", output_directory)
    st.session_state.PDBhasLigandSortisdone = True

st.button('Sort structures', key='sort_btn', on_click=sortPDBs)

if st.session_state.PDBhasLigandSortisdone == True:
    with st.container(border=True, height=300):
        st.write('Structures have been sorted into respective subfolders '
                 '"structures_with_ligand" and "structures_without_ligand" based on the following '
                 'list of identified ligands per structure:')
        df = pd.read_csv(str(output_directory) + '/pdbs_have_ligands.csv')
        st.table(df)


st.markdown('''##### Liganded Structures - Input Directory''')
st.write('Please select your input directory (where your liganded structure files are located) or just continue with the default location:')
st.write(f'(If you have used the tool above, this should be: "{str(output_directory)}/structures_with_ligand")')
# input_directory_liganded = st_directory_picker_input_liganded(initial_path=output_directory, key='input_directory_liganded')
# st.session_state.input_directory_liganded = str(input_directory_liganded)
# If output directory is selected, set the default liganded directory
if st.session_state.output_directory:
    default_liganded_dir = os.path.join(st.session_state.output_directory, "structures_with_ligand")
    if os.path.isdir(default_liganded_dir):
        st.session_state.input_directory_liganded = default_liganded_dir
    else:    # If "structures_with_ligand" directory does not exist in output_directory, set output_directory as default
        st.session_state.input_directory_liganded = st.session_state.output_directory

    # if st.session_state.input_directory_liganded is None:
    #     st.session_state.input_directory_liganded = default_liganded_dir

    st.write(f"Default Liganded Structures Directory: `{st.session_state.input_directory_liganded}`")

    # Custom Liganded Directory Selection
    if not is_running_in_docker():
        if st.button("Select Custom Directory for Liganded Structures"):
            select_folder_non_docker("input_directory_liganded")
    else:
        select_folder_in_docker("input_directory_liganded", "Select Custom Liganded Directory")

    st.write(f"📂 Selected Custom Directory of Liganded Structures: `{st.session_state.input_directory_liganded}`")

else:
    st.warning("⚠️ Please select an output directory first.")


if st.session_state.input_directory_liganded:
    input_directory_liganded = st.session_state.input_directory_liganded


st.markdown('''##### Variables''')
cutoff = st.number_input("Distance cutoff for binding site detection (in angstrom):", value=3.5, step=0.1, format="%.1f")
cutoff = round(cutoff,1) # round to one decimal, because when using the +/- buttons it goes out of 0.1 steps
st.write('Distance cutoff for binding site detection is ', cutoff, ' A.')

outputdir_BindingSite_ident = str(output_directory) + '/BindingSite_ident_Bio3D'


def run_bs_ident_Bio3D():
    result = subprocess.run(
        ['Rscript', str(parentfilepath)+'/identify_binding_site_bio3d.R', '-i', input_directory_liganded, '-o', outputdir_BindingSite_ident, '-d', str(cutoff)])
    st.write("Structures are taken from ", input_directory_liganded)
    st.write("Bio3D calculations are running...")
    st.write("Output files are saved in: ", outputdir_BindingSite_ident)
    st.session_state.BSidentifyisdone = True
    try:
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath)+'/tools/combine_pdb_ligands.py', '-i', input_directory_liganded, '-o', outputdir_BindingSite_ident, '-p', outputdir_BindingSite_ident+'/binding_site_interface_labelled_frequency.pdb'])
    finally:
        pass

st.markdown("##### Calculation and outputs:\n"
            " - Binding site residues are identified based on distance cutoff from ligand atoms.\n"
            " - Radius of gyration is calculated for the identified binding site residues per structure.\n"
            " - Output table 'binding_site_residues.tsv' contains the following information per structure: "
            "pdbID, ligIDs, number of residues, C-alpha radius of gyration, residue names of identified binding site.\n"
            " - The frequency of identification within the binding site is calculated for each residue and saved as "
            "data table 'binding_site_residue_occurrence_frequency.csv' and visualised in histograms.\n"
            " - Identified residues are labelled for each structure using the b-factor column of the given structure "
            " and pdb files are saved in subdirectory 'structures_labeled_binding_site'.\n"
            " - Three further pdb files are saved containing the coordinates of the reference structure and in the "
            "b-factor column either whether a residue was identified as being part of the binding site in any "
            "structure (values of 0 or 1) 'binding_site_interface_labelled_occurrence.pdb', or the calculated "
            "occurrence frequencies 'binding_site_interface_labelled_frequency.pdb', or as percentage "
            "'binding_site_interface_labelled_percentage.pdb'.\n"
            " - Ligand-Residue interactions are plotted as heatmap with pattern clustering to allow detection of "
            "ligand dependant characteristics, such as eventual sub-pockets 'interaction_heatmap.png'.\n"
            " - The output list 'binding_site_residue_numbers.txt' contains only the binding site residue numbers "
            "and will be used for binding site superpositioning and flexibility analysis (see next steps).")

st.button('Run Binding Site Identification', key='ident_bs_btn', on_click=run_bs_ident_Bio3D)

if st.session_state.BSidentifyisdone == True:
    initial_residue_numbers_file = outputdir_BindingSite_ident + '/binding_site_residue_numbers.txt'
    changed_residue_numbers_file = outputdir_BindingSite_ident + '/binding_site_residue_numbers_edited.txt'
    if (st.session_state.loadinitialresidues is False) and (os.path.exists(changed_residue_numbers_file)):
        st.session_state.residue_numbers_file = changed_residue_numbers_file
    else:
        st.session_state.residue_numbers_file = initial_residue_numbers_file
    with st.container(border=True, height=600):
        st.subheader("Binding Site Identification Results")
        st.write("More data is provided in structural and table format in the output directory.")
        st.write("How many residues are involved in ligand binding per structure?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_count.png', caption='binding residues counts')
        st.write("Which residues are are most frequently involved in ligand binding?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_frequency.png', caption='binding residues frequency')
        st.write("Reference structure residues colored by frequency (percentage) of involvement in binding site")
        b_fac_on_structure_vis(outputdir_BindingSite_ident + '/binding_site_interface_labelled_percentage.pdb')
        st.write("Which ligands share the same residue binding pattern?")
        try:
            st.image(outputdir_BindingSite_ident + '/interaction_heatmap.png', caption='ligand-residue interactions')
        except:
            st.write("ERROR: No output available.")
        st.write("How extended are binding sites (approximated through radius of gyration of binding residues C-alphas)?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_ca_rgyr.png', caption='Ca radius of gyration')

st.divider()
#######################################################################################################################
st.markdown("##### Optional: Edit the identified binding site residues\n"
            "Identified binding site residues are displayed below and ordered from most frequent to least frequent. "
            "You can edit them if you want (Format: One residue number per line). "
            "Make sure to only include residues that are present in all your structures.")

def load_file(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return file.read()
    return ""

def save_file(file_path, content):
    with open(file_path, 'w') as file:
        file.write(content)


if st.button("Reload initially identified residues", key="load_initial_residues_btn"):
    st.session_state.loadinitialresidues = True
    st.session_state.residue_numbers_file = initial_residue_numbers_file

# Load the file contents
file_contents = load_file(st.session_state.residue_numbers_file)
# Display the file contents in an editable text area
edited_contents = st.text_area("Edit residues to be used:", value=file_contents, height=200)

# Add a button to save the changes
if st.button("Save Changes", key="save_session_state_btn"):
    st.session_state.loadinitialresidues = False
    # changed_residue_numbers_file = outputdir_BindingSite_ident+'/binding_site_residue_numbers_edited.txt'
    save_file(changed_residue_numbers_file, edited_contents)
    st.success("File saved at "+changed_residue_numbers_file)
    st.session_state.residue_numbers_file = changed_residue_numbers_file

st.divider()
#######################################################################################################################

toc.subheader("2 - Superimpose using only Binding Site Residues")
st.markdown('''
_Superimpose all your structures (apo or liganded) locally using only Binding Site residues with one of the provided methods OR use the globally superimposed structures by selecting "None"._
''')
superimp_method_bs = st.radio(
    "Select the method you would like to use",
    ["None", "Bio3D", "ProDy"],
    captions = ["Skip and use globally superimposed structures.",
                "Superimpose using only binding site residues with Bio3D.",
                "Superimpose using only binding site residues with ProDy. (not recommended, slow)"],
    key="superimp_method_bs")
st.write("Your choice is: ", superimp_method_bs)

def superimpose_bs(superimp_method_bs):
    if superimp_method_bs == "Bio3D":
        result = subprocess.run(
            ["Rscript", str(parentfilepath)+"/superimpose_binding_site_bio3d.R", '-i', str(input_directory),
             '-o', str(output_directory), '-b', st.session_state.residue_numbers_file])
        st.write("Bio3D calculations are running...")
        st.session_state.superimposed_bs = str(output_directory)+'/superimposed_on_bs'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed_bs)
    if superimp_method_bs == "ProDy":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath)+"/superimpose_binding_site_prody.py", '-i', str(input_directory),
             '-o', str(output_directory), '-b', st.session_state.residue_numbers_file])
        st.write("ProDy calculations are running...")
        st.session_state.superimposed_bs = str(output_directory)+'/superimposed_on_bs/split_ensemble'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed_bs)
    if superimp_method_bs == "None":
        st.session_state.superimposed_bs = str(input_directory_liganded)
        st.write("Superimposed structures from input directory/liganded will be used.")


# default setting is "None" (no superpositioning required)
if st.session_state.superimposed_bs == "":
    st.session_state.superimposed_bs = str(input_directory_liganded)
    st.write("Superimposed structures from input directory/liganded will be used.")

st.button('Confirm/Go!', key="superimp_bs_btn", on_click=superimpose_bs, args=[superimp_method_bs])

superimposed_bs = st.session_state.superimposed_bs

st.divider()
#######################################################################################################################

toc.subheader("3 - Flexibility Analysis of Binding Site Residues")
st.markdown('''
_Investigate structural flexibility/variability locally._
''')

#toc.subheader("A) Focused variability analysis")
st.markdown('''#### A) Focused variability analysis''')
st.markdown("[Used package/tool: Bio3D (R)]")
st.markdown(" - Root Mean Square Fluctuations (RMSF)\n"
            " - Principal Component Analysis (PCA) on coordinates (all-atom)")

outputdir_BindingSite_analysis_Bio3D = str(output_directory) + '/BindingSite_analysis_Bio3D'

def run_analysis_BindingSite_Bio3D():
    print(superimposed_bs)
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/analyse_flex_binding_site_bio3d.R",
         '-i', superimposed_bs, '-o', outputdir_BindingSite_analysis_Bio3D,
         '-b', str(st.session_state.residue_numbers_file)])
    st.write("Superimposed structures are taken from ", superimposed_bs)
    st.write("Bio3D calculations are running...")
    st.write("Output files are saved in: ", outputdir_BindingSite_analysis_Bio3D)
    st.session_state.BSanalysisdone = True
    #show_report(parentpath=outputdir_BindingSite_analysis, filename="/bindingsite_analysis_bio3d_html_report.html")


st.button('Run', key="analyse_BindingSite_Bio3D_btn", on_click=run_analysis_BindingSite_Bio3D)

if st.session_state.BSanalysisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Binding Site Analysis Results")
        st.write("More data is provided in structural format in the output directory.")
        tab1, tab2 = st.tabs(["All-atom RMSF & PCA", "Backbone RMSF & PCA"])
        with tab1:
            st.markdown("#### - RMSF analysis (all-atom)")
            try:
                st.image(outputdir_BindingSite_analysis_Bio3D + '/RMSF_bsite_ascending.png', caption='binding residues RMSFs')
            except:
                st.write("ERROR: No output available.")
            st.markdown("#### - PCA on coordinates (all-atom)")
            try:
                st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_bsite_allatom.png', caption='binding residues PCA')
                st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_atom_contribution_bsite_allatom.png', caption='PCA loadings')
                st.write("Sampling along PC1")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC1_bsite_allatom.pdb')
                st.write("Sampling along PC2")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC2_bsite_allatom.pdb')
                st.write("Sampling along PC3")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC3_bsite_allatom.pdb')
            except:
                st.write("ERROR: No output available.")
        with tab2:
            st.markdown("#### - RMSF analysis (backbone)")
            try:
                st.image(outputdir_BindingSite_analysis_Bio3D + '/RMSF_bsite_backbone_ascending.png', caption='binding residues RMSFs')
            except:
                st.write("ERROR: No output available.")
            st.markdown("#### - PCA on coordinates (backbone)")
            try:
                st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_bsite_backbone.png', caption='binding residues PCA')
                st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_atom_contribution_bsite_backbone.png', caption='PCA loadings')
                st.write("Sampling along PC1")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC1_bsite_backbone.pdb')
                st.write("Sampling along PC2")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC2_bsite_backbone.pdb')
                st.write("Sampling along PC3")
                multimodel_animation(outputdir_BindingSite_analysis_Bio3D + '/PC3_bsite_backbone.pdb')
            except:
                st.write("ERROR: No output available.")

# toc.subheader("4 - SASA Analysis of Binding Site Residues")
# st.markdown('''
# _Investigate Solvent Accessible Surface Area (SASA) analysis of binding site residues._
# ''')


st.divider()
#######################################################################################################################

toc.header(":five: Conserved Water Analysis")
st.markdown('''
_Investigate conserved water molecules in your ensemble._
''')
# toc.subheader("4 - Conserved Water Analysis")
st.markdown("[Used package/tool: vanddraabe (R - modified) and Bio3D (R)]")
st.markdown(" - Identification of conserved water molecules by clustering\n ")
st.markdown("Reference for vanddraabe's water analysis can be found [here](https://github.com/exeResearch/vanddraabe), [here](https://www.rdocumentation.org/packages/vanddraabe/versions/1.1.1) "
            "and [here](http://www.exeresearch.com/uploads/6/8/2/0/6820495/thrombin10_conservedwater_vignette.pdf).\n  "
            "The script `Conserved_Waters.R` from the vanddraabe package has been modified to enable analysis of custom ensembles (without metadata from the RCSB-PDB).\n"
            "")

outputdir_water_analysis = str(output_directory) + '/Analysis_Waters'

if st.button('Run', key="run_water_analysis_btn"):
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/analyse_conserved_waters.R", '-i', str(superimposed), '-o',
         str(outputdir_water_analysis)])
    st.write("Superimposed structures are taken from ", superimposed)
    st.write("Output files are saved in: ", outputdir_water_analysis)
    st.session_state.wateranalysisdone = True

if st.session_state.wateranalysisdone == True:
    with st.container(border=True, height=600):
        try:
            st.subheader("Conserved Water Analysis Results")
            st.write("Superimposed structures were taken from ", superimposed)
            st.write("Output files can be found in: ", outputdir_water_analysis)
            st.write("Detailed water analysis results are provided as Excel file `WaterAnalysis_DATA_RESULTS.xlsx` in the output directory.\n"
                     "More data is provided in structural format in the output directory:\n"
                     "Two pdb structure files, `WaterAnalysis_ConservedWaters_ALL.pdb` and `WaterAnalysis_ConservedWaters_PASSED.pdb` "
                     "contain the water molecules with conservation saved as occupancy value. "
                     "A PyMol script highlighting the conserved waters at the ligand binding site of the reference structure "
                     "can be found at `WaterAnalysis_ConservedWaters_PASSED_PyMOL_white_background_.pml`.")
            st.image(outputdir_water_analysis + '/ConservationPlot.png', caption='ConservationPlot')
            st.image(outputdir_water_analysis + '/OccupancyBarplot.png', caption='OccupancyBarplot')
            st.image(outputdir_water_analysis + '/MobilityBarplot.png', caption='MobilityBarplot')
            st.image(outputdir_water_analysis + '/BvalueBarplot.png', caption='BvalueBarplot')
            st.image(outputdir_water_analysis + '/BvalueBarplot_calculated_values.png', caption='BvalueBarplot_calculated_values')
            st.image(outputdir_water_analysis + '/nBvalueBarplot.png', caption='nBvalueBarplot')
            st.image(outputdir_water_analysis + '/MobNormBvalEvalPlots.png', caption='MobNormBvalEvalPlots')
        except:
            st.write("ERROR: No output available.")



st.divider()
#######################################################################################################################

toc.header(":six: Compare with predicted flexibility")
st.markdown('''
_Compare differences in structural flexibility of your ensemble with other sources and simulation based results._
''')

toc.subheader("A) Normal Mode Analysis of elastic network models")
st.markdown("##### a) with Bio3D")
st.markdown("[Used package/tool: Bio3D (R)]")
st.markdown(" - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - Comparison of ensemble PCs and NMA modes of reference structure\n"
            " - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - [optional] Ensemble NMA (eNMA) based on ANM on all structures + clustering based on fluctuation "
            "similarity including RMSIP and Bhattacharyya coefficient comparison"
            )
eNMA = st.toggle('Activate eNMA calculation - only recommended for small-to-medium-size ensembles (<=50), as computationally expensive')
if eNMA:
    st.write('eNMA calculation activated!')

outputdir_NMA_Bio3D = str(output_directory) + '/Prediction_NMA_Bio3D'

if st.button('Run', key="run_NMA_bio3d_btn"):
    if eNMA:
        result = subprocess.run(
            ["Rscript", str(parentfilepath) + "/predict_flex_nma_bio3d.R", '-i', str(input_directory), '-o',
             str(outputdir_NMA_Bio3D), '-e'])
    else:
        result = subprocess.run(
            ["Rscript", str(parentfilepath) + "/predict_flex_nma_bio3d.R", '-i', str(input_directory), '-o',
             str(outputdir_NMA_Bio3D)])
    st.write("Output files are saved in: ", outputdir_NMA_Bio3D)
    st.session_state.NMABio3Disdone = True

if st.session_state.NMABio3Disdone == True:
    with st.container(border=True, height=600):
        try:
            st.subheader("Normal Mode Analysis (NMA) results")
            st.write("More data is provided in structural format in the output directory.")
            st.markdown("#### - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure")
            st.image(outputdir_NMA_Bio3D + '/ANM_NMA_reference_pdb.png', caption='ANM NMA')
            st.image(outputdir_NMA_Bio3D + '/ANM_NMA_dynamic_cross_correlations_reference_pdb.png',
                     caption='ANM NMA residue cross correlations')
            st.write("Reference structure colored by ANM fluctuations")
            b_fac_on_structure_vis(outputdir_NMA_Bio3D + '/ANM_fluctuations_onReference.pdb')
            st.write("First non-trivial ANM mode as interpolated trajectory")
            multimodel_animation(outputdir_NMA_Bio3D + '/ANM_NMA_reference_pdb_mode7_traj.pdb')
            st.markdown("#### - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure")
            st.image(outputdir_NMA_Bio3D + '/GNM_NMA_reference_pdb.png', caption='GNM NMA')
            st.image(outputdir_NMA_Bio3D + '/GNM_NMA_dynamic_cross_correlations_reference_pdb.png',
                     caption='GNM NMA residue cross correlations')
            st.write("Reference structure colored by GNM fluctuations")
            b_fac_on_structure_vis(outputdir_NMA_Bio3D + '/GNM_fluctuations_onReference.pdb')
        except:
            st.write("No NMA output available.")
        try:
            st.subheader("Ensemble Normal Mode Analysis (eNMA) results")
            st.image(outputdir_NMA_Bio3D + '/eNMA_fluctuations.png', caption='Ensemble Normal Mode fluctuations')
            st.write("Interpolated structures along eNMA modes as trajectory")
            multimodel_animation(outputdir_NMA_Bio3D + '/eNMA_traj.pdb')
        except:
            st.write("No eNMA output available.")



st.markdown("##### b) with ProDy")
st.markdown("[Used package/tool: ProDy (Python)]")
st.markdown(" - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - Dynamical Domain Decomposition of reference structure (using GNM)")

outputdir_NMA_ProDy = str(output_directory) + '/Prediction_NMA_ProDy'

if st.button('Run', key="run_NMA_prody_btn"):
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/predict_flex_nma_prody.py", '-i', str(input_directory), '-o',
         str(outputdir_NMA_ProDy)])
    st.write("Output files are saved in: ", outputdir_NMA_ProDy)
    st.session_state.NMAProDyisdone = True

if st.session_state.NMAProDyisdone == True:
    with st.container(border=True, height=600):
        try:
            st.subheader("Normal Mode Analysis (NMA) results")
            st.write("More data is provided in structural format in the output directory.")
            st.markdown("#### - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure")
            st.image(outputdir_NMA_ProDy + '/B-factors_vs_ANM-msfs.png', caption='ANM NMA RMSFs vs B-factors')
            st.write("Reference structure colored by ANM fluctuations")
            b_fac_on_structure_vis(outputdir_NMA_ProDy + '/ANM_msf_first20_on_bfactor.pdb')
            st.image(outputdir_NMA_ProDy + '/CrossCorrelations_ANM.png', caption='ANM NMA residue cross correlations')
            st.image(outputdir_NMA_ProDy + '/ContactMap_ANM.png', caption='ANM NMA residue contact map')
            st.markdown("#### - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure")
            st.image(outputdir_NMA_ProDy + '/B-factors_vs_GNM-msfs.png', caption='GNM NMA RMSFs vs B-factors')
            st.write("Reference structure colored by GNM fluctuations")
            b_fac_on_structure_vis(outputdir_NMA_ProDy + '/GNM_msf_first20_on_bfactor.pdb')
            st.image(outputdir_NMA_ProDy + '/CrossCorrelations_GNM.png', caption='GNM NMA residue cross correlations')
            st.image(outputdir_NMA_ProDy + '/ContactMap_GNM.png', caption='GNM NMA residue contact map')
        except:
            st.write("ERROR: No output available.")


toc.subheader("B) *All-atom* Normal Mode Analysis of elastic network models")
#st.markdown("##### a) with Bio3D")
st.markdown("[Used package/tool: Bio3D (R)]")
st.markdown(" - Anisotropic Network Model (ANM) all-atom Normal Mode Analysis (aaNMA) for reference structure\n"
            # " - Comparison of ensemble PCs and aaNMA modes of reference structure\n"
            " - [optional] all-atom Ensemble NMA (aaeNMA) based on ANM on all structures + clustering based on fluctuation "
            "similarity including RMSIP and Bhattacharyya coefficient comparison"
            )
aaeNMA = st.toggle('Activate aaeNMA calculation - only recommended for mini ensembles (<=5), as computationally very expensive')
if aaeNMA:
    st.write('All-atom eNMA calculation activated!')

outputdir_aaNMA_Bio3D = str(output_directory) + '/Prediction_aaNMA_Bio3D'

if st.button('Run', key="run_aaNMA_bio3d_btn"):
    if aaeNMA:
        result = subprocess.run(
            ["Rscript", str(parentfilepath) + "/predict_flex_aanma_bio3d.R", '-i', str(input_directory), '-o',
             str(outputdir_aaNMA_Bio3D), '-e'])
    else:
        result = subprocess.run(
            ["Rscript", str(parentfilepath) + "/predict_flex_aanma_bio3d.R", '-i', str(input_directory), '-o',
             str(outputdir_aaNMA_Bio3D)])
    st.write("Output files are saved in: ", outputdir_aaNMA_Bio3D)
    st.session_state.aaNMABio3Disdone = True

if st.session_state.aaNMABio3Disdone == True:
    with st.container(border=True, height=600):
        try:
            st.subheader("All-atom Normal Mode Analysis (aaNMA) results")
            st.write("More data is provided in structural format in the output directory.")
            st.markdown("#### - Anisotropic Network Model (ANM) All-atom Normal Mode Analysis (aaNMA) for reference structure")
            st.image(outputdir_aaNMA_Bio3D + '/ANM_aaNMA_reference_pdb.png', caption='ANM NMA')
            st.image(outputdir_aaNMA_Bio3D + '/ANM_aaNMA_dynamic_cross_correlations_reference_pdb.png',
                     caption='ANM aaNMA residue cross correlations')
            st.write("Reference structure colored by ANM all-atom fluctuations")
            b_fac_on_structure_vis(outputdir_aaNMA_Bio3D + '/ANM_aa_fluctuations_data_on_structure.pdb')
            st.write("First non-trivial ANM mode as interpolated trajectory")
            multimodel_animation(outputdir_aaNMA_Bio3D + '/ANM_aaNMA_reference_pdb_mode7_traj.pdb')
            # st.markdown("#### - Gaussian Network Model (GNM) All-atom Normal Mode Analysis (aaNMA) for reference structure")
            # st.image(outputdir_aaNMA_Bio3D + '/GNM_aaNMA_reference_pdb.png', caption='GNM NMA')
            # st.image(outputdir_aaNMA_Bio3D + '/GNM_aaNMA_dynamic_cross_correlations_reference_pdb.png',
            #          caption='GNM aaNMA residue cross correlations')
            # st.write("Reference structure colored by GNM all-atom fluctuations")
            # b_fac_on_structure_vis(outputdir_aaNMA_Bio3D + '/GNM_aa_fluctuations_onReference.pdb')
        except:
            st.write("No all-atom NMA output available.")
        try:
            st.subheader("All-atom Ensemble Normal Mode Analysis (eNMA) results")
            st.image(outputdir_NMA_Bio3D + '/aaeNMA_fluctuations.png', caption='All-atom Ensemble Normal Mode fluctuations')
            st.write("Interpolated structures along all-atom eNMA modes as trajectory")
            multimodel_animation(outputdir_NMA_Bio3D + '/aaeNMA_traj.pdb')
        except:
            st.write("No all-atom eNMA output available.")


toc.subheader("C) Essential Site Scanning Analysis")
st.markdown("[Used package/tool: ProDy (Python)]")
st.markdown('''Essential Site Scanning Analysis (ESSA) - an elastic network model (ENM)-based method - is performed 
            for all residues of the reference protein structure. 
            "Essential sites are here defined as residues that would significantly alter the protein’s global dynamics if bound to a ligand." 
            ([ESSA method](https://www.sciencedirect.com/science/article/pii/S2001037020303093 "ESSA publication")). 
            The previously identified binding site residues (as provided through the binding site identification step) 
            are also annotated on the output plot.
            ''')

outputdir_ESSA_ProDy = str(output_directory) + '/Prediction_ESSA_ProDy'

if st.button('Run', key="run_ESSA_prody_btn"):
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/predict_flex_binding_site_essa_prody.py",
         '-i', str(input_directory),
         '-o', str(outputdir_ESSA_ProDy),
         '-b', st.session_state.residue_numbers_file])
    st.write("Output files are saved in: ", outputdir_ESSA_ProDy)
    st.session_state.ESSAisdone = True

if st.session_state.ESSAisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Essential Site Scanning Analysis (ESSA) results")
        try:
            st.write("More data is provided in structural format in the output directory.")
            st.image(outputdir_ESSA_ProDy + '/ESSA_profile_of_reference_structure.png',
                     caption='ESSA profile with annotated binding residues')
            st.write("Reference structure colored by ESSA z-scores")
            b_fac_on_structure_vis(outputdir_ESSA_ProDy + '/reference_structure_gnm_zs.pdb')
        except:
            st.write("ERROR: No output available.")


# # Future projects
#
# toc.subheader("B) Monte-Carlo Sampling")
# if st.button('Run', key="run_monte-carlo_btn"):
#     st.write("Output files are saved in: ", output_directory)
#
#
# toc.subheader("C) Compare with PDBFlex cluster")
# st.markdown(" - Download PDBFlex cluster\n"
#             " - Root Mean Square Fluctuations (RMSF)\n"
#             " - Principal Component Analysis (PCA)")
#
# if st.button('Run', key="run_PDBFlex_btn"):
#     result = subprocess.run(
#         [f"{sys.executable}", str(parentfilepath) + "/compare_PDBFlex.py", '-i', str(input_directory), '-o',
#          str(output_directory)])
#     st.write("Output files are saved in: ", output_directory)
#
#
# toc.subheader("D) Compare with AlphaFold prediction")
# st.markdown(" - Download AlphaFold prediction\n"
#             " - Root Mean Square Fluctuations (RMSF)\n"
#             " - Principal Component Analysis (PCA)")
#
# if st.button('Run', key="run_AlphaFold_btn"):
#     result = subprocess.run(
#         [f"{sys.executable}", str(parentfilepath) + "/compare_AlphaFold.py", '-i', str(input_directory), '-o',
#          str(output_directory)])
#     st.write("Output files are saved in: ", output_directory)


st.divider()
#######################################################################################################################
st.write("***Funding***")
st.write("*This work was supported by the EMBL ARISE Fellowship from the European Union’s Horizon 2020 research "
         "and innovation program under the Marie Skłodowska-Curie grant agreement No 945405.*")


toc.toc()
