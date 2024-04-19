
import sys
import os
import subprocess

import pandas as pd
import streamlit as st
from stoc import stoc
#from directorypicker import st_directory_picker
from directorypicker import st_directory_picker_input, st_directory_picker_output, st_directory_picker_input_liganded
from pathlib import Path
from PIL import Image
#from pdf2image import convert_from_path
#import pypdfium2 as pdfium

## Personalisation

#icon = Image.open("icon.png")
st.set_page_config(
    page_title="EnsembleFlex", # the title of the page
    #page_icon=icon, # the page’s icon
    layout="centered", # ["centered", "wide"] # set Streamlit app to use wide mode
    initial_sidebar_state="auto", # whether sidebar will be initially loaded
)

# make app accessible only through localhost - Not working!
# st.set_option('browser.serverAddress', 'localhost')

filepath = Path(__file__).parent.resolve()
parentfilepath = Path(__file__).parent.parent.resolve()

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

# Session states initialisation
if 'superimposed' not in st.session_state:
    st.session_state.superimposed = ""
if 'Bio3Danalysisdone' not in st.session_state:
    st.session_state.Bio3Danalysisdone = False
if 'SASAanalysisdone' not in st.session_state:
    st.session_state.SASAanalysisdone = False
if 'PDBhasLigandSortdone' not in st.session_state:
    st.session_state.PDBhasLigandSortdone = False
if 'BSidentifydone' not in st.session_state:
    st.session_state.BSidentifydone = False
if 'superimposed_bs' not in st.session_state:
    st.session_state.superimposed_bs = ""
if 'BSanalysisdone' not in st.session_state:
    st.session_state.BSanalysisdone = False
if 'NMABio3Disdone' not in st.session_state:
    st.session_state.NMABio3Disdone = False
if 'NMAProDyisdone' not in st.session_state:
    st.session_state.NMAProDyisdone = False
if 'ESSAisdone' not in st.session_state:
    st.session_state.ESSAisdone = False

# def main():

#st.title("EnsembleFlex - Flexibility Analysis of Structure Ensembles")
toc.title("EnsembleFlex - Flexibility Analysis of Structure Ensembles")
st.markdown('''
_An interactive tool for analysing ensemble structures._
''')

toc.header("Introduction")
st.markdown('''
Welcome to an exploration of the dynamic world within proteins! 

Proteins, the fundamental molecular machines orchestrating the intricate dance of life, 
exhibit a remarkable diversity in structure and function. While classical structural biology 
has provided invaluable insights into the static three-dimensional architectures of proteins, 
the dynamic nature of these macromolecules is increasingly recognized as a critical determinant 
of their biological activity. Protein structural flexibility, the ability of proteins to undergo 
conformational changes while maintaining their functional integrity, lies at the heart of this 
dynamic tapestry. 

The inherent flexibility of proteins manifests at various levels, spanning from local fluctuations 
of individual amino acid side chains to large-scale domain movements. Understanding the dynamic 
behavior of proteins is essential for deciphering the mechanisms underlying their diverse functions, 
including enzymatic catalysis, molecular recognition, and signal transduction. Moreover, protein 
flexibility plays a pivotal role in modulating interactions with other biomolecules, enabling 
proteins to adapt to different cellular environments and respond to external stimuli.

Understanding protein structural flexibility is not just about unraveling molecular intricacies; 
it holds the key to deciphering diseases, designing targeted therapies, and engineering proteins 
with tailored functionalities.

This webpage serves as a gateway into the fascinating field of protein structural flexibility in an 
automated manner, 
where we unravel the dynamic tapestry of these molecular machines.
''')

st.markdown('### Output')
st.markdown('*The generated output directory will contain the following subdirectories:*\n  '
            'Uppercase directory names contain analysis output files of diverse format and '
            'lowercase directory names only contain .pdb structure files.')
st.code('''
EnsemblFlex  
├── superimposed  
├── Analysis_Bio3D  
│   └── pymol_pdbs  
├── Analysis_SASA_Biopython  
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

st.markdown('''#### Input Directory''')
st.write('Please select your input directory (where your structure files are located):')
input_directory = st_directory_picker_input(key="input_directory")

st.markdown('''#### Output Directory''')
st.write('Please select your output directory (where all calculation results and files will be saved):')
output_directory = st_directory_picker_output(key="output_directory")

st.divider()
#######################################################################################################################

toc.header(":two: Superimpose (global)")
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
        #sys.which("Rscript")
        # result = subprocess.Popen(["Rscript", str(parentfilepath)+"/superimpose_bio3d.R",
        #                            '-i', str(input_directory), '-o', str(output_directory)],
        #                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # result1 = result.communicate()
        result = subprocess.run(
            ["Rscript", str(parentfilepath)+"/superimpose_bio3d.R", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("Bio3D calculations are running...")
        # # Show stdout for external command
        # for line in iter(lambda: result.stdout.readline(), b""):
        #     st.text(line.decode("utf-8"))
        st.session_state.superimposed = str(output_directory)+'/superimposed'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed)
    if superimp_method == "ProDy":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath)+"/superimpose_prody.py", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("ProDy calculations are running...")
        st.session_state.superimposed = str(output_directory)+'/superimposed/split_ensemble'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed)
    if superimp_method == "None":
        st.session_state.superimposed = str(input_directory)
        st.write("Superimposed structures from input directory will be used.")

st.button('Confirm/Go!', key="superimp_btn", on_click=superimpose, args=[superimp_method])

superimposed = st.session_state.superimposed

st.divider()
#######################################################################################################################

toc.header(":three: Flexibility Analysis (global)")
st.markdown('''
_Investigate structural flexibility globally using the provided methods._
''')

st.markdown('''##### Variables''')
number_of_groups = st.number_input("Desired number of clusters\n\n"
                                   "(You may want to run it first with the default and adjust the value based on "
                                   "the output in subsequent runs.)", value=3)


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
        ["Rscript", str(parentfilepath) + "/analyse_flex_bio3d.R", '-i', superimposed, '-o', outputdirBio3D])
    st.write("Superimposed structures are taken from ", superimposed)
    st.write("Bio3D calculations are running...")
    st.write("Files are saved in: ", outputdirBio3D)
    st.session_state.Bio3Danalysisdone = True

def show_report(parentpath, filename):
    # Show generated html report
    path_to_html = str(parentpath) + filename
    if os.path.isfile(path_to_html):
        # Read file and keep in variable
        with open(path_to_html, 'r') as f:
            html_data = f.read()
        ## Show in app
        #st.header("Show an external HTML")
        st.components.v1.html(html_data, height=200)
    else:
        st.write("HTML report could not be generated.")

st.button('Run', key="Analyse_Bio3D_btn", on_click=run_analysis_Bio3D)
st.markdown("Results in form of images will be displayed below.  "
            "Please note that more information is available in form of data connected to structures (.pdb .pml) "
            "in the output folder.")

if st.session_state.Bio3Danalysisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Flexibility Analysis Results")
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Base info", "RMSF & RMSD analysis", "PCA", "UMAP analysis",
                                                      "All-atom PCA & DDM analysis", "Clustering"])
        with tab1:
            st.markdown("#### Base info: Alignment overview and B-factors of reference")
            st.image(outputdirBio3D + '/alignment_overview.png', caption='alignment overview')
            st.image(outputdirBio3D + '/B-factors.png', caption='B-factors of reference structure')
            #st.image("https://static.streamlit.io/examples/cat.jpg", width=200)
        with tab2:
            st.markdown("#### RMSF and RMSD analysis")
            st.markdown("#### - on backbone coordinates")
            st.image(outputdirBio3D + '/RMSF.png', caption='RMSF')
            st.image(outputdirBio3D + '/RMSD_hist.png', caption='RMSD histogram')
            st.image(outputdirBio3D + '/RMSD_heatmap.png', caption='RMSD heatmap')
            st.markdown("#### - on all-atom coordinates (including side chains)")
            st.image(outputdirBio3D + '/RMSF_allatom.png', caption='All-atom RMSF')
            st.image(outputdirBio3D + '/RMSD_hist_allatom.png', caption='All-atom RMSD histogram')
            st.image(outputdirBio3D + '/RMSD_heatmap_allatom.png', caption='All-atom RMSD heatmap')
        with tab3:
            st.markdown("#### Principal Component Analysis (PCA) on coordinates and on torsion data")
            st.markdown("#### - on coordinates (backbone)")
            st.image(outputdirBio3D + '/PCA.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_residue_contribution.png', caption='PCA residue contribution')
            st.image(outputdirBio3D + '/PCA_dendrogram.png', caption='PCA dendrogram')
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
            st.markdown("#### All-atom PCA")
            st.markdown("#### - on coordinates (all-atom)")
            st.image(outputdirBio3D + '/PCA_allatom.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_loadings_allatom.png', caption='PCA loadings')
            st.image(outputdirBio3D + '/PCA_dendrogram_allatom.png', caption='PCA dendrogram')
            st.markdown("#### - on difference distance matrices (all-atom)")
            st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix.png', caption='PCA overall info')
            st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix_loadings.png', caption='PCA loadings')
            st.image(outputdirBio3D + '/PCA_on_allatom_DifferenceDistanceMatrix_dendrogram.png', caption='PCA dendrogram')
            # st.markdown("#### All-atom Uniform Manifold Approximation and Projection (UMAP) analysis")
            # st.markdown("#### - on coordinates (all-atom)")
            # st.image(outputdirBio3D + '/UMAP_allatom.png', caption='2D UMAP plot')
            # st.image(outputdirBio3D + '/UMAP_dendrogram_allatom.png', caption='2D UMAP dendrogram')
            st.markdown("#### Difference Distance Matrix (DDM) analysis between structures (only all-atom)")
            try:
                st.image(outputdirBio3D + '/eDDM_complete.png', caption='Ensemble Difference Distance Matrix')
            except: pass
        with tab6:
            st.markdown("#### Overall clustering results")
            st.image(outputdirBio3D + '/cluster_attributions_heatmap.png', caption='Cluster attributions')


    show_report(parentpath=outputdirBio3D, filename="/analysis_bio3d_html_report.html")


st.markdown("##### Helper Tool: Sort structures based on clustering")
st.markdown("[Used package/tool:  (Python)]")

cluster_option = st.selectbox(
    'Which clustering result would you like to use to sort structures into subdirectories?',
    ("Consensus_Cluster", "backbone_RMSD", "backbone_PCA_onCoords", "backbone_UMAP_onCoords", "backbone_PCA_onTorsion",
     "allatom_RMSD", "allatom_PCA_onCoords", "allatom_UMAP_onCoords", "allatom_PCA_onDist"))
st.write('Your choice is:', cluster_option)

if st.button('Run', key="Subset_clusters"):
    outputdirBio3D_clusters = str(output_directory) + '/Analysis_Bio3D/' + cluster_option
    cluster_dataframe = str(output_directory) + '/Analysis_Bio3D/cluster_attributions_with_consensus.csv'
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
# st.button('Run', key="Analyse_ProDy_btn", on_click=run_analysis_ProDy)


toc.subheader("B) SASA Difference Analysis")
st.markdown("[Used package/tool: Biopython (Python)]")
st.markdown('''
_Investigate differences in Solvent Accessible Surface Area (SASA) of protein residues across the whole ensemble._  
_This highlights residues that experience differences in solvent accessibility across the ensemble._  
''')

outputdirSASA = str(output_directory) + '/Analysis_SASA_Biopython'

def run_analysis_SASA():
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/analyse_flex_sasa_biopython.py", '-i', superimposed, '-o', outputdirSASA])
    st.write("Superimposed structures are taken from ", superimposed)
    st.write("SASA calculations are done.")
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/tools/data_on_structure.R", '-i', superimposed, '-o', outputdirSASA,
         '-d', outputdirSASA+"/SASA_global.csv", '-c', "sd"])
    st.write("SASA standard deviation per residue saved in b-factor on reference structure.")
    st.write("Files are saved in: ", outputdirSASA)
    st.session_state.SASAanalysisdone = True

st.button('Run', key="Analyse_SASA_btn", on_click=run_analysis_SASA)

if st.session_state.SASAanalysisdone == True:
    with st.container(border=True, height=600):
        st.markdown("#### SASA difference analysis result")
        st.write("More data is provided in structural and table format in the output directory.")
        st.image(outputdirSASA + '/SASA_sd.png', caption='SASA standard deviation per residue')

st.divider()
#######################################################################################################################


toc.header(":four: Binding Site Analysis")
#st.header(":four: Binding Site Analysis")

toc.subheader("1 - Identify the Binding Site")
st.markdown("[Used package/tool: Bio3D (R)]")
st.markdown('''
_Please provide a directory of structures with ligands (no apo-structures)._  
_Be aware that any molecule that is not protein, nucleic acid or water will be considered as "ligand"._
''')


st.markdown('''##### Helper Tool: Sort liganded structures''')
st.markdown("[Used package/tool: Bio3D (R)]")
st.write('This tool automatically sorts copies of your superimposed structures into the subfolders "structures_with_ligand" and "structures_without_ligand".')

def sortPDBs():
    result = subprocess.run(
        ['Rscript', str(parentfilepath) + '/tools/sort_pdbs_has_ligand.R', '-i', superimposed, '-o', output_directory])
    liganded = str(output_directory) + '/structures_with_ligand'
    not_liganded = str(output_directory) + '/structures_without_ligand'
    st.write("The subfolders are located in ", output_directory)
    st.session_state.PDBhasLigandSortdone = True

st.button('Sort structures', key='sort_btn', on_click=sortPDBs)
if st.session_state.PDBhasLigandSortdone == True:
    with st.container(border=True, height=300):
        st.write('Structures have been sorted into respective subfolders '
                 '"structures_with_ligand" and "structures_without_ligand" based on the following '
                 'list of identified ligands per structure:')
        df = pd.read_csv(str(output_directory) + '/pdbs_have_ligands.csv')
        st.table(df)


st.markdown('''##### Liganded Structures - Input Directory''')
st.write('Please select your input directory (where your liganded structure files are located):')
input_directory_liganded = st_directory_picker_input_liganded(key='input_directory_liganded')

st.markdown('''##### Variables''')
cutoff = st.number_input("Distance cutoff for binding site detection (in angstrom):", value=4.0, step=0.1, format="%.1f")
cutoff = round(cutoff,1) # round to one decimal, because when using the +/- buttons it goes out of 0.1 steps
st.write('Distance cutoff for binding site detection is ', cutoff, ' A.')

outputdir_BindingSite_ident = str(output_directory) + '/BindingSite_ident_Bio3D'


def run_bs_ident_Bio3D():
    result = subprocess.run(
        ['Rscript', str(parentfilepath) + '/identify_binding_site_bio3d.R', '-i', input_directory_liganded, '-o', outputdir_BindingSite_ident, '-d', str(cutoff)])
    st.write("Structures are taken from ", input_directory_liganded)
    st.write("Bio3D calculations are running...")
    st.write("Files are saved in: ", outputdir_BindingSite_ident)
    st.session_state.BSidentifydone = True
    # show binding site csv file
    #csvfile = pd.read_csv(outputdir_BindingSite_ident+"/binding_site_residues.csv")  # path folder of the data file
    #st.write(csvfile)
    # show binding frequency plot
    image = Image.open(outputdir_BindingSite_ident+'/Histogram_binding_residues_percentage_colored.png')
    st.image(image, caption='Histogram of binding residues')
    #show_report(parentpath=outputdir_BindingSite_ident, filename="/analysis_bio3d_html_report.html")

st.markdown("##### Calculation and outputs:\n"
            " - Binding site residues are identified based on distance cutoff from ligand atoms.\n"
            " - Radius of gyration is calculated for the identified binding site residues per structure.\n"
            " - Output table 'binding_site_residues.tsv' contains the following information per structure: "
            "pdbID, ligIDs, number of residues, C-alpha radius of gyration, residue names of identified binding site.\n"
            " - The frequency of identification within the binding site is calculated for each residue and saved as "
            "data table 'binding_site_residue_occurrence_frequency.csv' and visualised in histograms.\n"
            " - Identified residues are labeled for each structure using the b-factor column of the given structure "
            " and pdb files are saved in subdirectory 'structures_labeled_binding_site'.\n"
            " - Three further pdb files are saved containing the coordinates of the reference structure and in the "
            "b-factor column either whether a residue was identified as being part of the binding site in any "
            "structure (values of 0 or 1) 'binding_site_interface_labelled_occurrence.pdb', or the calculated "
            "occurrence frequencies 'binding_site_interface_labelled_frequency.pdb', or as percentage "
            "'binding_site_interface_labelled_percentage.pdb'.\n"
            " - Ligand-Residue interactions are plotted as heatmap with pattern clustering (to allow detection of "
            "eventual sub pockets) 'interaction_heatmap.png'.\n"
            " - The output list 'binding_site_residue_numbers.txt' contains only the binding site residue numbers "
            "and will be used for binding site superpositioning and flexibility analysis (see next steps).")

st.button('Run Binding Site Identification', key='BS_ident_btn', on_click=run_bs_ident_Bio3D)

if st.session_state.BSidentifydone == True:
    with st.container(border=True, height=600):
        st.subheader("Binding Site Identification Results")
        st.write("More data is provided in structural and table format in the output directory.")
        st.write("How many residues are involved in ligand binding per structure?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_count.png', caption='binding residues counts')
        st.write("Which residues are are most frequently involved in ligand binding?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_frequency.png', caption='binding residues frequency')
        st.write("Which ligands share the same residue binding pattern?")
        st.image(outputdir_BindingSite_ident + '/interaction_heatmap.png', caption='ligand-residue interactions')
        st.write("How extended are binding sites (approximated through radius of gyration of binding residues C-alphas)?")
        st.image(outputdir_BindingSite_ident + '/histogram_binding_residues_ca_rgyr.png', caption='Ca radius of gyration')

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
            ["Rscript", str(parentfilepath)+"/superimpose_binding_site_bio3d.R", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("Bio3D calculations are running...")
        st.session_state.superimposed_bs = str(output_directory)+'/superimposed_binding_site'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed_bs)
    if superimp_method_bs == "ProDy":
        result = subprocess.run(
            [f"{sys.executable}", str(parentfilepath)+"/superimpose_binding_site_prody.py", '-i', str(input_directory), '-o', str(output_directory)])
        st.write("ProDy calculations are running...")
        st.session_state.superimposed_bs = str(output_directory)+'/superimposed_binding_site/split_ensemble'
        st.write("Superimposed structures are saved in: ", st.session_state.superimposed_bs)
    if superimp_method_bs == "None":
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
    result = subprocess.run(
        ["Rscript", str(parentfilepath) + "/analyse_flex_binding_site_bio3d.R",
         '-i', superimposed_bs, '-o', outputdir_BindingSite_analysis_Bio3D,
         '-b', outputdir_BindingSite_ident+'/binding_site_residue_numbers.txt'])
    st.write("Superimposed structures are taken from ", superimposed_bs)
    st.write("Bio3D calculations are running...")
    st.write("Files are saved in: ", outputdir_BindingSite_analysis_Bio3D)
    st.session_state.BSanalysisdone = True
    #show_report(parentpath=outputdir_BindingSite_analysis, filename="/bindingsite_analysis_bio3d_html_report.html")


st.button('Run', key="Analyse_BindingSite_Bio3D_btn", on_click=run_analysis_BindingSite_Bio3D)

if st.session_state.BSanalysisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Binding Site Analysis Results")
        st.write("More data is provided in structural format in the output directory.")
        st.markdown("#### - RMSF analysis")
        st.image(outputdir_BindingSite_analysis_Bio3D + '/RMSF_bsite2.png', caption='binding residues RMSFs')
        st.markdown("#### - PCA on coordinates (all-atom)")
        st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_bsite_allatom.png', caption='binding residues PCA')
        st.image(outputdir_BindingSite_analysis_Bio3D + '/PCA_atom_contribution_bsite_allatom.png', caption='PCA loadings')

# toc.subheader("4 - SASA Analysis of Binding Site Residues")
# st.markdown('''
# _Investigate Solvent Accessible Surface Area (SASA) analysis of binding site residues._
# ''')


st.divider()
#######################################################################################################################

toc.header(":five: Compare with predicted flexibility")
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
eNMA = st.toggle('Activate eNMA calculation - only recommended for small-to-medium-size ensembles, as computationally expensive')
if eNMA:
    st.write('eNMA calculation activated!')

outputdir_NMA_Bio3D = str(output_directory) + '/Prediction_NMA_Bio3D'

if st.button('Run', key="NMA_bio3d_btn"):
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
        st.subheader("Normal Mode Analysis (NMA) results")
        st.write("More data is provided in structural format in the output directory.")
        st.markdown("#### - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure")
        st.image(outputdir_NMA_Bio3D + '/ANM_NMA_reference_pdb.png', caption='ANM NMA')
        st.image(outputdir_NMA_Bio3D + '/ANM_NMA_dynamic_cross_correlations_reference_pdb.png',
                 caption='ANM NMA residue cross correlations')
        st.markdown("#### - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure")
        st.image(outputdir_NMA_Bio3D + '/GNM_NMA_reference_pdb.png', caption='GNM NMA')
        st.image(outputdir_NMA_Bio3D + '/GNM_NMA_dynamic_cross_correlations_reference_pdb.png',
                 caption='GNM NMA residue cross correlations')



st.markdown("##### b) with ProDy")
st.markdown("[Used package/tool: ProDy (Python)]")
st.markdown(" - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure\n"
            " - Dynamical Domain Decomposition of reference structure (using GNM)")

outputdir_NMA_ProDy = str(output_directory) + '/Prediction_NMA_ProDy'

if st.button('Run', key="NMA_prody_btn"):
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/predict_flex_nma_prody.py", '-i', str(input_directory), '-o',
         str(outputdir_NMA_ProDy)])
    st.write("Output files are saved in: ", outputdir_NMA_ProDy)
    st.session_state.NMAProDyisdone = True

if st.session_state.NMAProDyisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Normal Mode Analysis (NMA) results")
        st.write("More data is provided in structural format in the output directory.")
        st.markdown("#### - Anisotropic Network Model (ANM) C-alpha Normal Mode Analysis (NMA) for reference structure")
        st.image(outputdir_NMA_ProDy + '/B-factors_vs_ANM-msfs.png', caption='ANM NMA RMSFs vs B-factors')
        st.image(outputdir_NMA_ProDy + '/CrossCorrelations_ANM.png', caption='ANM NMA residue cross correlations')
        st.image(outputdir_NMA_ProDy + '/ContactMap_ANM.png', caption='ANM NMA residue contact map')
        st.markdown("#### - Gaussian Network Model (GNM) C-alpha Normal Mode Analysis (NMA) for reference structure")
        st.image(outputdir_NMA_ProDy + '/B-factors_vs_GNM-msfs.png', caption='GNM NMA RMSFs vs B-factors')
        st.image(outputdir_NMA_ProDy + '/CrossCorrelations_GNM.png', caption='GNM NMA residue cross correlations')
        st.image(outputdir_NMA_ProDy + '/ContactMap_GNM.png', caption='GNM NMA residue contact map')


toc.subheader("B) Essential Site Scanning Analysis")
st.markdown("[Used package/tool: ProDy (Python)]")
st.markdown('''Essential Site Scanning Analysis (ESSA) - an elastic network model (ENM)-based method - is performed 
            for all residues of the reference protein structure. 
            "Essential sites are here defined as residues that would significantly alter the protein’s global dynamics if bound to a ligand." 
            ([ESSA method](https://www.sciencedirect.com/science/article/pii/S2001037020303093 "ESSA publication")). 
            The previously identified binding site residues (as provided through the binding site identification step) 
            are also annotated on the output plot.
            ''')

outputdir_ESSA_ProDy = str(output_directory) + '/Prediction_ESSA_ProDy'

if st.button('Run', key="ESSA_prody_btn"):
    result = subprocess.run(
        [f"{sys.executable}", str(parentfilepath) + "/predict_flex_binding_site_essa_prody.py",
         '-i', str(input_directory),
         '-o', str(outputdir_ESSA_ProDy),
         '-b', outputdir_BindingSite_ident+'/binding_site_residue_numbers.txt'])
    st.write("Output files are saved in: ", outputdir_ESSA_ProDy)
    st.session_state.ESSAisdone = True

if st.session_state.ESSAisdone == True:
    with st.container(border=True, height=600):
        st.subheader("Essential Site Scanning Analysis (ESSA) results")
        st.write("More data is provided in structural format in the output directory.")
        st.image(outputdir_ESSA_ProDy + '/ESSA_profile_of_reference_structure.png',
                 caption='ESSA profile with annotated binding residues')

# toc.subheader("B) Monte-Carlo Sampling")
# if st.button('Run', key="monte-carlo_btn"):
#     st.write("Output files are saved in: ", output_directory)
#
#
# toc.subheader("C) Compare with PDBFlex cluster")
# st.markdown(" - Download PDBFlex cluster\n"
#             " - Root Mean Square Fluctuations (RMSF)\n"
#             " - Principal Component Analysis (PCA)")
#
# if st.button('Run', key="PDBFlex_btn"):
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
# if st.button('Run', key="AlphaFold_btn"):
#     result = subprocess.run(
#         [f"{sys.executable}", str(parentfilepath) + "/compare_AlphaFold.py", '-i', str(input_directory), '-o',
#          str(output_directory)])
#     st.write("Output files are saved in: ", output_directory)


st.divider()
#######################################################################################################################


# toc.generate()
toc.toc()

# if __name__ == "__main__":
#     main()