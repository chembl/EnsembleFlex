#!/usr/bin/env python

"""
Generate HTML report file of Bio3D Flexibility Analysis.

Usage:
    python3 analysis_flex_bio3d_reporting.py <input_directory> <output_directory> <clustering_groups>

Example:
    python3 analysis_flex_bio3d_reporting.py EnsemblFlex/superimposed EnsemblFlex/Analysis_Bio3D 3
"""

import glob
import os
import sys


def generate_report(str_input_path, output_path, clustering_groups):

    os.chdir(output_path) # change working directory to output_path

    # 1. Set up multiple variables to store the titles, text within the report
    page_title_text = 'Bio3D Analysis Report'
    title_text = 'Bio3D Analysis Report'
    text_intro = '...'
    # list of pdb files
    pdbfiles = glob.glob(str_input_path+"/*.pdb")
    filenames = [os.path.basename(x) for x in pdbfiles]
    #pdbfiles = args.pdbfiles
    text_rmsd = '...'
    text_pca = '...'
    #clustering_groups = clustering_groups

    # 2. Combine them together using a long f-string
    html = f'''
        <html>
            <head>
                <title>{page_title_text}</title>
            </head>
            <body>
                <h1>{title_text}</h1>
                <p>The following {len(pdbfiles)} files were included in the analysis:</p>
                <div style="height:200px;width:700px;border:1px solid #ccc;font:16px/26px Georgia, Garamond, Serif;overflow:auto;">
                {filenames}
                </div>
                <p>Structure {filenames[0]} served as reference structure for superpositioning and for all calculations \
                based on atomic coordinates.</p>
                
                B-factor fluctuations of the C&#945; atoms of the ensemble are shown below. 
                This may serve as comparison with the further flexibility analysis results.
                <img src='B-factors.png' width="700">
                
                <p>Besides root mean square fluctuations (RMSFs) Principal component analysis (PCA) and normal mode analysis (NMA) have emerged as two invaluable tools \
                for studying conformational changes in proteins.</p>
                
                <h2> 1. Root Mean Square Deviation (RMSD) Analysis </h2>
                RMSD is a standard measure of structural distance between coordinate sets. \
                The distribution of pairwise RMSD values for the structural ensemble is shown in the histogram
                <img src='RMSD_hist.png' width="700">
                and the RMSD value differences between structures are visualised in the heatmap, \
                where structures are ordered based on the RMSD values:
                <img src='RMSD_heatmap.png' width="700">
                <p>Hierarchical clustering is performed on the RMSD distances with the number of \
                groups = {clustering_groups}, provided by the user, and the structures are colored based on their \
                cluster attribution.</p>
                <img src='RMSD_clust.png' width="700">
                
                
                <h2> 2. Root Mean Square Fluctuation (RMSF) Analysis </h2>
                <p>Root mean square fluctuation (RMSF) is a numerical measurement similar to Root mean square deviation \
                (RMSD), but instead of indicating positional differences between entire structures over time, RMSF is a \
                calculation of individual residue flexibility, or how much a particular residue moves (fluctuates) within \
                an ensemble of structures. <br>
                RMSFs shown in the plot are calculated on the C&#945; atoms per residue.</p>
                <img src='RMSF.png' width="700">
                <p>RMSF values are saved in B-factor column of RMSFonReference.pdb for visualisation. \
                Just open the structure in a molecular visualisation program and color by b-factor.</p>
                
                
                <h2> 3. Principal Component Analysis (PCA) </h2>
                <p>Principal component analysis (PCA) is a mathematical technique, used to find patterns in high-dimensional \
                datasets, such as protein structures. It enables us to find relationships/patterns, which would be invisible \
                from a pure visual examination. It can be applied to ensembles of protein structures to detect the global, \
                correlated motions of the system (the principal components). PCA is a method that identifies the components \
                which account for the greatest amount of variability in the ensemble.</p>
                <p>In order to get a first overview on on the clustering of structures and their variability \
                PCA is performed on C&#945; coordinates of the superimposed ensemble and all ensemble structures are \
                plotted in PC space below. The proportion of variance shows how significant \
                the PCs are in explaining structural differences. \
                Explained variance in PCA is a statistical measure of how much variation in a dataset can be \
                attributed to each of the principal components (eigenvectors) generated by the PCA. To compute the \
                percentage of variance (information) accounted for by each principal component, the eigenvalue of each \
                component is divided by the sum of eigenvalues.</p>
                <img src='PCA.png' width="700">
                
                <p>Hierarchical clustering is performed on the distances in the PC1-PC2 subspace with the number of \
                groups = {clustering_groups}, provided by the user, and the structures are colored based on their \
                cluster attribution.</p>
                <img src='PCA_clust.png' width="700">

                <p>The following plot shows how much each residue contributes to PC1, PC2, and PC3.</p>
                <img src='PCA_residue_contribution.png' width="700">
                
                
                <h2> 4. Ensemble Normal Mode Analysis (eNMA) </h2>
                <p>Normal Mode Analysis (NMA) is an alternative method to study dynamics of molecules, and does not require \
                an ensemble of structures, but is calculated on a single structure instead. It is based on the theory of \
                vibration and conformational fluctuation is given by a superposition of normal modes. <br>
                The normal modes of each protein structure of the ensemble was calculated and aligned providing \
                aligned atomic fluctuations and aligned eigenvectors, describing the intrinsic dynamics of the protein structure.</p>
                <img src='eNMA_fluctuations.png' width="700">
                <p>Interpolated structures along eNMA modes are saved in the file "eNMA.pdb".</p>
                
                <h3> 4.1 Normal Mode Analysis (NMA) on the reference structure </h3>
                <p>NMA fluctuations of the reference structure are shown for comparison.</p>
                <img src='NMA_fluctuations_reference_pdb.png' width="700">
                <p>Interpolated structures along the NMA mode 7 are saved in the file "NMA_reference_pdb_mode7_traj.pdb" \
                and the respective visualisation of the vectors is saved in the PyMol script "NMA_reference_pdb_mode7.pml".</p>
                
                
                <h2> 5. Ensemble Difference Distance Matrix (eDDM) Analysis </h2>
                <p>The eDDM analysis compares structural ensembles under distinct ligation, activation, etc. conditions. \
                (At least two groups of structures are required.
                The following procedure is applied: The distance matrices are calculated on all heavy-atom coordinates, \
                then, PCA of the distance matrices is performed, followed by a conventional hierarchical clustering \
                in the PC1-PC2 subspace to get the intrinsic grouping of the structures. \
                Finally, the difference mean distance between groups are calculated for each residue pair \
                and statistical significance is assessed using a two-sample Wilcoxon test. \
                Long-distance pairs in all structures are omitted. 
                Those groups are then used for a </p>

            </body>
        </html>
        '''
    # 3. Write the html string as an HTML file
    with open('analysis_bio3d_html_report.html', 'w') as f:
        f.write(html)

    # # 4. Convert the HTML file to PDF
    # from weasyprint import HTML, CSS
    # css = CSS(string='''
    #     @page {size: A4; margin: 1cm;}
    #     th, td {border: 1px solid black;}
    #     ''')
    # HTML('analysis_bio3d_html_report.html').write_pdf('analysis_bio3d_pdf_report.pdf', stylesheets=[css])


## Running the function using arguments

str_input_path = sys.argv[1]
output_path = sys.argv[2]
clustering_groups = sys.argv[3]

print(str_input_path, output_path, clustering_groups)

generate_report(str_input_path, output_path, clustering_groups)