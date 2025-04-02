#!/usr/bin/env Rscript

# """
# Ensemble flexibility analysis using mainly the R package Bio3D.
#
# Usage:
#     Rscript analyse_flex_bio3d.R -i <input_directory> -o <output_directory>
#
# Example:
#     Rscript analyse_flex_bio3d.R -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_Bio3D
# """

args = commandArgs(trailingOnly=TRUE)

projectdir = getwd()
#print(projectdir)

## Installation of dependencies for R-only usage (without Conda environment)
#install.packages("devtools")
#library(devtools)
#devtools::install_bitbucket("Grantlab/bio3d", subdir = "bio3d-core", ref="core")
#devtools::install_bitbucket("Grantlab/bio3d-eddm")
#install.packages("bio3d", dependencies=TRUE)
#install.packages("optparse", dependencies=TRUE)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")
#install.packages("pheatmap")
#install.packages("ggplot2")

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(bio3d)
library(msa)
library(pheatmap)
library(umap) # UMAP dimension reduction
library(cluster) # clustering
library(clValid) # cluster validation


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="Bio3D_Analysis",
              help="output directory [default=%default]", metavar="character"),
  make_option(c("-n", "--ngroups"), type="integer", default=3,
              help="number of groups for clustering [default=%default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check provided output directory <opt$outdir>

if (opt$outdir == "Bio3D_Analysis"){
  # output directory is supposed to be located in input directory
  outdir <- file.path(opt$indir, opt$outdir)
  dir.create(outdir, showWarnings = FALSE)
} else if (startsWith(opt$outdir, "./")){
  # if only relative output foldername is provided use current working directory as basepath
  sub_dir <- strsplit(opt$outdir, split='./', fixed=TRUE)
  outdir <- file.path(getwd(), sub_dir)
  dir.create(outdir, showWarnings = FALSE)
} else if (length(unlist(strsplit(opt$outdir, split='/', fixed=TRUE))) == 1){
  # if only output foldername is provided use current working directory as basepath
  outdir <- file.path(getwd(), opt$outdir)
  dir.create(outdir, showWarnings = FALSE)
} else {
  # full path is provided
  outdir <- file.path(opt$outdir)
  dir.create(outdir, showWarnings = FALSE)
}

## Check provided input directory <opt$indir>

if (dir.exists(file.path(opt$indir))){ #if (!is.null(opt$indir)){
  if (isAbsolutePath(opt$indir)){
    # full path to existing directory is provided
    indir <- file.path(opt$indir)
  } else {
    # only subdirectory is provided, but full path exists
    indir <- file.path(getwd(), opt$indir)
  }
  # get pdb files in list
  files <- list.files(path = indir, pattern = "*.pdb", full.names = T, recursive = F)
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}

#print(indir)
#print(outdir)
setwd(outdir)
abs_output_path = getAbsolutePath(".")
#print(abs_output_path)

##-------------------------------------
# Protect spaces in path names with gsub(" ","\\\\ ",pathname)
scriptpath = gsub(" ","\\\\ ",paste(projectdir,'/src/analysis_bio3d_reporting.py', sep=''))
str_input_path = gsub(" ","\\\\ ",opt$indir)
output_path = gsub(" ","\\\\ ",abs_output_path)


## The actual program...
#--------------------------------------------------------------------

# loading pdb files
#pdbs <- pdbaln(files)
pdbs <- pdbaln(files, exefile='msa')

ids <- sub("[.].*", "", basename(pdbs$id)) # get filenames and drop any extensions
#ids <- sapply(strsplit(basename(pdbs$id), "[.]"), head, 1)
#ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))
#ids <- unlist(substr(basename(pdbs$id), 1, 7))
print(ids)

# Plot the Multiple Sequence Alignment
png("alignment_overview.png", units="in", width=5, height=5, res=300)
plot(pdbs)
dev.off()
print("Plot saved to file alignment_overview.png")


resno = pdbs$resno[1, !is.gap(pdbs)]
resid = aa123(pdbs$ali[1, !is.gap(pdbs)])


### Set number_of_groups for clustering
number_of_groups <- opt$ngroups


## RMSD
##-------------------------------------
rd <- rmsd(pdbs, fit=FALSE)
png("RMSD_hist.png", units="in", width=5, height=5, res=300)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of backbone RMSDs")
dev.off()
print("Plot saved to file RMSD_hist.png")

## RMSD heatmap
png("RMSD_heatmap.png", units="in", width=5, height=5, res=300)
#heatmap(rd, labCol=ids, main="RMSD Heatmap")
pheatmap(rd, main="Backbone RMSD Heatmap", fontsize = 6, show_colnames = FALSE) #annotation_row = ids
dev.off()
print("Plot saved to file RMSD_heatmap.png")

## RMSD Hierarchical clustering
hc_rmsd <- hclust(as.dist(rd))
png("RMSD_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_rmsd, labels=ids, cex=0.5, k=number_of_groups,
           ylab="RMSD (Å)", main="Backbone RMSD Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file RMSD_dendrogram.png")

## k-means clustering in RMSD space
#---------
# Choice of number_of_groups - Elbow method
if (length(files)>=5) {
    wss <- NULL
    # kmeans clustering
    for (i in 1:5) wss[i] <- sum(kmeans(rd, centers=i)$withinss)
    png(filename="cluster_number_elbow_plot_kmeans_onRMSD.png", width=750, height=750, units="px", res=120)
    plot(1:5, wss, type="b", xlab="Number of Clusters", ylab="Within-clusters sum of squares",
        main="Within-cluster dissimilarity based on k-means clustering")
    dev.off()
}
set.seed(11)
kmeans_rmsd <- kmeans(rd, centers=number_of_groups, nstart=5) # use the nstart argument for kmeans() to force it make a number of random starts and pick the best clustering in the end:
# Visualization of clusters
#plotting them in the 2D space of the scaled distance matrix
ddim_rmsd <- cmdscale(dist(rd), k=2)
# clusplot(ddim_rmsd, kmeans_rmsd$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
#          main="Bivariate k-means cluster plot\nof the scaled RMSD matrix")
png(filename="RMSD_bivariate_kmeans_cluster_plot.png", width=750, height=750, units="px", res=120)
clusplot(ddim_rmsd, kmeans_rmsd$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
         main="Bivariate k-means cluster plot\nof the scaled RMSD matrix (on backbone coordinates)", sub = NA)
dev.off()
# shade: the ellipses are shaded in relation to their density. The density is the number of points in the cluster divided by the area of the ellipse.

## Cluster validation metrics
if (length(files)>=5) {
    cluster_validation_test_rmsd <- clValid(rd, c(2:5), # cluster sizes 2,3,4,5
                      clMethods = c("hierarchical", "kmeans",  "pam" ),
                      validation = c("internal"), # "stability"
                      maxitems = nrow(rd) # needed for very large datasets
    )
    summary(cluster_validation_test_rmsd)
    writeLines(capture.output(summary(cluster_validation_test_rmsd)), paste0(output_path, "/cluster_validation_RMSD.txt"))
    print(paste0("Cluster validation statistics saved to file cluster_validation_RMSD.txt"))
}


## RMSF
##-------------------------------------
## Ignore gap containing positions
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)

## Tailor the PDB structure to exclude gap positions for SSE annotation
print(files[1])
pdb1 <- read.pdb(files[1])
pdb1$id <- files[1]
pdb2 <- read.pdb(files[length(files)])
pdb2$id <- files[length(files)]

ref_pdb <- trim.pdb(pdb1, inds = atom.select(pdb1, resno = pdbs$resno[pdb1$id, gaps.res$f.inds]))

residue_vec <- pdbs$resno[pdb1$id, gaps.res$f.inds]
residue_vec_with_gaps <- pdbs$resno[pdb1$id]

data_on_structure <- function(pdbX, residue_vec, data_vec, dataname){
    dataframe <- cbind(data.frame(residue_vec), data_vec)
    #pdbX <- read.pdb(pdbfile)
    # set all b-factors to 0
    pdbX$atom$b <- 0

    for(i in 1:length(dataframe[,1])){
        # get residue number from first column in dataframe
        residue_num <- dataframe[i,1]
        #print(residue_num)
        # get data value from provided column in dataframe
        value <- dataframe[i,2]
        #print(value)
        # get indices of residue
        residue_inds <- atom.select(pdb1, resno=c(residue_num))
        #print(residue_inds)
        # store data value in b-factor column for each residue in PDB file
        #print(pdb1$atom$b[ residue_inds$atom ])
        pdbX$atom$b[ residue_inds$atom ] <- value
    }

    # write to file
    write.pdb(pdbX, file=paste0(dataname,"_data_on_structure.pdb"))
}

## Plot RMSF with SSE annotation and labeled with residue numbers
rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
png("RMSF.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf, resno=ref_pdb, sse=ref_pdb, ylab="RMSF (Å)",
           xlab="Residue No.", col="gray", main="Backbone RMSFs per Residue") #, typ="l"
dev.off()
print("Plot saved to file RMSF.png")

rf_withgaps <- rmsf(pdbs$xyz)
png("RMSF_including_gaps.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf_withgaps, rm.gaps=FALSE, resno=pdb1, sse=pdb1, sse.min.length=0, ylab="RMSF (Å)",
           xlab="Residue No.", col="gray", main="Backbone RMSFs per Residue") #, typ="l"
dev.off()
print("Plot saved to file RMSF_including_gaps.png")

# Save reference structure with RMSF in B-factor column
ca.ref_pdb <- trim.pdb(ref_pdb, "calpha")
tryCatch(
    #try to do this
    {
    write.pdb(ca.ref_pdb, b=rf, file="RMSFonReference.pdb")
    print("PDB saved to file RMSFonReference.pdb")
    },
    #if an error occurs, tell me the error
    error=function(e) {
        message('An Error Occurred. This is probably related to the reference structure.')
        print(e)
    }
)

# save data in b-factor column on structure
data_on_structure(pdb1, residue_vec, rf, "RMSF")


## B-factors
##-------------------------------------
png("B-factors.png", units="in", width=5, height=5, res=300)
plot.bio3d(pdbs$b, rm.gaps=TRUE, ylab="B-factor", #, resno=ref_pdb, sse=ref_pdb
           xlab="Residue No.", main="Backbone B-factors", col="gray") #, typ="l"
dev.off()
print("Plot saved to file B-factors.png")



## PCA
##-------------------------------------
# PCA on coordinates (backbone)
pc_xyz <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])
pc_xyz

## Residue contribution to PCA
png("PCA_residue_contribution.png", units="in", width=5, height=5, res=300)
par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc_xyz$au[,1], resno=ref_pdb, sse=ref_pdb, ylab="PC1")
plot.bio3d(pc_xyz$au[,2], resno=ref_pdb, sse=ref_pdb, ylab="PC2")
plot.bio3d(pc_xyz$au[,3], resno=ref_pdb, sse=ref_pdb, ylab="PC3")
dev.off()
print("Plot saved to file PCA_residue_contribution.png")

## Interpolated structures along PC1/2/3 produced by the mktrj.pca() function
mktrj.pca(pc_xyz, pc=1, resno=resno, resid=resid, file="PC1.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz, pc=2, resno=resno, resid=resid, file="PC2.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz, pc=3, resno=resno, resid=resid, file="PC3.pdb") #, mag = 1, step = 0.125
print("Interpolated trajectory structures saved to files PC1.pdb, PC2.pdb, PC3.pdb")

# PC1 Vector field representation
pymol(pc_xyz, pdb=pdb1, pc=1, resno=resno, resid=resid, as="ribbon", file="PC1vectors.pml", type="script")

### Hierarchical clustering in PC space
# Perform structural clustering in the PC1-PC2 subspace.
hc_pc12 <- hclust(dist(pc_xyz$z[, 1:2]))
grps_pc12 <- cutree(hc_pc12, k=number_of_groups)

# Plot PCs
png("PCA.png", units="in", width=5, height=5, res=300)
plot(pc_xyz, col=grps_pc12) #, col=annotation[, "color"]
dev.off()
print("Plot saved to file PCA.png")

# Plot PCA dendrogram
png("PCA_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="Backbone PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_dendrogram.png")


# Save PyMOL scripts for coloring by group
grps_rmsd <- cutree(hc_rmsd, k=number_of_groups)
# color by clustering (based on RMSD)
pymol(pdbs, col=grps_rmsd, as="cartoon", file="col_by_grps_RMSD.pml", type="script")
# color by clustering (based on PCA)
pymol(pdbs, col=grps_pc12, as="cartoon", file="col_by_grps_PC.pml", type="script")


## k-means clustering in PC space
#---------
if (length(files)>=5) {
    wss <- NULL
    # kmeans clustering
    for (i in 1:5) wss[i] <- sum(kmeans(pc_xyz$z[,1:2], centers=i)$withinss)
    png(filename="cluster_number_elbow_plot_kmeans_onPCA.png", width=750, height=750, units="px", res=120)
    plot(1:5, wss, type="b", xlab="Number of Clusters", ylab="Within-clusters sum of squares",
        main="Within-cluster dissimilarity based on k-means clustering of PC1-2")
    dev.off()
}
set.seed(11)
kmeans_pc_xyz <- kmeans(pc_xyz$z[,1:2], centers=number_of_groups, nstart=5) # use the nstart argument for kmeans() to force it make a number of random starts and pick the best clustering in the end:
# Visualization of clusters
#plotting them in the 2D space of the scaled distance matrix
# ddim_pc_xyz <- cmdscale(dist(pc_xyz$z[,1:2]), k=2)
# png(filename="PCA_bivariate_kmeans_cluster_plot.png", width=750, height=750, units="px", res=120)
# clusplot(ddim_pc_xyz, kmeans_pc_xyz$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
#          main="Bivariate k-means cluster plot\nin scaled PC1-2 space (based on backbone coordinates)")
png(filename="PCA_kmeans_cluster_plot.png", width=750, height=750, units="px", res=120)
clusplot(pc_xyz$z[,1:2], kmeans_pc_xyz$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
         main="Bivariate k-means cluster plot\nin PC1-2 space (based on backbone coordinates)")
dev.off()

## Cluster validation metrics
if (length(files)>=5) {
    cluster_validation_test_pc_xyz <- clValid(pc_xyz$z[,1:2], c(2:5), # cluster sizes 2,3,4,5
                          clMethods = c("hierarchical", "kmeans",  "pam" ),
                          validation = c("internal"), # "stability"
                          maxitems = nrow(rd) # needed for very large datasets
    )
    summary(cluster_validation_test_pc_xyz)
    writeLines(capture.output(summary(cluster_validation_test_pc_xyz)), paste0(output_path, "/cluster_validation_PCA.txt"))
    print(paste0("Cluster validation statistics saved to file cluster_validation_PCA.txt"))
}

## Torsion/Dihedral analysis
##-------------------------------------
# Calculate torsion/dihedral angles
# atm.inc: a numeric value indicating the number of atoms to increment by between successive torsion evaluations
tor <- t(apply( pdbs$xyz[, gaps.pos$f.inds], 1, torsion.xyz, atm.inc=3))

#-- PCA on torsion data
pc_tor <- pca.tor(tor)

# Perform structural clustering in the PC1-PC2 subspace.
hc_pc12_tor <- hclust(dist(pc_tor$z[, 1:2]))
grps_pc12_tor <- cutree(hc_pc12_tor, k=number_of_groups)

# Plot PCs
png("PCA_on_Torsion.png", units="in", width=5, height=5, res=300)
plot.pca(pc_tor, col=grps_pc12_tor)
dev.off()
print("Plot saved to file PCA_on_Torsion.png")

png("PCA_on_Torsion_loadings.png", units="in", width=5, height=5, res=300)
plot.pca.loadings(pc_tor)
dev.off()
print("Plot saved to file PCA_on_Torsion_loadings.png")

png("PCA_on_Torsion_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12_tor, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="Backbone Torsion PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_on_Torsion_dendrogram.png")



# Construct a Contact Map for Structures
##-------------------------------------
# Calculate and color by averaged contact density around each residue
# binary=FALSE: the raw matrix containing fraction of frames that two residues are in contact is returned
cm <- cmap(pdbs, binary=FALSE, all.atom=FALSE)

png(filename="contact_map.png", width=900, height=750, units="px", res=120)
plot.cmap(cm)
dev.off()
print("Plot saved to file contact_map.png")

vec <- rowSums(cm, na.rm=TRUE)
pymol(pdbs, col="user", user.vec=vec, as="cartoon", file="col_by_averaged_contact_density.pml", type="script")

png(filename="averaged_contact_density.png", width=900, height=750, units="px", res=120)
plot.bio3d(vec, resno=ref_pdb, sse=ref_pdb, ylab="Averaged Contact Density",
           xlab="Residue No.", col="gray", main="Averaged Contact Density per Residue") #, typ="l"
dev.off()
print("Plot saved to file averaged_contact_density.png")

# save data in b-factor column on structure
data_on_structure(pdb1, residue_vec_with_gaps, vec, "averaged_contact_density")


# UMAP Analysis
##-------------------------------------
# UMAP, short for “Uniform Manifold Approximation and Projection” is a non-linear dimension reduction technique.
# Important parameters are n_neighbors (15) and min_dist (0.1)
# n_components defines the number of dimensions to project to, default is 2, as needed for 2D plots
# umap configuration parameters: print(umap.defaults)
#            n_neighbors: 15
#           n_components: 2
#                 metric: euclidean
#               n_epochs: 200
#                  input: data
#                   init: spectral
#               min_dist: 0.1
#       set_op_mix_ratio: 1
#     local_connectivity: 1
#              bandwidth: 1
#                  alpha: 1
#                  gamma: 1
#   negative_sample_rate: 5
#                      a: NA
#                      b: NA
#                 spread: 1
#           random_state: NA
#        transform_state: NA
#                    knn: NA
#            knn_repeats: 1
#                verbose: FALSE
#        umap_learn_args: NA

# create a new settings object with n_neighbors set to 3 to be able to run it with just 3 structures
custom.settings = umap.defaults
custom.settings$n_neighbors = 3
custom.settings$random_state = 123 # for reproducibility

# UMAP on coordinates (backbone)
# returns: umap embedding of N items in 2 dimensions;
# object components: layout, data, knn, config
umap_fit <- umap(pdbs$xyz[, gaps.pos$f.inds], config=custom.settings)
print(umap_fit)
#print(umap_fit$layout) # contains the 2D projections
#print(umap_fit$data) # contains the initially supplied data matrix
#print(umap_fit$knn) # object components: indexes, distances
#print(umap_fit$knn$indexes) # NxN matrix
#print(umap_fit$knn$distances) # NxN matrix

### Hierarchical clustering in UMAP space
# Perform structural clustering in the UMAP1-2 space.
hc_umap <- hclust(dist(umap_fit$layout))
grps_umap <- cutree(hc_umap, k=number_of_groups)

# Plot UMAP
png("UMAP.png", units="in", width=5, height=5, res=300)
plot(umap_fit$layout, col=grps_umap, xlab="UMAP1", ylab="UMAP2", main="UMAP on backbone coordinates") #, col=annotation[, "color"]
dev.off()
print("Plot saved to file UMAP.png")

png("UMAP_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_umap, labels=ids, cex=0.5, k=number_of_groups,
           ylab="UMAP1-2 distance", main="UMAP Cluster Dendrogram (backbone coordinates)", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file UMAP_dendrogram.png")

# umap_df <- umap_fit$layout %>%
#   as.data.frame()%>%
#   rename(UMAP1="V1",
#          UMAP2="V2")
# umap_df %>%
#   ggplot(aes(x = UMAP1,
#              y = UMAP2))+
#   geom_point()+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#       subtitle = "UMAP plot")
# ggsave("UMAP_plot.png")

## k-means clustering in UMAP space
#---------
if (length(files)>=5) {
    wss <- NULL
    # kmeans clustering
    for (i in 1:5) wss[i] <- sum(kmeans(umap_fit$layout, centers=i)$withinss)
    png(filename="cluster_number_elbow_plot_kmeans_onUMAP.png", width=750, height=750, units="px", res=120)
    plot(1:5, wss, type="b", xlab="Number of Clusters", ylab="Within-clusters sum of squares",
        main="Within-cluster dissimilarity based on k-means clustering on backbone UMAP")
    dev.off()
}
set.seed(11)
kmeans_umap_xyz <- kmeans(umap_fit$layout, centers=number_of_groups, nstart=5) # use the nstart argument for kmeans() to force it make a number of random starts and pick the best clustering in the end:
# Visualization of clusters
#plotting them in the 2D space of the scaled distance matrix
# ddim_umap_xyz <- cmdscale(dist(umap_fit$layout), k=2)
# png(filename="UMAP_bivariate_kmeans_cluster_plot.png", width=750, height=750, units="px", res=120)
# clusplot(ddim_umap_xyz, kmeans_umap_xyz$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
#          main="Bivariate k-means cluster plot\nin scaled 2D UMAP space (based on backbone coordinates)")
png(filename="UMAP_kmeans_cluster_plot.png", width=750, height=750, units="px", res=120)
clusplot(umap_fit$layout, kmeans_umap_xyz$cluster, color=T, labels=4, lines=0, plotchar=T, shade=T,
         main="Bivariate k-means cluster plot\nin 2D UMAP space (based on backbone coordinates)")
dev.off()

## Cluster validation metrics
if (length(files)>=5) {
    cluster_validation_test_umap_xyz <- clValid(umap_fit$layout, c(2:5), # cluster sizes 2,3,4,5
                          clMethods = c("hierarchical", "kmeans",  "pam" ),
                          validation = c("internal"), # "stability"
                          maxitems = nrow(rd) # needed for very large datasets
    )
    summary(cluster_validation_test_umap_xyz)
    writeLines(capture.output(summary(cluster_validation_test_umap_xyz)), paste0(output_path, "/cluster_validation_UMAP.txt"))
    print(paste0("Cluster validation statistics saved to file cluster_validation_UMAP.txt"))
}



##-------------------------------------
## All-atom Analysis
##-------------------------------------

# read all atoms for each residue
#pdbs_allatoms <- read.all(pdbaln(files, exefile='msa'))
# update the aligned structures to include all heavy atoms:
pdbs_allatoms <- read.all(pdbs)
##-------------------------------------


## all-atom RMSD
##-------------------------------------
rd_allatom <- rmsd(pdbs_allatoms$all, fit=FALSE)
png("RMSD_hist_allatom.png", units="in", width=5, height=5, res=300)
hist(rd_allatom, breaks=40, xlab="RMSD (Å)", main="Histogram of all-atom RMSDs")
dev.off()
print("Plot saved to file RMSD_hist_allatom.png")

## all-atom RMSD heatmap
png("RMSD_heatmap_allatom.png", units="in", width=5, height=5, res=300)
#heatmap(rd_allatom, labCol=ids, main="RMSD Heatmap")
pheatmap(rd_allatom, main="All-atom RMSD Heatmap", fontsize = 6, show_colnames = FALSE) #annotation_row = ids
dev.off()
print("Plot saved to file RMSD_heatmap_allatom.png")

## all-atom RMSD Hierarchical clustering
hc_rmsd_allatom <- hclust(as.dist(rd_allatom))
png("RMSD_dendrogram_allatom.png", units="in", width=5, height=5, res=300)
hclustplot(hc_rmsd_allatom, labels=ids, cex=0.5, k=number_of_groups,
           ylab="RMSD (Å)", main="All-atom RMSD Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file RMSD_dendrogram_allatom.png")


# all-atom RMSF
##-------------------------------------
rf_allatom <- rmsf(pdbs_allatoms$all[, gaps.pos$f.inds])
png(filename="RMSF_allatom.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf_allatom, resno=ref_pdb, sse=ref_pdb, ylab="RMSF (Å)",
           xlab="Residue No.", col="gray", main="All-atom RMSFs per Residue") #, typ="l"
dev.off()
print("Plot saved to file RMSF_allatom.png")

# Save reference structure with RMSF in B-factor column
tryCatch(
    #try to do this
    {
    write.pdb(ca.ref_pdb, b=rf_allatom, file="RMSFonReference_allatom.pdb")
    print("PDB saved to file RMSFonReference_allatom.pdb")
    },
    #if an error occurs, tell me the error
    error=function(e) {
        message('An Error Occurred. This is probably related to the reference structure.')
        print(e)
    }
)

# save data in b-factor column on structure
data_on_structure(pdb1, residue_vec, rf_allatom, "RMSF_allatom")


## PCA on all-atom coordinates
##-------------------------------------
# Perform PCA on all-atom coordinates
# In input xyz (MxN),  N > 3000 and M < N
# Singular Value Decomposition (SVD) approach is faster and is recommended (set 'use.svd = TRUE')
# Here (use.svd=TRUE) singular value decomposition (SVD) is called instead of eigenvalue decomposition
pc_xyz_allatom <- pca(pdbs_allatoms$all, use.svd = TRUE, rm.gaps=TRUE, fit=FALSE)
pc_xyz_allatom

# ## Residue contribution to all-atom PCA - not possible to calculate straight away, as each residue has different number of atoms
# png("PCA_residue_contribution_allatom.png", units="in", width=5, height=5, res=300)
# par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
# plot.bio3d(pc_xyz_allatom$au[,1], resno=pdbs_allatoms[1], sse=pdbs_allatoms[1], ylab="PC1")
# plot.bio3d(pc_xyz_allatom$au[,2], resno=pdbs_allatoms[1], sse=pdbs_allatoms[1], ylab="PC2")
# plot.bio3d(pc_xyz_allatom$au[,3], resno=pdbs_allatoms[1], sse=pdbs_allatoms[1], ylab="PC3")
# dev.off()
# print("Plot saved to file PCA_residue_contribution_allatom.png")

## Interpolated structures along PC1/2/3 produced by the mktrj.pca() function
mktrj.pca(pc_xyz_allatom, pc=1, file="PC1_allatom.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz_allatom, pc=2, file="PC2_allatom.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz_allatom, pc=3, file="PC3_allatom.pdb") #, mag = 1, step = 0.125
print("Interpolated trajectory structures saved to files PC1_allatom.pdb, PC2_allatom.pdb, PC3_allatom.pdb")

# PC1-allatom Vector field representation
pymol(pc_xyz_allatom, pdb=pdbs_allatoms[1]$all, pc=1, as="lines", file="PC1vectors_allatom.pml", type="script")

# Perform structural clustering in the PC1-PC2 subspace (for allatom data).
hc_pc12_allatom <- hclust(dist(pc_xyz_allatom$z[, 1:2]))
grps_pc12_allatom <- cutree(hc_pc12_allatom, k=number_of_groups)

# Plot PCA
png("PCA_allatom.png", units="in", width=5, height=5, res=300)
plot(pc_xyz_allatom, col=grps_pc12_allatom)
dev.off()
print("Plot saved to file PCA_allatom.png")

# Plot PC loadings
png("PCA_loadings_allatom.png", units="in", width=5, height=5, res=300)
plot.pca.loadings(pc_xyz_allatom)
dev.off()
print("Plot saved to file PCA_loadings_allatom.png")

# Plot PCA dendrogram
png("PCA_dendrogram_allatom.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12_allatom, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="All-atom PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_dendrogram_allatom.png")

# Save PyMOL scripts for coloring by group
grps_rmsd_allatom <- cutree(hc_rmsd_allatom, k=number_of_groups)
# color by clustering (based on RMSD on all-atom)
pymol(pdbs_allatoms, col=grps_rmsd_allatom, as="cartoon", file="col_by_grps_RMSD_allatom.pml", type="script")
# color by clustering (based on PCA on all-atom)
pymol(pdbs_allatoms, col=grps_pc12_allatom, as="cartoon", file="col_by_grps_PC_allatom.pml", type="script")


# Construct a Contact Map for Structures
##-------------------------------------
# Calculate and color by averaged contact density around each residue
# binary=FALSE: the raw matrix containing fraction of frames that two residues are in contact is returned
cm_allatom <- cmap(pdbs_allatoms, binary=FALSE, all.atom=TRUE)

png(filename="contact_map_allatom.png", width=900, height=750, units="px", res=120)
plot.cmap(cm_allatom)
dev.off()
print("Plot saved to file contact_map_allatom.png")

vec_allatom <- rowSums(cm_allatom, na.rm=TRUE)
pymol(pdbs, col="user", user.vec=vec_allatom, as="cartoon", file="col_by_averaged_contact_density_allatom.pml", type="script")

png(filename="averaged_contact_density_allatom.png", width=900, height=750, units="px", res=120)
plot.bio3d(vec_allatom, resno=ref_pdb, sse=ref_pdb, ylab="Averaged Contact Density",
           xlab="Residue No.", col="gray", main="All-atom Averaged Contact Density per Residue") #, typ="l"
dev.off()
print("Plot saved to file averaged_contact_density_allatom.png")

# save data in b-factor column on structure
data_on_structure(pdb1, residue_vec_with_gaps, vec_allatom, "averaged_contact_density_allatom")

# # UMAP on all-atom coordinates - NOT working!
# ##-------------------------------------
# umap_fit_allatom <- umap(pdbs_allatoms$all, config=custom.settings)
# print(umap_fit_allatom)
#
# ### Hierarchical clustering in UMAP space
# # Perform structural clustering in the UMAP1-2 space.
# hc_umap_allatom <- hclust(dist(umap_fit_allatom$layout))
# grps_umap_allatom <- cutree(hc_umap_allatom, k=number_of_groups)
#
# # Plot UMAP
# png("UMAP_allatom.png", units="in", width=5, height=5, res=300)
# plot(umap_fit_allatom$layout, col=grps_umap_allatom, xlab="UMAP1", ylab="UMAP2", main="UMAP plot") #, col=annotation[, "color"]
# dev.off()
# print("Plot saved to file UMAP_allatom.png")
#
# png("UMAP_dendrogram_allatom.png", units="in", width=5, height=5, res=300)
# hclustplot(hc_umap_allatom, labels=ids, cex=0.5, k=number_of_groups,
#            ylab="UMAP1-2 distance", main="UMAP Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
# dev.off()
# print("Plot saved to file UMAP_dendrogram_allatom.png")


tryCatch(
    #try to do this
    {
    ## Ensemble Difference Distance Matrix (eDDM) Analysis
    ##-------------------------------------
    # The eDDM analysis compares structural ensembles under distinct ligation, activation, etc. conditions.
    # At least two groups of structures are required.
    # we use PCA of the distance matrices (through the bio3d function, pca.array()) followed by a conventional hierarchical
    # clustering in the PC1-PC2 subspace to get the intrinsic grouping of the structures.
    library(bio3d.eddm)
    eDDM_analysis <- TRUE
    # need to update the aligned structures to include all heavy atoms:
    #pdbs_allatoms <- read.all(pdbs) # This is already done at the beginning

    # Calculate distance matrices.
    # The option `all.atom=TRUE` tells that all heavy-atom coordinates will be used.
    dm <- dm(pdbs_allatoms, all.atom=TRUE)
    # Perform PCA of distance matrices.
    pc_dm <- pca.array(dm)

    # Perform structural clustering in the PC1-PC2 subspace.
    hc_dm_pc12 <- hclust(dist(pc_dm$z[, 1:2]))
    grps_dm_pc12 <- cutree(hc_dm_pc12, k=number_of_groups)

    # Plot PCs
    png("PCA_on_allatom_DifferenceDistanceMatrix.png", units="in", width=5, height=5, res=300)
    plot(pc_dm, col=grps_dm_pc12) # for only PC1-2 subplot add:, pc.axes=c(1,2)
    dev.off()
    print("Plot saved to file PCA_on_allatom_DifferenceDistanceMatrix.png")
    # Plot PC loadings
    png("PCA_on_allatom_DifferenceDistanceMatrix_loadings.png", units="in", width=5, height=5, res=300)
    plot.pca.loadings(pc_dm)
    dev.off()
    print("Plot saved to file PCA_on_allatom_DifferenceDistanceMatrix_loadings.png")
    # Plot PCA dendrogram
    png("PCA_on_allatom_DifferenceDistanceMatrix_dendrogram.png", units="in", width=5, height=5, res=300)
    hclustplot(hc_dm_pc12, labels=ids, cex=0.5, k=number_of_groups,
               ylab="PC1-2 distance", main="Difference Distance Matrix PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
    dev.off()
    print("Plot saved to file PCA_on_allatom_DifferenceDistanceMatrix_dendrogram.png")

    # Calculating difference distance matrices between groups
    # the difference mean distance between groups for each residue pair are calculated and statistical significance assessed
    # using a two-sample Wilcoxon test. Long-distance pairs in all structures are omitted.
    tbl <- eddm(pdbs_allatoms, grps=grps_dm_pc12, dm=dm, mask="cmap")
    # Plot
    png("eDDM_complete.png", units="in", width=5, height=5, res=300)
    plot(tbl, pdbs=pdbs_allatoms, full=TRUE, resno=NULL, sse=pdbs$sse[1, ], type="tile") #,labels=TRUE, labels.ind=c(1:3, 17:19)
    dev.off()
    print("Plot saved to file eDDM_complete.png")
    # Generate PyMol script to visualise all residue pairs showing any distance changes
    pymol(tbl, pdbs=pdbs_allatoms, grps=grps_dm_pc12, as="sticks", file="eDDM_any.pml", type="script")
    # Identifying significant distance changes
    #keys <- subset.eddm(tbl, alpha=0.005, beta=1.0, switch.only=TRUE)
    keys <- subset.eddm(tbl, alpha=0.005, beta=0.8, switch.only=FALSE)
    if (length(keys)) {
        # Plot
        png("eDDM_significant.png", units="in", width=5, height=5, res=300)
        plot(keys, pdbs=pdbs_allatoms, full=TRUE, resno=NULL, sse=pdbs$sse[1, ], type="tile") #,labels=TRUE, labels.ind=c(1:3, 17:19)
        dev.off()
        print("Plot saved to file eDDM_significant.png")
        # Generate PyMol script to visualise all identified key residue pairs showing significant distance changes
        tryCatch(
            #try to do this
            {
            pymol(keys, pdbs=pdbs_allatoms, grps=grps_dm_pc12, as="sticks", file="eDDM_significant.pml", type="script")
            },
            #if an error occurs, tell me the error
            error=function(e) {
                print(e)
            }
        )
    } else {print("No key switch residues identified.")}

    },
    #if an error occurs, tell me the error
    error=function(e) {
        message('R package "bio3d.eddm" is not installed. Proceeding without eDDM analysis.')
        print(e)
        eDDM_analysis <<- FALSE
    }
)

# ### Two-structure comparisons
# ###-------------------------------------
#
# # Locate the two structures in pdbs
# #print(files[1], files[length(files)])
# ind.a <- grep(pdb1$id, pdbs$id)
# ind.b <- grep(pdb2$id, pdbs$id)
#
# # Exclude gaps in the two structures to make them comparable
# gaps.xyz2 <- gap.inspect(pdbs$xyz[c(ind.a, ind.b), ])
# a.xyz <- pdbs$xyz[ind.a, gaps.xyz2$f.inds]
# b.xyz <- pdbs$xyz[ind.b, gaps.xyz2$f.inds]
#
# ## Difference distance matrix analysis (DDM)
# ##-------------------------------------
# a <- dm.xyz(a.xyz)
# b <- dm.xyz(b.xyz)
# # Calculate DDM
# ddm <- a - b
# # Plot DDM
# png("DDM.png", units="in", width=5, height=5, res=300)
# plot.dmat( ddm, nlevels=10, grid.col="gray", xlab=basename(pdb1$id), ylab=basename(pdb2$id), main="Difference distance matrix")
# dev.off()
#
# ## Torsion/Dihedral analysis
# ##-------------------------------------
# # Compare CA based pseudo-torsion angles between the two structures
# a <- torsion.xyz(a.xyz, atm.inc=1)
# b <- torsion.xyz(b.xyz, atm.inc=1)
# d.ab <- wrap.tor(a-b)
# d.ab[is.na(d.ab)] <- 0
#
# # Plot results with SSE annotation
# png("Ca_torsion_diff_between_first&last_structure.png", units="in", width=5, height=5, res=300)
# plot.bio3d(abs(d.ab), resno=pdb1, sse=pdb1, typ="h", xlab="Residue No.",
#            ylab = "Difference Angle", main="Ca torsion difference - first&last structure")
# dev.off()


# # Calculate pair-wise RMSD values
# rmsd.map <- rmsd(pdbs$xyz, a.inds=gaps.pos$f.inds, fit=TRUE)
# png("pair-wise_RMSDs.png", units="in", width=5, height=5, res=300)
# heatmap(rmsd.map, labCol=ids, symm=TRUE, main="Pair-wise RMSDs") #, labRow=annotation[, "state"]
# dev.off()


##-------------------------------------
## Cluster attributions
##-------------------------------------
# store cluster attributions in dataframe
if (eDDM_analysis==TRUE) {
    clusters_df <- as.data.frame(list(ids, grps_rmsd, grps_pc12, grps_umap, grps_pc12_tor,
                                 grps_rmsd_allatom, grps_pc12_allatom, grps_dm_pc12), #, grps_umap_allatom
                                 col.names = c("Structure", "backbone_RMSD", "backbone_PCA_onCoords",
                                 "backbone_UMAP_onCoords", "backbone_PCA_onTorsion", "allatom_RMSD",
                                 "allatom_PCA_onCoords", "allatom_PCA_onDist")) #, "allatom_UMAP_onCoords"
} else {
    clusters_df <- as.data.frame(list(ids, grps_rmsd, grps_pc12, grps_umap, grps_pc12_tor,
                                 grps_rmsd_allatom, grps_pc12_allatom), #, grps_dm_pc12, grps_umap_allatom
                                 col.names = c("Structure", "backbone_RMSD", "backbone_PCA_onCoords",
                                 "backbone_UMAP_onCoords", "backbone_PCA_onTorsion", "allatom_RMSD",
                                 "allatom_PCA_onCoords")) #, "allatom_PCA_onDist", "allatom_UMAP_onCoords"
}

# convert first column to rownames
clusters_df <- data.frame(clusters_df, row.names = 1)

# safe dataframe
write.csv(clusters_df, "cluster_attributions.csv", row.names=TRUE, quote=FALSE)

# library(RColorBrewer)
# default_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))#(100)
# mycolors <- default_color(length(unique(clusters_consensus_df$Consensus_Cluster)))
# names(mycolors) <- unique(clusters_consensus_df$Consensus_Cluster)
# mycolors <- list(mycolors = mycolors)

png("cluster_attributions_heatmap.png", units="in", width=5, height=5, res=300)
pheatmap(clusters_df, main="Cluster Attributions",
        fontsize = 6,
        angle_col = 45,
        cutree_rows = number_of_groups,
        legend_breaks = c(1:number_of_groups),
        #annotation_row = clusters_consensus_df["Consensus_Cluster"],
        #annotation_colors = mycolors,
        #annotation_legend = FALSE
        )
dev.off()
print("Plot saved to file cluster_attributions_heatmap.png")

#clusters_consensus_df <- clusters_df
# get consensus cluster attribution (most common) and add to df
#clusters_consensus_df["Consensus_Cluster"] <- apply(clusters_consensus_df, 1, function(x) names(which.max(table(x))))


## Finding Consensus Cluster with co-assignment matrix

# Compute co-assignment matrix from clusters_df
compute_coassignment_matrix <- function(clusters_df) {
  n <- nrow(clusters_df)
  coassign_matrix <- matrix(0, n, n)
  rownames(coassign_matrix) <- rownames(clusters_df)
  colnames(coassign_matrix) <- rownames(clusters_df)

  for (method in colnames(clusters_df)) {
    clustering <- clusters_df[[method]]
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if (clustering[i] == clustering[j]) {
          coassign_matrix[i, j] <- coassign_matrix[i, j] + 1
          coassign_matrix[j, i] <- coassign_matrix[i, j]
        }
      }
    }
  }

  coassign_matrix <- coassign_matrix / ncol(clusters_df)  # Normalize by number of methods
  diag(coassign_matrix) <- 1  # Each entry is always in the same cluster as itself
  return(coassign_matrix)
}

# Compute co-assignment matrix
coassign_matrix <- compute_coassignment_matrix(clusters_df)

# Perform hierarchical clustering on co-assignment matrix to get consensus clustering
hc <- hclust(as.dist(1 - coassign_matrix), method = "average")  # 1 - similarity for distance

# Cut tree into predefined number of groups
consensus_clusters <- cutree(hc, k = number_of_groups)

# Convert to dataframe for annotation
consensus_df <- data.frame(Consensus_Cluster = factor(consensus_clusters))
rownames(consensus_df) <- rownames(clusters_df)

# safe dataframe
write.csv(consensus_df, "consensus_cluster.csv", row.names=TRUE, quote=FALSE)

# Save and plot heatmap of co-assignment matrix
png("coassignment_heatmap.png", units="in", width=8, height=5, res=300)
pheatmap(coassign_matrix,
         main = "Co-assignment Matrix",
         fontsize = 6,
         angle_col = 45,
         cluster_rows = hc,
         cluster_cols = hc,
         annotation_row = consensus_df,
         annotation_col = consensus_df,
         show_rownames = TRUE,  # Show structure names
         show_colnames = TRUE,
         color = colorRampPalette(c("white", "black"))(100))
dev.off()
print("Plot saved to file coassignment_heatmap.png")




### Save R session info
##-------------------------------------
writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info_", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
#writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info.txt"))
print(paste0("R session info saved to file R_session_info_", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))

### Generate Report
##-------------------------------------
#system(paste('python3', scriptpath, str_input_path, output_path, number_of_groups, sep=' '), wait=FALSE)
