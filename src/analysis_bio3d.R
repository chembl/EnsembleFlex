#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

projectdir = getwd()
#print(projectdir)

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

library("optparse")
library("bio3d")
library("msa")

option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="input directory path", metavar="character"),
  # make_option(c("-f", "--filenames"), type="character", default=NULL,
  #             help="dataset file names - example: [file1.pdb,file2.pdb]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="Bio3D_Analysis",
              help="output directory [default=%default]", metavar="character"),
  make_option(c("-n", "--ngroups"), type="integer", default=3,
              help="number of groups for clustering [default=%default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (opt$outdir == "Bio3D_Analysis"){
  outdir <- file.path(opt$indir, opt$outdir)
  dir.create(outdir, showWarnings = FALSE)
} else {
  outdir <- file.path(opt$outdir)
  dir.create(outdir, showWarnings = FALSE)
}

setwd(outdir)
print(outdir)

if (!is.null(opt$indir)){
  setwd(outdir)
  files <- list.files(path = opt$indir, pattern = "*.pdb", full.names = T, recursive = F)
  # } else if (!is.null(opt$filenames)) {
  #   files <- as.list(strsplit(opt$filenames, ",")[[1]])
  #   #files <- opt$filenames
  #   setwd(paste(dirname(files[1]),opt$outdir))
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}


## The actual program...
#--------------------------------------------------------------------

# loading pdb files
# option 1.
#pdbs <- pdbaln(files)
pdbs <- pdbaln(files, exefile='msa')
ids <- sub("[.].*", "", basename(pdbs$id)) # get filenames and drop any extensions
#ids <- sapply(strsplit(basename(pdbs$id), "[.]"), head, 1)
#ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))
#ids <- unlist(substr(basename(pdbs$id), 1, 7))
print(ids)
# # option 2. if they are of identical composition and you only want xyz coordinates
# xyz <- NULL
# for(i in files) {
#   xyz = rbind(xyz, read.pdb(files)$xyz)
# }

### Set number_of_groups for clustering
number_of_groups <- opt$ngroups

## RMSD
##-------------------------------------
rd <- rmsd(pdbs, fit=FALSE)
png("RMSD_hist.png", units="in", width=5, height=5, res=300)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
dev.off()

## RMSD heatmap
library(pheatmap)
png("RMSD_heatmap.png", units="in", width=5, height=5, res=300)
#heatmap(rd, labCol=ids, main="RMSD Heatmap")
pheatmap(rd, main="RMSD Heatmap", fontsize = 6, show_colnames = FALSE) #annotation_row = ids
dev.off()

## RMSD Hierarchical clustering
hc_rmsd <- hclust(as.dist(rd))
png("RMSD_clust.png", units="in", width=5, height=5, res=300)
hclustplot(hc_rmsd, labels=ids, cex=0.5, k=number_of_groups,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()


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

## Plot RMSF with SSE annotation and labeled with residue numbers
rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
png("RMSF.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf, resno=ref_pdb, sse=ref_pdb, ylab="RMSF (Å)",
           xlab="Residue No.", main="RMSFs", col="gray") #, typ="l"
dev.off()

# Save reference structure with RMSF in B-factor column
ca.ref_pdb <- trim.pdb(ref_pdb, "calpha")
write.pdb(ca.ref_pdb, b=rf, file="RMSFonReference.pdb")


## B-factors
##-------------------------------------
png("B-factors.png", units="in", width=5, height=5, res=300)
plot.bio3d(pdbs$b, rm.gaps=TRUE, resno=ref_pdb, sse=ref_pdb, ylab="B-factor",
           xlab="Residue No.", main="B-factors", col="gray") #, typ="l"
dev.off()


## PCA
##-------------------------------------
# Do PCA on coordinates
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])
pc.xray

## Residue contribution to PCA
png("PCA_residue_contribution.png", units="in", width=5, height=5, res=300)
par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno=ref_pdb, sse=ref_pdb, ylab="PC1")
plot.bio3d(pc.xray$au[,2], resno=ref_pdb, sse=ref_pdb, ylab="PC2")
plot.bio3d(pc.xray$au[,3], resno=ref_pdb, sse=ref_pdb, ylab="PC3")
dev.off()

## Interpolated structures along PC1/2/3 produced by the mktrj.pca() function
mktrj.pca(pc.xray, pc=1, file="PC1.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc.xray, pc=2, file="PC2.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc.xray, pc=3, file="PC3.pdb") #, mag = 1, step = 0.125

# Vector field representation
pymol(pc.xray, pdb=ref_pdb, pc=1, as="cartoon", file="PC1vectors.pml", type="script")


### Hierarchical clustering in PC space
# Perform structural clustering in the PC1-PC2 subspace.
hc_pc12 <- hclust(dist(pc.xray$z[, 1:2]))
grps_pc12 <- cutree(hc_pc12, k=number_of_groups)
# Plot PCs
png("PCA.png", units="in", width=5, height=5, res=300)
plot(pc.xray, col=grps_pc12) #, col=annotation[, "color"]
dev.off()
png("PCA_clust.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()


## Cluster attributions
grps_rmsd <- cutree(hc_rmsd, k=number_of_groups)
# store cluster attributions in dataframe
clusters_df <- as.data.frame(list(ids, grps_rmsd, grps_pc12), col.names = c("Structure", "RMSD Cluster", "PC1/2 Cluster"))
# safe dataframe
write.csv(clusters_df, "cluster_attributions.csv", row.names=FALSE, quote=FALSE)


# Save PyMOL scripts for coloring by group
# color by clustering (based on PCA)
pymol(pdbs, col=grps_pc12, as="cartoon", file="col_by_grps_PC.pml", type="script")
# color by clustering (based on RMSD)
pymol(pdbs, col=grps_rmsd, as="cartoon", file="col_by_grps_RMSD.pml", type="script")



## eNMA
##-------------------------------------
# NMA on all structures; use 'rm.gaps=FALSE' to keep the gap containing columns, but note that this is not compatible
# with mktrj
modes <- nma.pdbs(pdbs, rm.gaps=TRUE)

# Make fluctuation plot
png("eNMA_fluctuations.png", units="in", width=5, height=5, res=300)
plot(modes, pdbs=pdbs, main="Normal Mode fluctuations") #, col=annotation[, "color"]
dev.off()

## Interpolated structures along eNMA modes produced by the mktrj() function
mktrj(enma = modes, pdbs = pdbs, mag = 10, step = 1.25, file = "eNMA.pdb", rock = TRUE)
# Vector field representation
#pymol(modes, pdb=ref_pdb, mode=7, file="eNMA_col_mode7.pml", type="script")


# NMA for single reference structure
##-------------------------------------
modes_ref_pdb <- nma(ref_pdb)
png(filename="NMA_fluctuations_reference_pdb.png", width=900, height=750, units="px", res=120)
plot.nma(modes_ref_pdb, resno=ref_pdb, sse=ref_pdb, sse.min.length=3)#, main="NMA on reference structure"
dev.off()
# Make a PDB trajectory
mktrj(modes_ref_pdb, mode=7, pdb=ref_pdb, file="NMA_reference_pdb_mode7_traj.pdb")
# Vector field representation
pymol(modes_ref_pdb, mode=7, pdb=ref_pdb, file="NMA_reference_pdb_mode7.pml", type="script")

# # Plot a heat map with clustering dendogram
# ## The similarity of structural dynamics is calculated by RMSIP based on the 10 lowest frequency normal modes.
# ## The RMSIP values are pre-calculated in the modes object and can be accessed through the attribute modes$rmsip
# png("pair-wise_RMSIPs.png", units="in", width=5, height=5, res=300)
# heatmap((1-modes$rmsip), labCol=ids, symm=TRUE, main="Pair-wise RMSIPs") #, labRow=annotation[, "state"]
# dev.off()


## Ensemble Difference Distance Matrix (eDDM) Analysis
##-------------------------------------
# The eDDM analysis compares structural ensembles under distinct ligation, activation, etc. conditions.
# At least two groups of structures are required.
# we use PCA of the distance matrices (through the bio3d function, pca.array()) followed by a conventional hierarchical
# clustering in the PC1-PC2 subspace to get the intrinsic grouping of the structures.
library(bio3d.eddm)
# need to update the aligned structures to include all heavy atoms:
pdbs.aa <- read.all(pdbs)
# Calculate distance matrices.
# The option `all.atom=TRUE` tells that all heavy-atom coordinates will be used.
dm <- dm(pdbs.aa, all.atom=TRUE)
# Perform PCA of distance matrices.
pc_dm <- pca.array(dm)
# # Plot importance of PCs
# png("PCA_on_all_atom_DistanceMatrix_PCimportance.png", units="in", width=5, height=5, res=300)
# plot.pca.scree(pc_dm)
# dev.off()
# Perform structural clustering in the PC1-PC2 subspace.
hc_dm_pc12 <- hclust(dist(pc_dm$z[, 1:2]))
grps_dm_pc12 <- cutree(hc_dm_pc12, k=number_of_groups)
# Plot PCs
png("PCA_on_all_atom_DistanceMatrix.png", units="in", width=5, height=5, res=300)
plot(pc_dm, col=grps_dm_pc12) # for only PC1-2 subplot add:, pc.axes=c(1,2)
dev.off()
# Calculating difference distance matrices between groups
# the difference mean distance between groups for each residue pair are calculated and statistical significance assessed
# using a two-sample Wilcoxon test. Long-distance pairs in all structures are omitted.
tbl <- eddm(pdbs.aa, grps=grps_dm_pc12, dm=dm, mask="cmap")
# Plot
png("eDDM_complete.png", units="in", width=5, height=5, res=300)
plot(tbl, pdbs=pdbs.aa, full=TRUE, resno=NULL, sse=pdbs$sse[1, ], type="tile") #,labels=TRUE, labels.ind=c(1:3, 17:19)
dev.off()
# Generate PyMol script to visualise all residue pairs showing any distance changes
pymol(tbl, pdbs=pdbs.aa, grps=grps_dm_pc12, as="sticks", file="eDDM_any.pml", type="script")
# Identifying significant distance changes
keys <- subset.eddm(tbl, alpha=0.005, beta=1.0, switch.only=TRUE)
if (length(keys)) {
    # Plot
    png("eDDM_significant.png", units="in", width=5, height=5, res=300)
    plot(keys, pdbs=pdbs.aa, full=TRUE, resno=NULL, sse=pdbs$sse[1, ], type="tile") #,labels=TRUE, labels.ind=c(1:3, 17:19)
    dev.off()
    # Generate PyMol script to visualise all identified key residue pairs showing significant distance changes
    pymol(keys, pdbs=pdbs.aa, grps=grps_dm_pc12, as="sticks", file="eDDM_significant.pml", type="script")
    } else {print("No key switch residues identified")}


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


# Protect spaces in path names with gsub(" ","\\\\ ",pathname)
scriptpath = gsub(" ","\\\\ ",paste(projectdir,'/src/analysis_bio3d_reporting.py', sep=''))
str_input_path = gsub(" ","\\\\ ",opt$indir)
output_path = gsub(" ","\\\\ ",outdir)

### Save R session info
writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info_", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
#writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info.txt"))

### Generate Report
system(paste('python3', scriptpath, str_input_path, output_path, number_of_groups, sep=' '), wait=FALSE)
