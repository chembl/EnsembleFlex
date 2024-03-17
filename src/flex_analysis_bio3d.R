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
library("R.utils") # for function "isAbsolutePath"
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
} else if (length(strsplit(opt$outdir, split='/', fixed=TRUE)) == 1){
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
  files <- list.files(path = indir, pattern = "*.pdb", full.names = T, recursive = F)
  # } else if (!is.null(opt$filenames)) {
  #   files <- as.list(strsplit(opt$filenames, ",")[[1]])
  #   #files <- opt$filenames
  #   setwd(paste(dirname(files[1]),opt$outdir))
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}

print(indir)
setwd(outdir)
print(outdir)



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
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
dev.off()
print("Plot saved to file RMSD_hist.png")

## RMSD heatmap
library(pheatmap)
png("RMSD_heatmap.png", units="in", width=5, height=5, res=300)
#heatmap(rd, labCol=ids, main="RMSD Heatmap")
pheatmap(rd, main="RMSD Heatmap", fontsize = 6, show_colnames = FALSE) #annotation_row = ids
dev.off()
print("Plot saved to file RMSD_heatmap.png")

## RMSD Hierarchical clustering
hc_rmsd <- hclust(as.dist(rd))
png("RMSD_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_rmsd, labels=ids, cex=0.5, k=number_of_groups,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file RMSD_dendrogram.png")


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
           xlab="Residue No.", col="gray", main="RMSFs per Residue") #, typ="l"
dev.off()
print("Plot saved to file RMSF.png")

rf_withgaps <- rmsf(pdbs$xyz)
png("RMSF_including_gaps.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf_withgaps, rm.gaps=FALSE, resno=pdb1, sse=pdb1, sse.min.length=0, ylab="RMSF (Å)",
           xlab="Residue No.", col="gray", main="RMSFs per Residue") #, typ="l"
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



## B-factors
##-------------------------------------
png("B-factors.png", units="in", width=5, height=5, res=300)
plot.bio3d(pdbs$b, rm.gaps=TRUE, ylab="B-factor", #, resno=ref_pdb, sse=ref_pdb
           xlab="Residue No.", main="B-factors", col="gray") #, typ="l"
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

png("PCA_dendrogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_dendrogram.png")


# Save PyMOL scripts for coloring by group
grps_rmsd <- cutree(hc_rmsd, k=number_of_groups)
# color by clustering (based on RMSD)
pymol(pdbs, col=grps_rmsd, as="cartoon", file="col_by_grps_RMSD.pml", type="script")
# color by clustering (based on PCA)
pymol(pdbs, col=grps_pc12, as="cartoon", file="col_by_grps_PC.pml", type="script")


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

png("PCA_on_Torsion_dendogram.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12_tor, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="PC Cluster Dendrogram (on Torsion)", fillbox=FALSE) #, colors=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_on_Torsion_dendogram.png")

# ## eNMA
# ##-------------------------------------
# # NMA on all structures; use 'rm.gaps=FALSE' to keep the gap containing columns, but note that this is not compatible
# # with mktrj
# modes <- nma.pdbs(pdbs, rm.gaps=TRUE)
#
# # Make fluctuation plot
# png("eNMA_fluctuations.png", units="in", width=5, height=5, res=300)
# plot(modes, pdbs=pdbs, main="Normal Mode fluctuations") #, col=annotation[, "color"]
# dev.off()
#
# ## Interpolated structures along eNMA modes produced by the mktrj() function
# mktrj(enma = modes, pdbs = pdbs, mag = 10, step = 1.25, file = "eNMA.pdb", rock = TRUE)
# # Vector field representation
# #pymol(modes, pdb=ref_pdb, mode=7, file="eNMA_col_mode7.pml", type="script")


# NMA for single reference structure
##-------------------------------------
modes_ref_pdb <- nma(ref_pdb)
png(filename="NMA_fluctuations_reference_pdb.png", width=900, height=750, units="px", res=120)
plot.nma(modes_ref_pdb, resno=ref_pdb, sse=ref_pdb, sse.min.length=3)#, main="NMA on reference structure"
dev.off()
print("Plot saved to file NMA_fluctuations_reference_pdb.png")
# Make a PDB trajectory
mktrj(modes_ref_pdb, mode=7, pdb=ref_pdb, file="NMA_reference_pdb_mode7_traj.pdb")
print("Interpolated trajectory structures saved to file NMA_reference_pdb_mode7_traj.pdb")
# Vector field representation
pymol(modes_ref_pdb, mode=7, pdb=ref_pdb, file="NMA_reference_pdb_mode7.pml", type="script")

# # Plot a heat map with clustering dendogram
# ## The similarity of structural dynamics is calculated by RMSIP based on the 10 lowest frequency normal modes.
# ## The RMSIP values are pre-calculated in the modes object and can be accessed through the attribute modes$rmsip
# png("pair-wise_RMSIPs.png", units="in", width=5, height=5, res=300)
# heatmap((1-modes$rmsip), labCol=ids, symm=TRUE, main="Pair-wise RMSIPs") #, labRow=annotation[, "state"]
# dev.off()


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
hist(rd_allatom, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
dev.off()
print("Plot saved to file RMSD_hist_allatom.png")

## all-atom RMSD heatmap
library(pheatmap)
png("RMSD_heatmap_allatom.png", units="in", width=5, height=5, res=300)
#heatmap(rd_allatom, labCol=ids, main="RMSD Heatmap")
pheatmap(rd_allatom, main="RMSD Heatmap", fontsize = 6, show_colnames = FALSE) #annotation_row = ids
dev.off()
print("Plot saved to file RMSD_heatmap_allatom.png")

## all-atom RMSD Hierarchical clustering
hc_rmsd_allatom <- hclust(as.dist(rd_allatom))
png("RMSD_dendrogram_allatom.png", units="in", width=5, height=5, res=300)
hclustplot(hc_rmsd_allatom, labels=ids, cex=0.5, k=number_of_groups,
           ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
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

png("PCA_allatom.png", units="in", width=5, height=5, res=300)
plot(pc_xyz_allatom, col=grps_pc12_allatom)
dev.off()
print("Plot saved to file PCA_allatom.png")

png("PCA_dendrogram_allatom.png", units="in", width=5, height=5, res=300)
hclustplot(hc_pc12_allatom, labels=ids, cex=0.5, k=number_of_groups,
           ylab="PC1-2 distance", main="PC Cluster Dendrogram", fillbox=FALSE) #, colors=annotation[, "color"]
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



## Ensemble Difference Distance Matrix (eDDM) Analysis
##-------------------------------------
# The eDDM analysis compares structural ensembles under distinct ligation, activation, etc. conditions.
# At least two groups of structures are required.
# we use PCA of the distance matrices (through the bio3d function, pca.array()) followed by a conventional hierarchical
# clustering in the PC1-PC2 subspace to get the intrinsic grouping of the structures.
library(bio3d.eddm)
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
    pymol(keys, pdbs=pdbs_allatoms, grps=grps_dm_pc12, as="sticks", file="eDDM_significant.pml", type="script")
    } else {print("No key switch residues identified.")}


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
## Cluster attributions for RMSD, PC, RMSD-allatom, PC-allatom
##-------------------------------------
# store cluster attributions in dataframe
clusters_df <- as.data.frame(list(ids, grps_rmsd, grps_pc12, grps_pc12_tor, grps_rmsd_allatom, grps_pc12_allatom),
                             col.names = c("Structure", "RMSD Cluster", "PC1/2 Cluster onCoords", "PC1/2 Cluster onTorsion", "RMSD-allatom Cluster", "PC1/2-allatom Cluster onCoords"))
# safe dataframe
write.csv(clusters_df, "cluster_attributions.csv", row.names=FALSE, quote=FALSE)


##-------------------------------------
# Protect spaces in path names with gsub(" ","\\\\ ",pathname)
scriptpath = gsub(" ","\\\\ ",paste(projectdir,'/src/analysis_bio3d_reporting.py', sep=''))
str_input_path = gsub(" ","\\\\ ",opt$indir)
output_path = gsub(" ","\\\\ ",outdir)

### Save R session info
##-------------------------------------
writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info_", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
#writeLines(capture.output(sessionInfo()), paste0(output_path, "/R_session_info.txt"))
print(paste0("R session info saved to file R_session_info_", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))

### Generate Report
##-------------------------------------
system(paste('python3', scriptpath, str_input_path, output_path, number_of_groups, sep=' '), wait=FALSE)
