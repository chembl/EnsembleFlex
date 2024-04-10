#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(bio3d)
library(msa)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="Bio3D_NMA",
              help="output directory [default=%default]", metavar="character"),
  make_option(c("-n", "--ngroups"), type="integer", default=3,
              help="number of groups for clustering [default=%default]", metavar="integer"),
  make_option(c("-e", "--eNMA"), type="logical", action="store_true", default=FALSE,
              help="if TRUE eNMA is computed on the whole ensemble [default=%default]", metavar="TRUE or FALSE")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check provided output directory <opt$outdir>

if (opt$outdir == "Bio3D_NMA"){
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


## The actual program...
#--------------------------------------------------------------------

# loading pdb files
pdbs <- pdbaln(files, exefile='msa')

ids <- sub("[.].*", "", basename(pdbs$id)) # get filenames and drop any extensions

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

ca.ref_pdb <- trim.pdb(ref_pdb, "calpha")


## Anisotropic Network Model (ANM) based NMA for single reference structure (default)
##-------------------------------------
modes_ref_pdb <- nma(ref_pdb)
png(filename="ANM_NMA_reference_pdb.png", width=900, height=750, units="px", res=120)
plot.nma(modes_ref_pdb, resno=ref_pdb, sse=ref_pdb, sse.min.length=3)#, main="NMA on reference structure"
dev.off()
print("Plot saved to file ANM_NMA_reference_pdb.png")

# Make a PDB trajectory
mktrj(modes_ref_pdb, mode=7, pdb=ref_pdb, file="ANM_NMA_reference_pdb_mode7_traj.pdb")
print("Interpolated trajectory structures saved to file ANM_NMA_reference_pdb_mode7_traj.pdb")
# Vector field representation
pymol(modes_ref_pdb, mode=7, pdb=ref_pdb, file="ANM_NMA_reference_pdb_mode7.pml", type="script")

# Dynamic Cross-Correlation from ANM
cm_anm <- dccm.nma(modes_ref_pdb)
# Plot correlation map
png(filename="ANM_NMA_dynamic_cross_correlations_reference_pdb.png", width=900, height=750, units="px", res=120)
plot(cm_anm, resno=ref_pdb, sse=ref_pdb, contour = FALSE, col.regions = bwr.colors(20),
     at = seq(-1, 1, 0.1))
dev.off()
# # DCCM PyMOL visualization: save a PDB file with CONECT records (when argument type='pdb')
# pymol.dccm(cm_anm, pdb=ref_pdb, step=0.2, omit=0.2, radius = 0.15, type="pdb",
#         file="ANM_NMA_dynamic_cross_correlations.pdb")

# Save reference structure with ANM_fluctuations in B-factor column
tryCatch(
    #try to ...
    {
    write.pdb(ca.ref_pdb, b=modes_ref_pdb$fluctuations, file="ANM_fluctuations_onReference.pdb")
    print("PDB saved to file ANM_fluctuations_onReference.pdb")
    },
    #if an error occurs, tell me the error
    error=function(e) {
        message('An Error Occurred. This is probably related to the reference structure.')
        print(e)
    }
)


## Gaussian Network Model (GNM) based NMA for single reference structure
##-------------------------------------
modes_gnm_ref_pdb <- gnm(ref_pdb)
png(filename="GNM_NMA_reference_pdb.png", width=900, height=750, units="px", res=120)
plot.nma(modes_gnm_ref_pdb, resno=ref_pdb, sse=ref_pdb, sse.min.length=3)#, main="NMA on reference structure"
dev.off()
print("Plot saved to file GNM_NMA_reference_pdb.png")

# Dynamic Cross-Correlation from Gaussian Network Model
cm_gnm <- dccm.gnm(modes_gnm_ref_pdb)
# Plot correlation map
png(filename="GNM_NMA_dynamic_cross_correlations_reference_pdb.png", width=900, height=750, units="px", res=120)
plot(cm_gnm, resno=ref_pdb, sse=ref_pdb, contour = FALSE, col.regions = bwr.colors(20),
     at = seq(-1, 1, 0.1))
dev.off()
# # DCCM PyMOL visualization: save a PDB file with CONECT records (when argument type='pdb')
# pymol(cm_gnm, ref_pdb, step=0.2, omit=0.2, radius = 0.15, type="pdb",
#         file="GNM_NMA_dynamic_cross_correlations.pdb")

# Save reference structure with GNM_fluctuations in B-factor column
tryCatch(
    #try to ...
    {
    write.pdb(ca.ref_pdb, b=modes_gnm_ref_pdb$fluctuations, file="GNM_fluctuations_onReference.pdb")
    print("PDB saved to file GNM_fluctuations_onReference.pdb")
    },
    #if an error occurs, tell me the error
    error=function(e) {
        message('An Error Occurred. This is probably related to the reference structure.')
        print(e)
    }
)


## Comparison of PCs and NMA modes of reference structure
##-------------------------------------
# PCA on coordinates (backbone)
pc_xyz <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])

## Root Mean Square Inner Product (RMSIP)
## RMSIP is a measure for the similarity between two set of modes obtained from principal component or normal modes analysis.
# Calculate the RMSIP between the PCs and the NMA-modes
r <- rmsip(modes_ref_pdb, pc_xyz, subset=10, row.name="NMA", col.name="PCA")
# Plot pairwise overlap values
png("RMSIP_ensemble_PC_NMA_reference_pdb.png", units="in", width=5, height=5, res=300)
plot(r, xlab="NMA", ylab="PCA", main="RMSIP between ensemble PCs and \nNMA-modes of the reference structure")
dev.off()


# ## eNMA
# ##-------------------------------------
if (opt$eNMA==TRUE){
    library("pheatmap")
    ## NMA on all structures;
    ## use 'rm.gaps=FALSE' to keep the gap containing columns, but note that this is not compatible with mktrj
    modes <- nma.pdbs(pdbs, fit=FALSE, rm.gaps=TRUE)

    ## Set number_of_groups for clustering
    number_of_groups <- opt$ngroups

    ## Cluster on Fluctuation similarity
    sip <- sip(modes)
    hc <- hclust(dist(sip))
    col <- cutree(hc, k=number_of_groups)

    ## Plot fluctuation data
    png("eNMA_fluctuations.png", units="in", width=5, height=5, res=300)
    plot(modes, pdbs=pdbs, col=col, main="Normal Mode fluctuations") #, col=annotation[, "color"]
    dev.off()

    ## Interpolated structures along eNMA modes produced by the mktrj() function
    mktrj(enma = modes, pdbs = pdbs, mag = 10, step = 1.25, file = "eNMA.pdb", rock = TRUE)
    # Vector field representation
    #pymol(modes, pdb=ref_pdb, mode=7, file="eNMA_col_mode7.pml", type="script")

#     # Dynamic Cross-Correlation from ANM based eNMA
#     cm <- dccm.enma(modes)
#     # Plot correlation map
#     png(filename="eNMA_dynamic_cross_correlations_reference_pdb.png", width=900, height=750, units="px", res=120)
#     plot(cm, resno=ref_pdb, sse=ref_pdb, contour = FALSE, col.regions = bwr.colors(20),
#          at = seq(-1, 1, 0.1))
#     dev.off()
#     # DCCM PyMOL visualization: save a PDB file with CONECT records (when argument type='pdb')
#     pymol(cm, ref_pdb, step=0.2, omit=0.2, radius = 0.15, type="pdb",
#             file="eNMA_dynamic_cross_correlations.pdb")

    ## eNMA RMSIP with clustering dendrogram
    ## The similarity of structural dynamics is calculated by RMSIP based on the 10 lowest frequency normal modes.
    ## The RMSIP values are pre-calculated in the modes object and can be accessed through the attribute modes$rmsip
    # Plot a heatmap with clustering dendrogram
    png("pair-wise_RMSIPs.png", units="in", width=5, height=5, res=300)
    #heatmap((1-modes$rmsip), labCol=ids, symm=TRUE, main="Pair-wise RMSIPs")
    pheatmap(modes$rmsip, labels_col=ids, symm=TRUE, main="Pair-wise RMSIPs")
    dev.off()

    ## Bhattacharyya coefficient
    ## The Bhattacharyya coefficient is a measure of the amount of overlap between two statistical samples or populations.
    bc <- bhattacharyya(modes)
    png("pair-wise_bhattacharyya_coeffs.png", units="in", width=5, height=5, res=300)
    #heatmap((1-bc), labCol=ids, symm=TRUE, main="Pair-wise Bhattacharyya coefficients")
    pheatmap(bc, labels_col=ids, symm=TRUE, main="Pair-wise Bhattacharyya coefficients")
dev.off()
}