#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install.packages("bio3d", dependencies=TRUE)
#install.packages("optparse", dependencies=TRUE)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")
#install.packages("pheatmap")

library("optparse")
library("bio3d")
library("msa")

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="dataset directory path", metavar="character"),
  # make_option(c("-f", "--filenames"), type="character", default=NULL, 
  #             help="dataset file names - example: [file1.pdb,file2.pdb]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Bio3D_Analysis",
              help="output directory [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$directory)
print(opt$out)


if (opt$out == "Bio3D_Analysis"){
  outdir <- file.path(opt$directory, opt$out)
  dir.create(outdir, showWarnings = FALSE)
} else {
  outdir <- file.path(opt$out)
  dir.create(outdir, showWarnings = FALSE)
}
setwd(outdir)
print(outdir)

if (!is.null(opt$directory)){
  setwd(outdir)
  files <- list.files(path = opt$directory, pattern = "*.pdb", full.names = T, recursive = F)
  # } else if (!is.null(opt$filenames)) {
  #   files <- as.list(strsplit(opt$filenames, ",")[[1]])
  #   #files <- opt$filenames
  #   setwd(paste(dirname(files[1]),opt$out))
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
number_of_groups <- 2
#ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))
ids <- unlist(substr(basename(pdbs$id), 1, 7))
print(ids)

# # option 2. if they are of identical composition and you only want xyz coordinates
# xyz <- NULL
# for(i in files) {
#   xyz = rbind(xyz, read.pdb(files)$xyz)
# }

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
pheatmap(rd, labCol=ids, main="RMSD Heatmap")
dev.off()

## RMSD clustering
hc.rd <- hclust(as.dist(rd))
png("RMSD_clust.png", units="in", width=5, height=5, res=300)
hclustplot(hc.rd, labels=ids, cex=0.5, k=number_of_groups,
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

ref.pdb <- trim.pdb(pdb1, inds = atom.select(pdb1, resno = pdbs$resno[pdb1$id, gaps.res$f.inds]))

## Plot RMSF with SSE annotation and labeled with residue numbers
rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
png("RMSF.png", units="in", width=5, height=5, res=300)
plot.bio3d(rf, resno=ref.pdb, sse=ref.pdb, ylab="RMSF (Å)",
           xlab="Residue No.", typ="l", main="RMSFs")
dev.off()


## Torsion/Dihedral analysis
##-------------------------------------
# Locate the two structures in pdbs
#print(files[1], files[length(files)])
ind.a <- grep(pdb1$id, pdbs$id)
ind.b <- grep(pdb2$id, pdbs$id)

# Exclude gaps in the two structures to make them comparable
gaps.xyz2 <- gap.inspect(pdbs$xyz[c(ind.a, ind.b), ])
a.xyz <- pdbs$xyz[ind.a, gaps.xyz2$f.inds]
b.xyz <- pdbs$xyz[ind.b, gaps.xyz2$f.inds]

# Compare CA based pseudo-torsion angles between the two structures
a <- torsion.xyz(a.xyz, atm.inc=1)
b <- torsion.xyz(b.xyz, atm.inc=1)
d.ab <- wrap.tor(a-b)
d.ab[is.na(d.ab)] <- 0

# Plot results with SSE annotation
png("Ca_torsion_diff_between_first&last_structure.png", units="in", width=5, height=5, res=300)
plot.bio3d(abs(d.ab), resno=pdb1, sse=pdb1, typ="h", xlab="Residue No.",
           ylab = "Difference Angle", main="Ca torsion difference - first&last structure")
dev.off()


## Difference distance matrix analysis (DDM)
##-------------------------------------
a <- dm.xyz(a.xyz)
b <- dm.xyz(b.xyz)

png("DDM.png", units="in", width=5, height=5, res=300)
plot.dmat( (a - b), nlevels=10, grid.col="gray", xlab=basename(pdb1$id), ylab=basename(pdb2$id), main="Difference distance matrix")
dev.off()


## PCA
##-------------------------------------
# Do PCA on coordinates
pc.xray <- pca.xyz(pdbs$xyz[, gaps.pos$f.inds])
pc.xray

png("PCA.png", units="in", width=5, height=5, res=300)
plot(pc.xray) #, col=annotation[, "color"]
dev.off()

## Residue contribution to PCA
png("PCA_residue_contribution.png", units="in", width=5, height=5, res=300)
par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno=ref.pdb, sse=ref.pdb, ylab="PC1")
plot.bio3d(pc.xray$au[,2], resno=ref.pdb, sse=ref.pdb, ylab="PC2")
plot.bio3d(pc.xray$au[,3], resno=ref.pdb, sse=ref.pdb, ylab="PC3")
dev.off()

## Interpolated structures along PC1 produced by the mktrj.pca() function
mktrj.pca(pc.xray, pc=1, file="pc1.pdb")


## NMA
##-------------------------------------
# NMA on all structures; use 'rm.gaps=FALSE' to keep the gap containing columns
modes <- nma.pdbs(pdbs, rm.gaps=FALSE)

# Make fluctuation plot
png("NMA_fluctuations.png", units="in", width=5, height=5, res=300)
plot(modes, pdbs=pdbs, main="Normal Mode fluctuations") #, col=annotation[, "color"]
dev.off()

# # Plot a heat map with clustering dendogram
# ## The similarity of structural dynamics is calculated by RMSIP based on the 10 lowest frequency normal modes.
# ## The RMSIP values are pre-calculated in the modes object and can be accessed through the attribute modes$rmsip
# png("pair-wise_RMSIPs.png", units="in", width=5, height=5, res=300)
# heatmap((1-modes$rmsip), labCol=ids, symm=TRUE, main="Pair-wise RMSIPs") #, labRow=annotation[, "state"]
# dev.off()


# # Calculate pair-wise RMSD values
# rmsd.map <- rmsd(pdbs$xyz, a.inds=gaps.pos$f.inds, fit=TRUE)
# png("pair-wise_RMSDs.png", units="in", width=5, height=5, res=300)
# heatmap(rmsd.map, labCol=ids, symm=TRUE, main="Pair-wise RMSDs") #, labRow=annotation[, "state"]
# dev.off()

#writeLines(capture.output(sessionInfo()), paste0(outdir, "/Data/session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
writeLines(capture.output(sessionInfo()), paste0(outdir, "/R_session_info.txt"))
