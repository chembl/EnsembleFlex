#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("optparse")
library("R.utils") # for function "isAbsolutePath"
library("bio3d")
library("msa")

option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir",
              help="output directory [default=%default]", metavar="character"),
  make_option(c("-b", "--bsresidues"), type="character", default=NULL,
              help="binding site residues as file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check provided output directory <opt$outdir>

if (opt$outdir == "outdir"){
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
  print(indir)
  files <- list.files(path = indir, pattern = "*.pdb", full.names = T, recursive = F)
  # } else if (!is.null(opt$filenames)) {
  #   files <- as.list(strsplit(opt$filenames, ",")[[1]])
  #   #files <- opt$filenames
  #   setwd(paste(dirname(files[1]),opt$outdir))
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}

## Check provided binding site residue file

if (file.exists(file.path(opt$bsresidues))){ #if (!is.null(opt$indir)){
  if (isAbsolutePath(opt$bsresidues)){
    # full path to existing directory is provided
    bsitefile <- file.path(opt$bsresidues)
  } else {
    # only subdirectory is provided, but full path exists
    bsitefile <- file.path(getwd(), opt$bsresidues)
  }
} else {
    print_help(opt_parser)
    stop("Binding site residue file must be provided.\n\n", call.=FALSE)
}

print(indir)
setwd(outdir)
print(outdir)


## program...

# loading pdb files
#pdbs <- pdbaln(files)
pdbs <- pdbaln(files, exefile='msa')
# read all atoms for each residue
pdbs_allatoms <- read.all(pdbs)
# read the reference pdb file (first file in filelist)
pdb <- read.pdb(files[[1]])


# Read in the binding_site_residue_numbers file
binding_residue_num <- scan(bsitefile, sep="\n")
print(binding_residue_num)

# get indices of all bindingsite residue backbone atoms
#binding.inds <- print(atom.select(pdb, "backbone", resno=c(binding_residue_num)))
#print(binding.inds)

# Tailor a reference pdb structure to only binding site residues
bsite.pdb = trim.pdb(pdb, inds=atom.select(pdb, resno=c(binding_residue_num)))

#-------------------------------
# RMSF of binding site residues including side chain contribution
#-------------------------------
# make shure that the residue is present in all structures

rmsf_bsite <- list()
for(i in 1:length(binding_residue_num)){
  resno = binding_residue_num[i]
  # Get indices for all atoms for resno
  bsiteres.pos = atom.select(bsite.pdb, resno=resno, string="protein")
  # residue numbers to group by
  resno_gr1 <- bsite.pdb$atom$resno[bsiteres.pos$atom]
  # mean rmsf value of all atoms of each residue
  rmsf_bsite[i] <- rmsf(pdbs_allatoms$all[, bsiteres.pos$xyz], grpby=resno_gr1)
}
rmsf_bsite <- unlist(rmsf_bsite)

#---------
# Same as above for all residues present in selection at once (instead of separate by residue)
# Disadvantage: rmsf results are ordered by increasing residue number
# Get indices for all atoms
bsite.pos = atom.select(bsite.pdb, string="protein")
# residue numbers to group by
resno_gr <- bsite.pdb$atom$resno[bsite.pos$atom]
# mean rmsf value of all atoms of each residue
rf_bsite <- rmsf(pdbs_allatoms$all[, bsite.pos$xyz], grpby=resno_gr)
# PCA on all atoms of binding site residues
pc_xyz_bsite_allatom <- pca(pdbs_allatoms$all[, bsite.pos$xyz], rm.gaps=TRUE, fit=F)
pc_xyz_bsite_allatom

#---------

# # check:
# unique(resno_gr) %>% length()
# sort(binding_residue_num) %>% length()
# rf_bsite %>% length()
# rmsf_bsite %>% length()

#---------
# Plotting
#---------
# Plot RMSF labeled with residue numbers
png(filename="RMSF_bsite.png", width=900, height=750, units="px", res=120)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,5,3,1)) # adjust margins.
barplot(rmsf_bsite, names.arg=binding_residue_num,
        main="All-atom mean RMSF of residues implied in ligand binding",
        cex.names=0.8, ylab="All-atom mean RMSF per residue [Å]") #, ylim=c(0, 4.5)
dev.off()

png(filename="RMSF_bsite2.png", width=1200, height=750, units="px", res=120)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,5,3,1)) # adjust margins.
barplot(rf_bsite, names.arg=unique(resno_gr),
        main="All-atom mean RMSF of residues implied in ligand binding",
        cex.names=0.8, ylab="All-atom mean RMSF per residue [Å]") #, ylim=c(0, 4.5)
dev.off()

# Save reference structure with RMSF in B-factor column
ca.bsite.pdb <- trim.pdb(bsite.pdb, "calpha")
write.pdb(ca.bsite.pdb, b=rf_bsite, file="RMSF_onBsite_allatom.pdb")
print("PDB saved to file RMSF_onBsite_allatom.pdb")


# Plot PCs
png("PCA_bsite_allatom.png", units="in", width=5, height=5, res=300)
plot(pc_xyz_bsite_allatom) #, col=annotation[, "color"]
dev.off()
print("Plot saved to file PCA_bsite_allatom.png")

## Atom contribution to all-atom bsite PCA
png("PCA_atom_contribution_bsite_allatom.png", units="in", width=5, height=5, res=300)
par(mfrow = c(3, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc_xyz_bsite_allatom$au[,1], resno=pdb, sse=pdb, ylab="PC1")
plot.bio3d(pc_xyz_bsite_allatom$au[,2], resno=pdb, sse=pdb, ylab="PC2")
plot.bio3d(pc_xyz_bsite_allatom$au[,3], resno=pdb, sse=pdb, ylab="PC3")
dev.off()

## Interpolated structures along PC1/2/3 produced by the mktrj.pca() function
mktrj.pca(pc_xyz_bsite_allatom, pc=1, file="PC1_bsite_allatom.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz_bsite_allatom, pc=2, file="PC2_bsite_allatom.pdb") #, mag = 1, step = 0.125
mktrj.pca(pc_xyz_bsite_allatom, pc=3, file="PC3_bsite_allatom.pdb") #, mag = 1, step = 0.125

# PC_bsite-allatom Vector field representation
pymol(pc_xyz_bsite_allatom, pdb=pdbs_allatoms[1]$all, pc=1, as="wire", file="PC1vectors_bsite_allatom.pml", type="script")
pymol(pc_xyz_bsite_allatom, pdb=pdbs_allatoms[1]$all, pc=2, as="wire", file="PC2vectors_bsite_allatom.pml", type="script")
pymol(pc_xyz_bsite_allatom, pdb=pdbs_allatoms[1]$all, pc=3, as="wire", file="PC3vectors_bsite_allatom.pml", type="script")

