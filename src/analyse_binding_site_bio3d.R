#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("optparse")
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


if (opt$outdir == "outdir"){
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
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}


## program...

# loading pdb files
#pdbs <- pdbaln(files)
pdbs <- pdbaln(files, exefile='msa')
# read all atoms for each residue
pdbs_allatoms <- read.all(pdbs)
# read the reference pdb file (first file in filelist)
pdb <- read.pdb(files[[1]])


# Read in the binding_site_residue_numbers file
binding_residue_num <- scan(opt$bsresidues, sep="\n")
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
#---------

# check:
unique(resno_gr) %>% length()
sort(binding_residue_num) %>% length()
rf_bsite %>% length()
rmsf_bsite %>% length()

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

