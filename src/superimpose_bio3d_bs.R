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
# read the reference pdb file (first file in filelist)
pdb <- read.pdb(files[[1]])


# Read in the binding_site_residue_numbers file
# binding_residue_num <- read.table(file = opt$bsresidues, header = F, stringsAsFactors = F)
# binding_residue_num <- as.vector(binding_residue_num)
binding_residue_num <- scan(opt$bsresidues, sep="\n")
print(binding_residue_num)

# get indices of all bindingsite residue backbone atoms
binding.inds <- print(atom.select(pdb, "backbone", resno=c(binding_residue_num)))
#print(binding.inds)

# pdbfit is a wrapper for the function fit.xyz
# The reference frame for supperposition (i.e. the fixed structure to which others are superposed) is the first entry
# in the input "pdbs" object (fixed=pdbs$xyz[1,]). For finer control use fit.xyz.
xyz <- pdbfit.pdbs(pdbs, inds = binding.inds, outpath = "superimposed_on_bs", pdbext = "")


superimp_dir <- paste(outdir, "/superimposed_on_bs", sep='')
# remove double file extensions "pdb_flsq." from file names
superimp_files <- list.files(path = superimp_dir, pattern = "*.pdb_flsq.pdb", full.names = T, recursive = F)
sapply(superimp_files,FUN=function(eachPath){
  file.rename(from=eachPath,to= sub(pattern="pdb_flsq.", paste0(""),eachPath))
})
