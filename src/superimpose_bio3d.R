#!/usr/bin/env Rscript

# """
# Ensemble structural superpositioning using the R package Bio3D.
#
# Usage:
#     Rscript superimpose_bio3d.R -i <input_directory> -o <output_directory>
#
# Example:
#     Rscript superimpose_bio3d.R -i pdbs -o EnsemblFlex
# """

args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(bio3d)
library(msa)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL, 
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir",
              help="output directory [default=%default]", metavar="character")
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
outdir <- getwd() # to get absolute path in variable


## program...
#--------------------------------------------------------------------

# loading pdb files
pdbs <- pdbaln(files, exefile='msa')

# Rigid core identification and structural superposition
core <- core.find(pdbs)
# use core inds for structural superposition
core.inds <- print(core, vol=1.0)
write.pdb(xyz = pdbs$xyz[1,core.inds$xyz], file = "quick_core.pdb")

# pdbfit is a wrapper for the function fit.xyz
# The reference frame for supperposition (i.e. the fixed structure to which others are superposed) is the first entry
# in the input "pdbs" object (fixed=pdbs$xyz[1,]). For finer control use fit.xyz.
xyz <- pdbfit.pdbs(pdbs, inds = core.inds, outpath = "superimposed", pdbext = "")


superimp_dir <- paste(outdir, "/superimposed", sep='')
# remove double file extensions "pdb_flsq." from file names
superimp_files <- list.files(path = superimp_dir, pattern = "*.pdb_flsq.pdb", full.names = T, recursive = F)
sapply(superimp_files,FUN=function(eachPath){
  file.rename(from=eachPath,to= sub(pattern="pdb_flsq.", paste0(""),eachPath))
})
