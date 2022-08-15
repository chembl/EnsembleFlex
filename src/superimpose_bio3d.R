#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#install.packages("bio3d", dependencies=TRUE)
#install.packages("optparse", dependencies=TRUE)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")

library("optparse")
library("bio3d")
library("msa")

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="dataset directory path", metavar="character"),
  # make_option(c("-f", "--filenames"), type="character", default=NULL, 
  #             help="dataset file names - example: [file1.pdb,file2.pdb]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="outdir", 
              help="output directory [default=%default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$directory)
print(opt$out)


if (opt$out == "outdir"){
  outdir <- file.path(opt$directory, opt$out)
  dir.create(outdir, showWarnings = FALSE)
} else {
  outdir <- file.path(opt$out)
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


## program...

# loading pdb files
# option 1.
#pdbs <- pdbaln(files)
pdbs <- pdbaln(files, exefile='msa')

# # option 2. if they are of identical composition and you only want xyz coordinates
# xyz <- NULL
# for(i in files) {
#   xyz = rbind(xyz, read.pdb(files)$xyz)
# }

# Rigid core identification and structural superposition
core <- core.find(pdbs)
# use core inds for structural superposition
core.inds <- print(core, vol=1.0)
write.pdb(xyz = pdbs$xyz[1,core.inds$xyz], file = "quick_core.pdb")

xyz <- pdbfit.pdbs(pdbs, inds = core.inds, outpath = "superimposed")



