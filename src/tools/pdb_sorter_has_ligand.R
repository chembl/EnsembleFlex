#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(bio3d)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir",
              help="output directory name [default=%default]", metavar="character")
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

if (!is.null(opt$indir)){
  setwd(outdir)
  # list of pdb files in specified folder (corefit structures containing one ligand)
  file_list <- list.files(path = opt$indir, pattern = "*.pdb", full.names = T, recursive = F)
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}


## program...

## Ligand identification and pdb sorting
#---------

# create sub-directories
dir.create(file.path(outdir, "structures_with_ligand"))
dir.create(file.path(outdir, "structures_without_ligand"))

# create empty table to save binding site data for all structures
pdb_table  <- data.frame(Filename = character(), Has_ligand = numeric(), Ligand_IDs = character(), Full_path = character(), stringsAsFactors = FALSE)

for(i in 1:length(file_list)){
  filename = file_list[[i]]
  #  pdbID <- tools::file_path_sans_ext(basename(filename))
  pdbID <- strsplit((basename(filename)), '[.]')[[1]][1]
  # read the pdb file
  one_pdb <- read.pdb(filename)

  # identify ligands
  # be aware that everything that is not protein, nucleic acid or water will be considered as ligand.
  lig.inds <- atom.select(one_pdb, "ligand")
  # Access the PDB data with the selected atom indices
  ligIDs <- one_pdb$atom[ lig.inds$atom, "resid" ]
  #print(basename(filename))
  #print(toString(unique(ligIDs)))
  # save number of identified ligands in file
  pdb_table[i,] <- list(basename(filename), length(unique(ligIDs)), paste(unique(ligIDs),collapse=' '), filename)

  # safe pdb structure in respective sub-directory
  if (length(ligIDs)>0) {
    write.pdb(one_pdb, file=paste(outdir,"/structures_with_ligand/",basename(filename),sep=''))
    print(paste(basename(filename),"has at least one ligand:",toString(unique(ligIDs)),sep=' '))
  } else {
    write.pdb(one_pdb, file=paste(outdir,"/structures_without_ligand/",basename(filename),sep=''))
    print(paste(basename(filename),"has no ligand",sep=' '))
  }
}

# safe dataframe
write.csv(pdb_table, "pdbs_have_ligands.csv", row.names=FALSE, quote=FALSE)

