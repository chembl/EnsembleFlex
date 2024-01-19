#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(bio3d)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir",
              help="output directory name [default=%default]", metavar="character")
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

print(indir)
setwd(outdir)
print(outdir)


## program...

## Ligand identification and pdb sorting
#---------

# create sub-directories
dir.create(file.path(outdir, "structures_with_ligand"))
dir.create(file.path(outdir, "structures_without_ligand"))

# create empty table to save binding site data for all structures
pdb_table  <- data.frame(Filename = character(), Has_ligand = numeric(), Ligand_IDs = character(), Full_path = character(), stringsAsFactors = FALSE)

for(i in 1:length(files)){
  filename = files[[i]]
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

