#!/usr/bin/env Rscript

# """
# Store data from data table in b-factor column of PDB structure using mainly the R package Bio3D.
# The residue number has to be provided in the first column of the data table and the name of the data column to be used
# has to be provided as argument.
#
# Usage:
#     Rscript data_on_structure.R -i <input_directory> -o <output_directory> -d <datafile> -c <column_name>
#
# Example:
#     Rscript data_on_structure.R -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_SASA_Biopython
#             -d EnsemblFlex/Analysis_SASA_Biopython/SASA_global.csv -c sd
# """

args = commandArgs(trailingOnly=TRUE)

projectdir = getwd()

library("optparse")
library("R.utils") # for function "isAbsolutePath"
library("bio3d")

option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="data_on_structure",
              help="output directory [default=%default]", metavar="character"),
  make_option(c("-d", "--datafile"), type="character", default=NULL,
              help="data as file", metavar="character"),
  make_option(c("-c", "--column"), type="character", default=NULL,
              help="data column name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check provided output directory <opt$outdir>

if (opt$outdir == "data_on_structure"){
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
  files <- list.files(path = indir, pattern = "*.pdb", full.names = T, recursive = F)
} else {
    print_help(opt_parser)
    stop("At least one input argument must be supplied (input filepath or files).\n\n", call.=FALSE)
}

## Check provided data file

if (file.exists(file.path(opt$datafile))){ #if (!is.null(opt$indir)){
  if (isAbsolutePath(opt$datafile)){
    # full path to existing directory is provided
    datafile <- file.path(opt$datafile)
  } else {
    # only subdirectory is provided, but full path exists
    datafile <- file.path(getwd(), opt$datafile)
  }
} else {
    print_help(opt_parser)
    stop("Data file must be provided.\n\n", call.=FALSE)
}

## Check provided column name

if (!is.null(opt$column)){ #if (!is.null(opt$indir)){
} else {
    print_help(opt_parser)
    stop("Column name must be provided.\n\n", call.=FALSE)
}


setwd(outdir)



## The actual program...
#--------------------------------------------------------------------

# read in the data file
# data <- scan(datafile, sep="\n")
dataframe <- read.csv(datafile)
print(dataframe)

# read in the pdb structure file
pdb1 <- read.pdb(files[1])

# set all b-factors to 0
pdb1$atom$b <- 0

for(i in 1:length(dataframe[,1])){
    # get residue number from first column in dataframe
    residue_num <- dataframe[i,1]
#     # get data value from second column in dataframe
#     value <- dataframe[i,2]
#     # get residue number from column 'ResNum' in dataframe
#     residue_num <- dataframe$ResNum[i]
#     # get data value from column 'sd' in dataframe
#     value <- dataframe$sd[i]
    print(residue_num)
    # get data value from provided column in dataframe
    value <- dataframe[i,opt$column]
    print(value)
    # get indices of residue
    residue_inds <- atom.select(pdb1, resno=c(residue_num))
    print(residue_inds)
    # store data value in b-factor column for each residue in PDB file
    print(pdb1$atom$b[ residue_inds$atom ])
    pdb1$atom$b[ residue_inds$atom ] <- value
}

# write to file
write.pdb(pdb1, file=paste0(opt$column,"_data_on_structure.pdb"))