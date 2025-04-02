#!/usr/bin/env Rscript

# """
# Conserved water analysis using mainly the R packages vanddraabe (modified version) and Bio3D.
#
# Usage:
#     Rscript analyse_conserved_waters.R -i <input_directory> -o <output_directory>
#
# Example:
#     Rscript analyse_conserved_waters.R -i EnsemblFlex/superimposed -o EnsemblFlex/Analysis_Bio3D
# """

args = commandArgs(trailingOnly=TRUE)

projectdir = getwd()
projectdir = "~/ARISE/ARISE-Project/EnsembleFlex"

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(vanddraabe)
library(bio3d)
library(msa)
library(reshape2)
library(ggplot2)
library(cowplot)


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="Bio3D_Analysis",
              help="output directory [default=%default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## Check provided output directory <opt$outdir>

if (opt$outdir == "Water_Analysis"){
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

setwd(outdir)
abs_output_path = getAbsolutePath(".")

##-------------------------------------
# Protect spaces in path names with gsub(" ","\\\\ ",pathname)
str_input_path = gsub(" ","\\\\ ",indir)
output_path = gsub(" ","\\\\ ",abs_output_path)


## The actual program...
#--------------------------------------------------------------------

# loading pdb files
#pdbs <- pdbaln(files, exefile='msa')

# load reference pdb and get ID
pdb1 <- read.pdb(files[1])

conservedWaters <- ConservedWaters(prefix = indir,
                                   cluster = 2.4, #2.4
                                   mobility = 1.5, #2.0
                                   nBvalue = 1.0, #1.0
                                   prot.h2o.dist.min = 5.10, #5.10
                                   chain = "all",
                                   cluster.method = "complete",
                                   filename = "WaterAnalysis"
                                   )
#print(conservedWaters)

png(filename="ConservationPlot.png", width=8, height=5, units="in", res=150)
ConservationPlot(data = conservedWaters, passed.waters = TRUE)
dev.off()
print("Plot saved to file ConservationPlot.png")

png(filename="OccupancyBarplot.png", width=8, height=5, units="in", res=150)
OccupancyBarplot(data = conservedWaters, passed.waters = TRUE)
dev.off()
print("Plot saved to file OccupancyBarplot.png")

png(filename="MobilityBarplot.png", width=8, height=5, units="in", res=150)
MobilityBarplot(data = conservedWaters, passed.waters = TRUE)
dev.off()
print("Plot saved to file MobilityBarplot.png")

png(filename="BvalueBarplot.png", width=8, height=5, units="in", res=150)
BvalueBarplot(data = conservedWaters, passed.waters = TRUE, calc.values = FALSE)
dev.off()
print("Plot saved to file BvalueBarplot.png")

png(filename="BvalueBarplot_calculated_values.png", width=8, height=5, units="in", res=150)
BvalueBarplot(data = conservedWaters, passed.waters = TRUE, calc.values = TRUE)
dev.off()
print("Plot saved to file BvalueBarplot_calculated_values.png")

png(filename="nBvalueBarplot.png", width=8, height=5, units="in", res=150)
nBvalueBarplot(data = conservedWaters, passed.waters = TRUE)
dev.off()
print("Plot saved to file nBvalueBarplot.png")

png(filename="MobNormBvalEvalPlots.png", width=8, height=5, units="in", res=150)
MobNormBvalEvalPlots(data = conservedWaters, passed.waters = TRUE,
                      title = "Mobility and Normalized B-value Evaluation")
dev.off()
print("Plot saved to file MobNormBvalEvalPlots.png")

png(filename="BoundWaterEnvPlots.png", width=8, height=5, units="in", res=150)
BoundWaterEnvPlots(data = conservedWaters, passed.waters = TRUE,
                   pct.conserved.gte = 50, num.clusters = 50)
print("Plot saved to file MobNormBvalEvalPlots.png")

png(filename="BoundWaterEnvSummaryPlot.png", width=8, height=5, units="in", res=150)
BoundWaterEnvSummaryPlot(data = conservedWaters, passed.waters = TRUE,
                         title = "Bound Water Environment per Conservation Percentage")
print("Plot saved to file BoundWaterEnvSummaryPlot.png")


# Get ligand IDs
# be aware that everything that is not protein, nucleic acid or water will be considered as ligand.
lig.inds <- atom.select(pdb1, "ligand")
# Access the PDB data with the selected atom indices
ligID <- pdb1$atom[ lig.inds$atom, "resid" ]
# Ensure ligID is unique
ligID <- unique(ligID)
print(paste0("Detected ligand(s) are: ", ligID, "."))

# # If there are several ligands, take the first detected one
# Handle multiple ligands or a single ligand
if (length(ligID) > 1) {
  ligID <- as.character(ligID[1])  # Take the first unique ligand ID and convert to string
  print(paste0("Multiple ligands detected. Using the first ligand with ID ", ligID, " for further processing."))
} else {
  ligID <- as.character(ligID)  # Convert the single ligand ID to string
  print(paste0("Ligand with ID ", ligID, " used for pymol script generation."))
}

pdb1$id <- sub("[.].*", "", basename(files[1])) # get filename and drop any extensions
print(paste0("Reference structure used for pymol script generation: ", pdb1$id))

CreatePyMOLscript(conservedWaters.data = conservedWaters,
                   passed.waters = TRUE,
                   PDBid.ref = pdb1$id,
                   LigResname.ref = ligID,
                   hbond = 3.75,
                   lig.carbon.color = "cyan",
                   filename = "WaterAnalysis"
                   )


# Generate a custom PyMOL script for global display of conserved waters with the same color-coding as the plots

cat(
  "# Load the protein structure (pdb1) from the file specified in files[1]\n",
  sprintf("cmd.load('%s', 'pdb1')\n", files[1]),
  "\n",
  "# Delete all non-protein atoms from pdb1\n",
  "cmd.remove('pdb1 and not polymer')\n",
  "\n",
  "# Color the protein (polymer) chain in gray70\n",
  "cmd.color('gray70', 'pdb1')\n",
  "\n",
  "# Load the conserved waters file from the current directory\n",
  "cmd.load('WaterAnalysis_ConservedWaters_PASSED.pdb', 'ConservedWaters')\n",
  "\n",
  "# Define and display conserved waters\n",
  "cmd.alter('ConservedWaters', 'vdw=((q) / (8 * 3.14**2))**0.5 * 9')\n",
  "cmd.show(\"spheres\", \"ConservedWaters\")\n",
  "\n",
  "# Only display conserved waters with 40% or greater conservation\n",
  "cmd.show('spheres', 'ConservedWaters and (q>0.399999) and not hydro')\n",
  "\n",
  "# Create the custom colors for percent conservation\n",
  "cmd.set_color('consLT50_grey', [0.851, 0.851, 0.851])\n",
  "cmd.set_color('cons50_69_red', [0.6471, 0.05882, 0.08235])\n",
  "cmd.set_color('cons70_79_mRed', [0.9373, 0.2314, 0.1725])\n",
  "cmd.set_color('cons80_89_lBlue', [0.6196, 0.7922, 0.8824])\n",
  "cmd.set_color('cons90_99_mBlue', [0.2588, 0.5725, 0.7765])\n",
  "cmd.set_color('cons100_dBlue', [0.03137, 0.31765, 0.61176])\n",
  "\n",
  "# Color conserved waters based on percent conservation\n",
  "cmd.color('consLT50_grey', 'ConservedWaters and (q<0.50)')\n",
  "cmd.color('cons50_69_red', 'ConservedWaters and (q>0.49999 and q<0.70)')\n",
  "cmd.color('cons70_79_mRed', 'ConservedWaters and (q>0.69999 and q<0.80)')\n",
  "cmd.color('cons80_89_lBlue', 'ConservedWaters and (q>0.79999 and q<0.90)')\n",
  "cmd.color('cons90_99_mBlue', 'ConservedWaters and (q>0.89999 and q<1.0)')\n",
  "cmd.color('cons100_dBlue', 'ConservedWaters and (q=1.00)')\n",
  "\n",
  file = "WaterAnalysis_ConservedWaters_PASSED_PyMOL_global.pml"
)
