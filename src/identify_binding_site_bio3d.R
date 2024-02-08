#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(optparse)
library(R.utils) # for function "isAbsolutePath"
library(dplyr) # for data wrangling
library(ggplot2) # for plotting
library(bio3d) # Bio3d - alignment, PCA, NMA


option_list = list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="input directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="outdir",
              help="output directory name [default=%default]", metavar="character"),
  make_option(c("-d", "--distance"), type="double", default=4.0,
              help="distance cutoff in A [default=%default]")
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

## Binding site identification
#---------

# create sub-directory
dir.create(file.path(outdir, "structures_labeled_binding_site"))

# create empty table to save binding site data for all structures
residue_table  <- data.frame(PDB_ID = character(), LigIDs = character(), Amount = numeric(), ResNames = character(), stringsAsFactors = FALSE)
#residue_numbers_table  <- data.frame(ResNumbers = numeric(), stringsAsFactors = FALSE)
residue_numbers_list  <- list()

for(i in 1:length(files)){
  filename = files[[i]]
  #  pdbID <- tools::file_path_sans_ext(basename(filename))
  pdbID <- strsplit((basename(filename)), '[.]')[[1]][1]

  # read the pdb file
  one_pdb <- read.pdb(filename)

  # identify binding site
  # automatically identify 'protein' and 'ligand'
  # be aware that everything that is not protein, nucleic acid or water will be considered as ligand.
  bs <- binding.site(one_pdb, cutoff = opt$distance)

  # Get ligand IDs
  # be aware that everything that is not protein, nucleic acid or water will be considered as ligand.
  lig.inds <- atom.select(one_pdb, "ligand")
  # Access the PDB data with the selected atom indices
  ligIDs <- one_pdb$atom[ lig.inds$atom, "resid" ]

  # save residue names of identified binding site
  #  print(length(bs$resnames))
  #  print(bs$resnames)
  residue_table[i,] <- list(pdbID, paste(unique(ligIDs),collapse=' '), length(bs$resnames), toString(bs$resnames))
  residue_numbers_list[[i]] <- bs$resno

  # use b-factor column to store interface in PDB file
  one_pdb$atom$b[ bs$inds$atom ] <- 1
  one_pdb$atom$b[ -bs$inds$atom ] <- 0

  # safe as pdb file
  write.pdb(one_pdb, file=paste(outdir,"/structures_labeled_binding_site/",pdbID,"-interface.pdb",sep=''))
}

# safe dataframe
# write.csv2() uses a comma (“,”) for the decimal point and a semicolon (“;”) for the separator.
write.csv(residue_table, "binding_site_residues.csv", row.names=FALSE, quote=FALSE)
write.table(residue_table, "binding_site_residues.tsv", row.names=FALSE, quote=FALSE, sep = "\t")

# collapse results from all structures using "," -> split at "," -> unlist list of lists -> trim leading and trailing whitespace
binding_residue_list_all <- residue_table$ResNames %>% paste(collapse = ",") %>% strsplit(",") %>% unlist() %>% trimws()
# cut the last 4 characters (the chain ID in brackets and a space)
binding_residue_list_all <- gsub('.{4}$', '', binding_residue_list_all)
# sort by frequency and store table
residue_occurance <- rev(sort(table(binding_residue_list_all)))
# store as data frame
residue_occurance_df <- data.frame(residue_occurance)
colnames(residue_occurance_df)[1] <- "Residue"

binding_residue_list_unique <- residue_occurance_df$Residue


#----------
# Bar plots / Histograms
#----------

number_of_structures <- as.character(length(files))

residue_occurance_freq <- mutate(residue_occurance_df, Percentage=Freq/length(files)*100)
# Safe dataframe
write.csv(residue_occurance_freq, "binding_site_residue_occurance_frequency.csv", row.names=FALSE, quote=FALSE)


# # Fitting Labels
#pdf(file="Histogram_binding_residues_frequency.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_frequency.png", width=8, height=5, units="in", res=150)
par(las=2) # make label text perpendicular to axis
#par(mar=c(5,5,3,1)) # adjust margins.
barplot(residue_occurance_freq$Freq, names.arg=residue_occurance_freq$Residue,
        main=sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance),
        cex.names=0.8, ylab="# of structures", ylim=c(0, 1.1*max(residue_occurance_freq$Freq)))
#barplot(rev(residue_occurance), main="Residues implied in ligand binding", horiz=TRUE, cex.names=0.8)
dev.off()

#pdf(file="Histogram_binding_residues_percentage.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_percentage.png", width=8, height=5, units="in", res=150)
par(las=2) # make label text perpendicular to axis
#par(mar=c(5,5,3,1)) # adjust margins.
barplot(residue_occurance_freq$Percentage, names.arg=residue_occurance_freq$Residue,
        main=sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance),
        cex.names=0.8, ylab="Frequency [%]", ylim=c(0, 100))
#barplot(rev(residue_occurance), main="Residues implied in ligand binding", horiz=TRUE, cex.names=0.8)
dev.off()

#pdf(file="Histogram_binding_residues_frequency_colored.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_frequency_colored.png", width=8, height=5, units="in", res=150)
ggplot(residue_occurance_freq, aes(x=Residue, y=Freq, fill=Freq)) +
  geom_bar(stat="identity") +
  ggtitle(sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance)) +
  ylab("Frequency [count]") +
  labs(fill='Frequency [count]') +
  xlab("Residue") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_continuous()
# , colors = c('firebrick','orange','green','lightblue','navy')
dev.off()

#pdf(file="Histogram_binding_residues_percentage_colored.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_percentage_colored.png", width=8, height=5, units="in", res=150)
ggplot(residue_occurance_freq, aes(x=Residue, y=Percentage, fill=Percentage)) +
  geom_bar(stat="identity") +
  ggtitle(sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance)) +
  ylab("Frequency [%]") +
  labs(fill='Frequency [%]') +
  xlab("Residue") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_continuous(limits=c(0, 100))
# , colors = c('firebrick','orange','green','lightblue','navy')
dev.off()


#----------
# Save binding site occurence frequency in reference pdb structure b-factor column
#----------

# get residue numbers sorted by occurancy
binding_residue_num <- gsub("[A-Z]", "", toString(binding_residue_list_unique)) %>% strsplit(",") %>% unlist() %>% trimws() %>% as.integer()
# safe binding_residue_num to file
# to be used for superimposing only on binding site residues
write.table(binding_residue_num, "binding_site_residue_numbers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# read the reference pdb file (first file in filelist)
pdb <- read.pdb(files[[1]])

#--------------
## Save residue occurance in b-factor column on reference structure

pdb_bindingsite_as_b_occ <- pdb
# get indices of all bindingsite residues
binding.inds <- atom.select(pdb_bindingsite_as_b_occ, resno=c(binding_residue_num))

# use b-factor column to store interface in PDB file
pdb_bindingsite_as_b_occ$atom$b <- 0
pdb_bindingsite_as_b_occ$atom$b[ binding.inds$atom ] <- 1
# write to file
write.pdb(pdb_bindingsite_as_b_occ, file=paste(outdir,"/binding_site_interface_labelled_occurance.pdb",sep=''))

#--------------
## Save residue occurance frequency in b-factor column on reference structure

pdb_bindingsite_as_b_freq <- pdb
# get indices of all bindingsite residues
binding.inds <- atom.select(pdb_bindingsite_as_b_freq, resno=c(binding_residue_num))

# use b-factor column to store interface in PDB file
pdb_bindingsite_as_b_freq$atom$b <- 0
#pdb_bindingsite_as_b_freq$atom$b[ binding.inds$atom ] <- 1
for(i in 1:length(binding_residue_num)){
  resnum = binding_residue_num[i]
  resnum.inds <- atom.select(pdb_bindingsite_as_b_freq, resno=resnum)
  pdb_bindingsite_as_b_freq$atom$b[ resnum.inds$atom ] <- as.numeric(residue_occurance_df$Freq[i])
}
# write to file
write.pdb(pdb_bindingsite_as_b_freq, file=paste(outdir,"/binding_site_interface_labelled_frequency.pdb",sep=''))

#--------------
## Save residue occurance frequency as percentage in b-factor column on reference structure

pdb_bindingsite_as_b_percent <- pdb
# get indices of all bindingsite residues
binding.inds <- atom.select(pdb_bindingsite_as_b_percent, resno=c(binding_residue_num))

# use b-factor column to store interface in PDB file
pdb_bindingsite_as_b_percent$atom$b <- 0
for(i in 1:length(binding_residue_num)){
  resnum = binding_residue_num[i]
  resnum.inds <- atom.select(pdb_bindingsite_as_b_percent, resno=resnum)
  pdb_bindingsite_as_b_percent$atom$b[ resnum.inds$atom ] <- as.numeric(residue_occurance_freq$Percentage[i])
}
# write to file
write.pdb(pdb_bindingsite_as_b_percent, file=paste(outdir,"/binding_site_interface_labelled_percentage.pdb",sep=''))



# Construct Matrix of Ligand-Residue interactions
# (with one row per structure and one column per residue)
# transform list of residue numbers to sparse binary matrix with 1s at residue number indices
l <- residue_numbers_list
unlist_l <- unlist(l)
M <- matrix(0, nrow = length(l), ncol = max(unique(unlist_l)))
ij <- cbind(rep(1:length(l), lengths(l)), unlist_l)
M[ij] <- 1
# remove all columns with only 0s
M <- M[, colSums(abs(M)) != 0]
#M

# Heatmap of Ligand-Protein interactions
png("interaction_heatmap.png", units="in", width=9, height=9, res=300)
heatmap(M, main=sprintf("Ligand-Protein interactions in %s structures (cutoff=%sA)",number_of_structures,opt$distance),
    labRow = residue_table$LigIDs,
    labCol = sort(unique(unlist_l)),
    xlab="Residue", ylab="Ligand",
    Colv = NA, Rowv = NA, scale="none") # no clustering, no scaling
dev.off()
print("Plot saved to file interaction_heatmap.png")
