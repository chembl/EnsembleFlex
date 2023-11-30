#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(optparse)
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

## Binding site identification
#---------

# create sub-directory
dir.create(file.path(outdir, "structures_labeled_binding_site"))

# create empty table to save binding site data for all structures
residue_table  <- data.frame(PDB_ID = character(), Amount = numeric(), ResNames = character(), stringsAsFactors = FALSE)

for(i in 1:length(file_list)){
  filename = file_list[[i]]
  #  pdbID <- tools::file_path_sans_ext(basename(filename))
  pdbID <- strsplit((basename(filename)), '[.]')[[1]][1]

  # read the pdb file
  one_pdb <- read.pdb(filename)

  # identify binding site
  # automatically identify 'protein' and 'ligand'
  # be aware that everything that is not protein, nucleic acid or water will be considered as ligand.
  bs <- binding.site(one_pdb, cutoff = opt$distance)

  # save residue names of identified binding site
  #  print(length(bs$resnames))
  #  print(bs$resnames)
  residue_table[i,] <- list(pdbID, length(bs$resnames), toString(bs$resnames))

  # use b-factor column to store interface in PDB file
  one_pdb$atom$b[ bs$inds$atom ] <- 1
  one_pdb$atom$b[ -bs$inds$atom ] <- 0

  # safe as pdb file
  write.pdb(one_pdb, file=paste(outdir,"/structures_labeled_binding_site/",pdbID,"-interface.pdb",sep=''))
}

# safe dataframe
write.csv(residue_table, "binding_site_residues.csv", row.names=FALSE, quote=FALSE)


# collapse results from all structures using "," -> split at "," -> unlist list of lists -> trim leading and trailing whitespace
binding_residue_list_all <- residue_table$ResNames %>% paste(collapse = ",") %>% strsplit(",") %>% unlist() %>% trimws()
# cut the last 4 characters (the chain ID in brackets and a space)
binding_residue_list_all <- gsub('.{4}$', '', binding_residue_list_all)

residue_occurance <- rev(sort(table(binding_residue_list_all)))

# as data frame
residue_occurance_df <- data.frame(residue_occurance)
colnames(residue_occurance_df)[1] <- "Residue"

binding_residue_list_unique <- residue_occurance_df$Residue


#----------
# Bar plots / Histograms
#----------

number_of_structures <- as.character(length(file_list))

# Fitting Labels
#pdf(file="Histogram_binding_residues_occurance.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_occurance.png", width=8, height=5, units="in", res=150)
par(las=2) # make label text perpendicular to axis
#par(mar=c(5,5,3,1)) # adjust margins.
barplot(residue_occurance, main=sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance),
        cex.names=0.8, ylab="# of structures", ylim=c(0, 1.1*max(residue_occurance)))
#barplot(rev(residue_occurance), main="Residues implied in ligand binding", horiz=TRUE, cex.names=0.8)
dev.off()

residue_occurance_freq <- mutate(residue_occurance_df, Percentage=Freq/length(file_list)*100)
# Safe dataframe
write.csv(residue_occurance_freq, "binding_site_residue_occurance.csv", row.names=FALSE, quote=FALSE)

#pdf(file="Histogram_binding_residues_percentage.pdf", width=8, height=5)
png(filename="Histogram_binding_residues_percentage.png", width=8, height=5, units="in", res=150)
par(las=2) # make label text perpendicular to axis
#par(mar=c(5,5,3,1)) # adjust margins.
barplot(residue_occurance_freq$Percentage, names.arg=residue_occurance_freq$Residue,
        main=sprintf("Residues implied in ligand binding in %s structures (cutoff=%sA)",number_of_structures,opt$distance),
        cex.names=0.8, ylab="Frequency [%]", ylim=c(0, 100))
#barplot(rev(residue_occurance), main="Residues implied in ligand binding", horiz=TRUE, cex.names=0.8)
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
pdb <- read.pdb(file_list[[1]])

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

