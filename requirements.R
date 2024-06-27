#!/usr/bin/env Rscript

# # install required R packages (also available via conda channel)
# install.packages("bio3d", dependencies=TRUE)
# install.packages("optparse", dependencies=TRUE)
#  if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install("msa")
# install.packages("ggplot2", dependencies=TRUE)
# install.packages("pheatmap", dependencies=TRUE)

# install required R packages (not available via conda channel)
# if setting CRAN mirror is required use the below.
# It will use the cloud mirror which automatically serves the files from a location near you.

# install.packages("devtools", dependencies=TRUE, repos = c(CRAN = "https://cloud.r-project.org"))
# library(devtools)
# devtools::install_bitbucket("Grantlab/bio3d", subdir = "bio3d-core", ref="core", dependencies=TRUE)
# devtools::install_bitbucket("Grantlab/bio3d-eddm", dependencies=TRUE)

#install.packages("remotes", dependencies=TRUE, repos = c(CRAN = "https://cloud.r-project.org"))
library(remotes) # remotes package is lighter than devtools and available on conda channels
remotes::install_bitbucket("Grantlab/bio3d/bio3d-core", ref="core", dependencies=TRUE)
remotes::install_bitbucket("Grantlab/bio3d-eddm", dependencies=TRUE)