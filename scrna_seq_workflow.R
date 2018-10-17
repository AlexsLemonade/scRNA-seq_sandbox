# C. Savonen
# CCDL ALSF 2018
# Purpose: Example workflow for scRNA-seq data
#
# scRNA-seq with Seraut
# 

# Problematic package hdf5r workaround: 
install.packages("remotes")
remotes::install_github("UCSF-TI/fake-hdf5r")

# Install Seraut
install.packages('Seurat')
library(Seurat)


library(devtools)
install_github("YosefLab/scone")