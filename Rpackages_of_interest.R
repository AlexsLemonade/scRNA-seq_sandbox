# C. Savonen 
# CCDL ALSF 2018
# 
# A collection of data
source("https://bioconductor.org/biocLite.R")
biocLite("scRNAseq") # A collection of data
biocLite("simpleSingleCell") # a different single cell workflow
biocLite("GEOquery")

library(GEOquery)
GSE <- c("GDS4495", "GSE102096", "GSE107405")
data <- lapply(GSE, getGEO)

library(scRNAseq)
