# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: Get metadata as a tsv using GEOquery. To be used for some
# types of normalization methods and also differential expression analyses.

# Options:
# "-g" - GEO ID for the dataset so the metadata can be downloaded.
# "-o" - Directory of where the output gene matrix RDS file should go.(Optional)
# 
# Command line example:
#
# Rscript scripts/4-metadata_setup.R \
# -g GSE84465 \
# -o data 


#-------------------------- Get necessary packages-----------------------------#
if (!("GEOquery" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GEOquery", suppressUpdates = TRUE)
}

# Attach needed libraries
library(GEOquery)
library(optparse)

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(opt_str = c("-g", "--gse"), type = "character", default = "",
                help = "GSE Id(s) of the metadata you would like to download.
                Separate by a space eg -g GSE84465 GSE57872",
                metavar = "character"),
    make_option(opt_str = c("-o", "--output"), type = "character",
                default = getwd(), help = "Directory where you would like the
                output to go", metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Get the different gse's into a list
gse <- strsplit(opt$gse, " ")

# Download metadata and save to tsv file for each dataset
metadata <- lapply(gse, function(x) {
    
                   # Get the metadata through GEOquery
                   meta <- getGEO(x, destdir = opt$output)
                   
                   # Write to a tsv file
                   readr::write_tsv(meta[[1]]@phenoData@data,
                   file.path(opt$output, paste0(x, "_meta.tsv")))
                   
                   # Print out a preview of the metadata so the user can see
                   # what variables are there
                   print(head(meta[[1]]@phenoData@data))
    })

