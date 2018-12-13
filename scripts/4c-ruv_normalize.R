# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: Using scnorm for post-processing of scRNA-seq Salmon processed data
#
# Options:
# "-d" - Directory of where individual samples' salmon folders are located.
# "-o" - Directory of where the output gene matrix tsv file should go.(Optional)
# "-l" - Optional label to add to output files. Generally necessary if processing
#        multiple datasets in the same pipeline.
# 
# Command line example:
#
# Rscript scripts/4b-scnorm_normalize.R \
# -d data/salmon_quants \
# -m results/GSE84465_meta.tsv plate.id
# -o data \
# -l "patel"

#-------------------------- Get necessary packages-----------------------------#
# Install RUVnormalize if it isn't
if (!("RUVnormalize" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RUVnormalize", suppressUpdates = FALSE)
}
library(RUVnormalize)
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
                help = "Path to gene matrix in tsv format, gene x sample",
                metavar = "character"),
    make_option(opt_str = c("-m", "--meta"), type = "character", default = getwd(),
                help = "Path to metadata to use for normalization in tsv format,
                followed by a space and the name of the column with batch info",
                metavar = "character"),
    make_option(opt_str = c("-o", "--output"), type = "character",
                default = getwd(), help = "Directory where you would like the
                output to go", metavar = "character"),
    make_option(opt_str = c("-l", "--label"), type = "character",
                default = "", help = "Optional label for output files",
                metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

#--------------------Import Salmon/tximport gene matrix----------------------#
# Import gene expression matrix data
tx.counts <- readr::read_delim(opt$data)

# Metadata info
meta <- unlist(strsplit(" ", opt$meta))

# Get the column name
meta.col <- meta[2]

# Read in meta data
meta <- readr::read_delim(meta)

# Retrieve colname 
batch.info <- meta[ , match(meta.col, colnames(meta))]



# Save this scnorm data to tsv file
readr::write_tsv(seurat, file = file.path(opt$output,
                                 paste0(opt$label, "scnorm_normalize.tsv")))
