# CCDL ALSF 2018
# C. Savonen
#
# Purpose: This script is used in a bash loop to download file(s) from command line.
#
# Options:
#  "-s" Run sample ID of the RNA-seq data you would like to download.
#  "-d" Path to directory where you would like the data downloaded to go
#  "-q" Path to downloaded metadata SQL file from GEO. Usually file is
#       called "SRAmetadb.sqlite"
#
# Example usage in bash:
#
# Rscript scripts/1-download_sra.R \
# -s SRR3934349 \
# -q ref_files/SRAmetadb.sqlite \
# -d raw_data

# Using SRAdb to download the file and optparse to determine options from bash use
library(optparse)
library(SRAdb)

# Get options from command line
option_list <- list( 
  make_option(opt_str = c("-s", "--sample"), type = "character",
              default = NULL, help = "SRR Run Sample ID of the RNA-seq
              data you would like to download. For multiple samples,
              separate ids by a comma and no spaces eg. SRR3934349,
              SRR3934348", metavar = "character"),
  make_option(opt_str = c("-d", "--dir"), type = "character",
              default = getwd(), help = "Directory where you would
              like the data downloaded to go", metavar = "character"),
  make_option(opt_str = c("-q", "--sql"), type = "character",
              default = NULL, help = "Directory path to sql",
              metavar = "character")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Separate sample ids
samples <- strsplit(opt$sample, ",")

# Get sra path
srafile <- file.path(opt$sql)

# Connect to SQL database
con <- dbConnect(RSQLite::SQLite(), srafile)

# Download each sample
lapply(samples, function(x) getFASTQfile(x, con, destDir = opt$dir))
