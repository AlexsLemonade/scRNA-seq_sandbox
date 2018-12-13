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
              default = NULL, help = "SRR Run sample ID of the RNA-seq
              data you would like to download.", metavar = "character"),
  make_option(opt_str = c("-d", "--dir"), type = "character",
              default = getwd(), help = "Directory where you would
              like the data downloaded to go", metavar = "character"),
  make_option(opt_str = c("-q", "--sql"), type = "character",
              default = NULL, help = "Directory path to sql",
              metavar = "character")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Get sra path
srafile <- file.path(opt$sql)

# Connect to SQL database
con <- dbConnect(RSQLite::SQLite(), srafile)

# Download each sample
getFASTQfile(opt$sample, con, destDir = opt$dir))
