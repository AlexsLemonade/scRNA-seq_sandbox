# CCDL ALSF 2018
# C. Savonen 
# 
# Get raw reads for single cell glioblastoma data from Darmanis et al, 2017
#
#
#-------------------------- Get necessary packages-----------------------------#
library(SRAdb)
library(DBI)
library(optparse)

option_list <- list( 
  make_option(opt_str = c("-d", "--data"), type = "character", default = NULL, 
              help = "SRP ID of the RNA-seq data you would like to download",
              metavar = "character"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help="If you'd like to only use a portion of the samples,
              state the max number you'd like to download.")
)
opt <- parse_args(OptionParser(option_list=option_list))

#------------------- Connect to NCBI's SRA SQL database------------------------#
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(), srafile)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile(opt$data, con)

# If we want to restrict the number of samples being processed:
if (!is.null(opt$number)){
  files <- files$run[runif(opt$number, min = 1, max = nrow(files))]
} else {
  files <- files$run
}
#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists("raw_data")){
  dir.create("raw_data")
}
getFASTQfile(files, con, destDir = "raw_data")
