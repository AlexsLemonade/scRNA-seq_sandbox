# CCDL ALSF 2018
# C. Savonen 
# 
# Get raw reads for the given SRP Id
#
#-------------------------- Get necessary packages-----------------------------#
library(SRAdb)
library(DBI)
library(optparse)

option_list <- list( 
  make_option(opt_str = c("-i", "--id"), type = "character", default = NULL,
              help = "SRP ID of the RNA-seq data you would like to download",
              metavar = "character"),
  make_option(opt_str = c("-d", "--dir"), type = "character", default = getwd(),
              help = "directory where you would like the data downloaded to go",
              metavar = "character"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help="If you'd like to only use a portion of the samples,
              state the max number you'd like to download.")
)
opt <- parse_args(OptionParser(option_list=option_list))
dat.dir <- file.path(opt$dir)

#------------------- Connect to NCBI's SRA SQL database------------------------#
if (!file.exists("SRAmetadb.sqlite") {
    srafile <- getSRAdbFile()
}
con <- dbConnect(RSQLite::SQLite(), srafile)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile(opt$id, con)

# If we want to restrict the number of samples being processed:
if (!is.null(opt$number)){
  set.seed(12345)
  files <- files$run[sample(nrow(files), opt$number)]
} else {
  files <- files$run
}
#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists(dat.dir)){
  dir.create(dat.dir)
}
getFASTQfile(files, con, destDir = dat.dir)
