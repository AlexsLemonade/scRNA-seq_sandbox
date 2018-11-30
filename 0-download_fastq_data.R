# CCDL ALSF 2018
# C. Savonen 
# 
# Get raw reads for the given SRP Id
# Magrittr pipe
`%>%` <- dplyr::`%>%`
#-------------------------- Get necessary packages-----------------------------#
library(SRAdb)
library(optparse)

option_list <- list( 
  make_option(opt_str = c("-i", "--id"), type = "character", default = NULL,
              help = "SRP ID of the RNA-seq data you would like to download
              eg. SRP079058",
              metavar = "character"),
  make_option(opt_str = c("-d", "--dir"), type = "character", default = getwd(),
              help = "directory where you would like the data downloaded to go",
              metavar = "character"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help="If you'd like to only use a portion of the samples,
              state the max number you'd like to download.")
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$dir <- "data/raw_data"
opt$id <- "SRP079058"
dat.dir <- file.path(opt$dir)

#------------------- Connect to NCBI's SRA SQL database------------------------#
srafile <- file.path("data", "SRAmetadb.sqlite")
if (!file.exists(srafile)) {
    getSRAdbFile()
}
con <- DBI::dbConnect(RSQLite::SQLite(), srafile)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile(opt$id, con)
write.csv(files, file.path("data", "SRA.files.csv"))

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
} else {
  if (is.null(opt$number)) {
    # Get a list of previously downloaded forward sequence files
    existing.for.files <- grep("_1.fastq.gz", dir(dat.dir), value = TRUE) 
    existing.for.files <- gsub("_1.fastq.gz", "", existing.for.files)

    # Get a list of previously downloaded reverse sequence files
    existing.rev.files <- grep("_2.fastq.gz", dir(dat.dir), value = TRUE) 
    existing.rev.files <- gsub("_1.fastq.gz", "", existing.for.files) 
  
    # Files that don't need to be downloaded
    existing.files <- existing.for.files[!is.na(match(existing.for.files,
                                                    existing.rev.files))]
    # Filter them out of the file list. 
    existing.files <- match(files, existing.files)
    files <- files[is.na(existing.files)]
  }
}
getFASTQfile(files, con, destDir = dat.dir)

