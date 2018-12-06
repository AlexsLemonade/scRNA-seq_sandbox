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
  make_option(c("-r", "--refresh"), action ="store_true", default = FALSE,
              help="Use this option if you don't want to re-process already existing files"),
  make_option(c("-n", "--number"), type = "numeric", default = TRUE,
              help="If you'd like to only use a portion of the samples,
              state the max number you'd like to download.")
)
opt <- parse_args(OptionParser(option_list = option_list))
dat.dir <- file.path(opt$dir)

#------------------- Connect to NCBI's SRA SQL database------------------------#
srafile <- file.path("data", "SRAmetadb.sqlite")
if (!file.exists(srafile)) {
    getSRAdbFile()
}
con <- DBI::dbConnect(RSQLite::SQLite(), srafile)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile(opt$id, con)
write.csv(files, file.path("data", "all.SRA.files.csv"))

# If we want to restrict the number of samples being processed:
if (!is.null(opt$number)){
  set.seed(12345)
  files <- files[sample(nrow(files), opt$number), ]
}
#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists(dat.dir)){
  dir.create(dat.dir)
} else {
  if (opt$refresh) {
    # Determine which files have already been downloaded and remove them from the list
    existing.files <- dir("data/salmon_quants")
    
    # Filter them out of the file list. 
    existing.files <- match(files$run, existing.files)
    files <- files[is.na(existing.files), ]
  }
}
# Write the table of the urls
write.table(files$run, file.path("files.2.download.txt"), col.names = FALSE,
            row.names = FALSE, sep = "\n", quote = FALSE)

# Can download them direct from here, but I wouldn't recommend doing this unless 
# you have a small number of samples or a large amount of space on your computer
# getFASTQfile( , con, destDir = dat.dir)

