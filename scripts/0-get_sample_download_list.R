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
              help = "Directory where you would like the data downloaded to go",
              metavar = "character"),
  make_option(c("-r", "--refresh"), action ="store_true", default = FALSE,
              help="Use this option if you don't want to re-process already existing files"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help="If you'd like to only use a portion of the samples,
              state the max number you'd like to download."),
  make_option(opt_str = c("-q", "--sql"), type = "character", default = NULL,
              help = "Directory path to sql", metavar = "character")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check if the sql file from GEO is where the option says.
if (!is.null(opt$sql)) {
  if (!file.exists(opt$sql)) {
      stop("The path to the GEO sql file does not exist. Double check that the file
       is at the path you specified relative to the current directory.")
  }
}
#------------------- Connect to NCBI's SRA SQL database------------------------#
# Get sra path
srafile <- file.path(opt$sql)
if (!file.exists(opt$sql)) {
    getSRAdbFile()
}
con <- DBI::dbConnect(RSQLite::SQLite(), opt$sql)

# Get a list of the samples associated with the project we are interested in
files <- listSRAfile(opt$id, con)
write.csv(files, file.path("data", "SRA.files.csv"))

# If we want to restrict the number of samples being processed:
if (!is.null(opt$number)){
  set.seed(12345)
  files <- files[sample(nrow(files), opt$number), ]
}
#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists(opt$dir)){
  dir.create(opt$dir)
} else {
  if (opt$refresh) {
    # Determine which files have already been downloaded and remove them from the list
    existing.files <- dir(opt$dir)
    
    # Filter them out of the file list. 
    existing.files <- match(files$run, existing.files)
    files <- files[is.na(existing.files), ]
  }
}
# Write the table of the urls
write.table(files$run, file.path("files.to.download.txt"), col.names = FALSE,
            row.names = FALSE, sep = "\n", quote = FALSE)

# Can download them direct from here, but I wouldn't recommend doing this unless 
# you have a small number of samples or a large amount of space on your computer
# getFASTQfile( , con, destDir = dat.dir)

