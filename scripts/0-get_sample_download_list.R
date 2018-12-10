# CCDL ALSF 2018
# C. Savonen 
# CCDL ALSF 2018
# C. Savonen
#
# Get list of all the samples for the given SRP Id. Optionally, refresh a sample
# list to not include already downloaded files. Makes "files.to.download.txt" file
# which run_pipeline.sh files use for getting the sample IDs to download.

# Options:
# "-i" - the project ID which you wish to download the samples for. eg "SRP079058"
# "-d" - the directory which you would like the downloaded samples to go to.
#        This directory will be created if it doesn't exist, but if it does
#        exist and the refresh option is used, then this directory will be
#        evaluated for what samples have already been downloaded.
# "-r" - If used, the directory given will be evaluated for what samples have
#        been downloaded and those samples will not be added to the list for
#        re-download.
# "-n" - This option is for extracting a random subset of the sample set. If
#        you want to do a test run, put how many samples you would like to be
#        randomly selected. "-n 20" for downloading 20 samples for example.
# "-q" - put the path to the GEO metadata sql, if the file already exists.
#        The sql file will be downloaded if this option isn't given.

# Example usage in bash:
# Rscript scripts/0-get_sample_download_list.R \
# -i SRP079058 \
# -d darmanis_data/salmon_quants \
# -q ref_files/SRAmetadb.sqlite \
# -r

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
              help = "Use this option if you don't want to re-process already existing files"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help = "If you'd like to only use a portion of the samples,
              state the max number you'd like to download."),
  make_option(opt_str = c("-q", "--sql"), type = "character", default = NULL,
              help = "Directory path to sql", metavar = "character")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check if the sql file from GEO is where the option says.
if (!is.null(opt$sql)) {
    stop("The path to the GEO sql file has not been specified.")
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
write.csv(files, file.path("SRA.files.csv"))

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

