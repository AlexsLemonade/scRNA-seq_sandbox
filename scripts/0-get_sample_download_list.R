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
# "-o" - Directory where you would like other output files (SRA.files.csv and
#        files.to.download.txt to go) default is current directory.(Optional)
# "-q" - Directory path to SRAmetadb.sqlite. The sql file will be downloaded
#        if the file does not exist at the given location.
# "-r" - If used, the directory given will be evaluated for what samples have
#        been downloaded and those samples will not be added to the list for
#        re-download. (Optional)
# "-n" - This option is for extracting a random subset of the sample set. If
#        you want to do a test run, put how many samples you would like to be
#        randomly selected. "-n 20" for downloading 20 samples for example.
#        (Optional)
# "-s" - Can put a number here to set as a seed for random selection of samples
#        Note that this option is irrelevant if you are not using option "-n"
#        (Optional)

# Example usage in bash:
# Rscript scripts/0-get_sample_download_list.R \
# -i SRP079058 \
# -d darmanis_data/salmon_quants \
# -q ref_files/SRAmetadb.sqlite

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
  make_option(opt_str = c("-o", "--output"), type = "character",
              default = getwd(), help = "Directory where you would like
              other output files to go",
              metavar = "character"),
  make_option(opt_str = c("-q", "--sql"), type = "character", default = NULL,
              help = "Directory path to SRAmetadb.sqlite. If the file does not
              exist at the given directory, it will be downloaded",
              metavar = "character"),
  make_option(c("-r", "--refresh"), action ="store_false", default = TRUE,
              help = "Use this option if you want to re-process already
              existing files. Default is to not re-download already
              existing files"),
  make_option(c("-n", "--number"), type = "numeric", default = NULL,
              help = "If you'd like to only use a portion of the samples,
              state the max number you'd like to download."),
  make_option(c("-s", "--seed"), type = "numeric", default = 12345,
              help = "Optional set the seed number for random if you are using
              option -n to randomly select a subset of samples.")
)
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check if the sql file path has been given. If not, stop.
if (is.null(opt$sql)) {
    stop("The path to the GEO's SRAmetadb.sqlite file has not been specified.")
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
write.csv(files, file.path(opt$output, "SRA.files.csv"))

# If we want to restrict the number of samples being processed:
if (!is.null(opt$number)){
  set.seed(opt$seed)
  files <- files[sample(nrow(files), opt$number), ]
}
#-------------------------Get all the FASTQ files------------------------------#
if (!dir.exists(opt$dir)){
  dir.create(opt$dir)
} else {
  if (!opt$refresh) {
    # Determine which files have already been downloaded and remove them from the list
    existing.files <- dir(opt$dir)
    
    # Filter them out of the file list. 
    existing.files <- match(files$run, existing.files)
    files <- files[is.na(existing.files), ]
  }
}
# Write the table of the urls
write.table(files$run, file.path(opt$output, "files.to.download.txt"),
            col.names = FALSE, row.names = FALSE, sep = "\n", quote = FALSE)

# Can download them direct from here, but I wouldn't recommend doing this unless 
# you have a small number of samples or a large amount of space on your computer
# getFASTQfile(files$run, con, destDir = dat.dir)
