# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: Using scnorm for post-processing of scRNA-seq Salmon processed data
#
# Options:
# "-d" - Directory of where individual samples' salmon folders are located.
# "-m" - Supply metadata .tsv file path that contains batch information (optional).
# "-v" - The column name in the metadata file that contains the batch information
#        you would like to normalize by.
# "-o" - Directory of where the output gene matrix tsv file should go.(Optional)
# "-l" - Optional label to add to output files. Generally necessary if processing
#        multiple datasets in the same pipeline.
# 
# Command line example:
#
# Rscript scripts/4b-scnorm_normalize.R \
# -d data/salmon_quants \
# -m results/GSE84465_meta.tsv \
# -v plate.id \
# -o data \
# -l "patel"

#-------------------------- Get necessary packages-----------------------------#
# Install SCnorm if it isn't
if (!("SCnorm" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("SCnorm", suppressUpdates = FALSE)
}
library(SCnorm)
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
                help = "Path to gene matrix in tsv format, gene x sample",
                metavar = "character"),
    make_option(opt_str = c("-m", "--meta"), type = "character", default = NA,
                help = "Path to metadata to use for normalization in tsv format,
                followed by a space and the name of the column with batch info",
                metavar = "character"),
    make_option(opt_str = c("-v", "--var"), type = "character", default = NA,
                help = "Column name in metadata of variable that has the batch 
                info to be used for normalization.", metavar = "character"),
    make_option(opt_str = c("-o", "--output"), type = "character",
                default = getwd(), help = "Directory where you would like the
                output to go", metavar = "character"),
    make_option(opt_str = c("-l", "--label"), type = "character",
                default = "", help = "Optional label for output files",
                metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
opt$data <- "darmanis_data/salmon_quants/darmanis_data_counts.tsv"
opt$meta <- "darmanis_data/GSE84465_meta.tsv"
opt$var <- "plate.id.ch1"
#--------------------Import Salmon/tximport gene matrix----------------------#
# Import gene expression matrix data
tx.counts <- read.table(opt$data, header = TRUE, stringsAsFactors = FALSE)

# Read in meta data with batch info if it is given
if (!is.na(opt$meta)) {
  meta <-  read.table(opt$meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if(is.null(opt$var)) { # if no column is specified, stop.
    stop("If normalizing by batch info, you need to specify the column 
       with batch info with -v option.")
  } else {
    # Retrieve the batch info based on the column name
    col.index <- match(opt$var, colnames(meta))
    if (!is.na(col.index)) {
      # If the column name matches, extract it    
      batch.info <- meta[ , match(opt$var, colnames(meta))]
      batch.info <- batch.info[match(colnames(tx.counts)[-1],meta$geo_accession)]
    } else {
      stop("Column name specified with -v does not match any column in metadata")
    }
  }
} else {
  # Set the batch info as all one batch. (Not ideal, but sometimes you don't
  # have batch info so this appears to be what SCnorm recommends.)
  batch.info <- rep(1, ncol(tx.counts) - 1)
}

# SCnorm wants the genes as the rownames
rownames(tx.counts) <- tx.counts$gene

# Do the normalization 
scnorm <- SCnorm(Data = tx.counts[, -which(colnames(tx.counts) == "gene")],
                 Conditions = batch.info, FilterCellNum = 10, NCores=3)

# Extract the data from the object
scnorm.data <- results(scnorm)

# Save this scnorm data to tsv file
readr::write_tsv(scnorm.data, file.path(opt$output,
                                 paste0(opt$label, "_scnorm_normalize.tsv")))
