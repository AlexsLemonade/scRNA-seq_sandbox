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
# -m results/GSE84465_meta.tsv
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

opt$data <- "patel_data/salmon_quants/patel_data_counts.tsv"
opt$meta <- "patel_data/GSE57872_meta.tsv"
#--------------------Import Salmon/tximport gene matrix----------------------#
# Import gene expression matrix data
tx.counts <- readr::read_tsv(opt$data) %>% as.data.frame()

# Read in meta data with batch info if it is given
if (!is.na(meta)) {
  meta <- readr::read_tsv(opt$meta)
  if(is.null(opt$var)) { # if no column is specified, stop.
    stop("If normalizing by batch info, you need to specify the column 
       with batch info with -v option.")
  } else {
    # Retrieve the batch info based on the columname
    batch.info <- meta[ , match(opt$var, colnames(meta))]
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

# Save this scnorm data to tsv file
readr::write_tsv(scnorm, file = file.path(opt$output,
                                 paste0(opt$label, "_scnorm_normalize.tsv")))
