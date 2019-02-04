# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: After Alevin has been run successfully on your samples, this script
# will turn all your samples into one gene matrix tsv file and run qc on each sample.

# Options:
# "-d" - Directory of where individual samples' alevin folders are located.
# "-o" - Directory of where the output gene matrix file should go.(Optional)
# "-q" - Directory of where the qc reports from alevinQC should go. If no
#        directory is specified, alevinQC is not run. (Optional)
# "-l" - Optional label to add to output files. Generally necessary if processing
#        multiple datasets in the same pipeline.
#
# Command line example:
#
# Rscript scripts/10x-pre-processing/1-make_gene_matrix-alevin.R \
#   -d data/alevin_output \
#   -o data \
#   -q qc_reports \
#   -l "pbmc"

#-------------------------- Get necessary packages-----------------------------#
# Attach needed libraries
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
  make_option(opt_str = c("-d", "--dir"), type = "character", default = getwd(),
              help = "Directory where alevin output folders are located",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character",
              default = getwd(), help = "Directory where you would like the
              output to go", metavar = "character"),
  make_option(opt_str = c("-q", "--qc"), type = "character",
              default = NULL, help = "Directory of where the qc reports from
              alevinQC should go. If no directory is specified, alevinQC is
              not run. (Optional)", metavar = "character"),
  make_option(opt_str = c("-l", "--label"), type = "character",
              default = "", help = "Optional label for output files",
              metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Add an underscore if opt$label is being used
if (opt$label != "") {
  opt$label <-  paste0("_", opt$label)
}

#----------------- Set up function from the COMBINE lab------------------------#
# This function is adapted from:
# https://combine-lab.github.io/alevin-tutorial/2018/running-alevin/
# Set up function from the COMBINE lab
ReadAlevin <- function(base_path = NULL){
  if (! dir.exists(base_path )){
    stop("Directory provided does not exist")
  }
  # Obtain paths to data
  barcode_loc <- file.path(base_path, "alevin", "quants_mat_rows.txt")
  gene_loc <- file.path(base_path, "alevin", "quants_mat_cols.txt")
  data_loc <- file.path(base_path, "alevin", "quants_mat.csv")

  if (!file.exists(barcode_loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene_loc)) {
    stop("Gene name file missing")
  }
  if (!file.exists(data_loc)) {
    stop("Expression matrix file missing")
  }
  # Read in the data from Alevin output
  expression_matrix <- readr::read_csv(data_loc, col_names = FALSE,
                                       progress = FALSE)

  # Transpose the matrix so it is gene x cell
  expression_matrix <- t(expression_matrix[, 1:ncol(expression_matrix)-1])

  # Apply the colnames and rownames to the dataset
  colnames(expression_matrix) <- readLines(barcode_loc)
  rownames(expression_matrix) <- readLines(gene_loc)

  # Make NA values into 0's
  expression_matrix[is.na(expression_matrix)] <- 0
  return(expression_matrix)
}

#------------------------------Run on each sample file-------------------------#
# Obtain sample list
alevin.files <- dir(opt$dir, full.names = TRUE)

# If the output directory for the qc reports, doesn't exist, make one
if (!is.null(opt$qc)) {
  if (!dir.exists(opt$qc)) {
    message(paste0("Can't find '", opt$qc,
                   "' in the current directory, making a directory."))
    dir.create(opt$qc)
  }
  # Produce a QC report
  lapply(alevin.files, function(file) {
         alevinQC::alevinQCReport(file,
                                  sampleId = basename(file),
                                  outputFile = paste0(basename(file),
                                                      "_qc_report.html"),
                                  outputFormat = "html_document",
                                  outputDir = opt$qc)
  })
}

# Run all the data and make it into one big matrix
all.data <- do.call("cbind", lapply(alevin.files, function(file) {
                    # Run this function on our files
                    alv.data <- ReadAlevin(file)

                    # We'll keep track of what sample these cells are coming from
                    colnames(alv.data) <- paste0(colnames(alv.data), ":",
                                                 basename(file))

                    # Return the gene matrix
                    return(alv.data)
                    })
)

# Extract sample info for each cell
samples <- stringr::word(colnames(all.data), sep = ":", -1)

# Take out the data and make genes a column so write_tsv will have it
gene.matrix <- data.frame("genes" = rownames(all.data), all.data)

# Save this overall gene matrix to a tsv file
readr::write_tsv(gene.matrix, file.path(opt$output, paste0("counts", opt$label,
                                                           ".tsv")))

# Save sample key to a tsv file
readr::write_tsv(data.frame(samples), file.path(opt$output,
                                                paste0("sample_key",
                                                       opt$label, ".tsv")))
