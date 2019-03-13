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
# "-r" - Store gene matrix as an RDS file instead of a tsv file. This is
#        advisable for particularly large datasets.
#        
# Filter options:
# Note: if any filter option is used, GeneMatrixFilter will be used for filtering,
# and filter options not specified will be set to default. If you want all filters
# to be automatically the default settings, use -f alone. 
# 
# "-f" - Option to use all default filters to be used as arguments to be used in 
#        data_prep_functions.R's GeneMatrixFilter function. If used: 
#        min_counts = 0.0001, num_samples = 3, num_genes = 50. Note that using this 
#        option makes the last three options listed here irrelevant, as they will
#        be override with the default options. 
# "-m" - Optional cutoff for number of counts for a data point to be considered 
#        'expressed'. Argument to be used in data_prep_functions.R's 
#        GeneMatrixFilter function.
# "-p" - Optional gene filter for number of samples that need to express a given
#        gene. Argument to be used in data_prep_functions.R's GeneMatrixFilter 
#        function.
# "-n" - Optional sample filter for number of genes that a particular sample 
#        needs to express to be kept. Argument to be used in data_prep_functions.R's 
#        GeneMatrixFilter function. 
#        
# Command line example:
#
# Rscript scripts/10x-pre-processing/1-make_gene_matrix-alevin.R \
#   -d data/alevin_output \
#   -o data \
#   -q qc_reports \
#   -l pbmc \
#   -m 1
#   -r

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
              metavar = "character"),
  make_option(opt_str = c("-r", "--rds"), action = "store_true",
              default = FALSE, help = "Store gene matrix as an RDS file instead
              of a tsv file. This is advisable for particularly large datasets."),
  make_option(opt_str = c("-f", "--filter_default"), action = "store_true",
              default = FALSE, help = "Option to use all default filters to be 
              used as arguments to be used in data_prep_functions.R's 
              GeneMatrixFilter function."),
  make_option(opt_str = c("-m", "--min_counts"), type = "numeric",
              default = NA, help = "Optional cutoff for number of counts for 
              a data point to be considered 'expressed'. Argument to be used in 
              data_prep_functions.R's GeneMatrixFilter function. Default is 0.0001"),
  make_option(opt_str = c("-s", "--num_samples"), type = "numeric",
              default = NA, help = "Optional gene filter for number of samples that 
              need to express a given gene. Argument to be used in 
              data_prep_functions.R's GeneMatrixFilter function. Default is 3."),
  make_option(opt_str = c("-g", "--num_genes"), type = "numeric",
              default = NA, help = "Optional sample filter for number of genes that 
              a particular sample needs to express to be kept. Argument to be used in 
              data_prep_functions.R's GeneMatrixFilter function. Default is 50.")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

opt$dir <- "tab_mur_data/alevin_output"
opt$output <- "tab_mur_data/normalized_tab_mur"
opt$qc <- "tab_mur_data/alevinqc_results"
opt$label <- "tab_mur" 
opt$rds <- TRUE

# Add an underscore if opt$label is being used
if (opt$label != "") {
  opt$label <-  paste0("_", opt$label)
}

# Make normalized data folder if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output)
}

#---------------------------Set up filter options------------------------------#
# Put all the filter options together
filt.opts <- c(opt$min_counts, opt$num_samples, opt$num_genes)

# Make vector of default filters
default.filts <- c(0.0001, 3, 50)

# If dfault filter option is used, make filt opts, the defaults
if (opt$filter_default) {
  filt.opts <- default.filts
  message("Default filter settings applied.")
} 

# Check on filter options
if (any(!is.na(filt.opts))) {
  # Use default settings for filters not specified:
  filt.opts[is.na(filt.opts)] <- default.filts[is.na(filt.opts)]
  
  # Before going further, if filters are being used, check for function source script
  data.prep <- file.path("scripts", "util", "data_prep_functions.R")
  if (!file.exists(data.prep)) {
    warning("Filter option was applied but can't find ", data.prep,
            "Make sure to run from the top directory and make sure data_prep_functions.R didn't move.")
  } else {
    source(data.prep)  
  }
  # Print out the filters being used
  message(cat("Filters being used:",
              "\n Minimum count considered detection:", filt.opts[1],
              "\n Number of samples expressing a gene:", filt.opts[2],
              "\n Number of genes expressed by a sample:", filt.opts[3]
              ))
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
                                       progress = FALSE, guess_max = 100000)

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

# Apply gene matrix filter if the options have been provided
if (all(!is.na(filt.opts))) {
  gene.matrix <- GeneMatrixFilter(gene.matrix, min_counts = filt.opts[1],
                                  num_samples = filt.opts[2],
                                  num_genes = filt.opts[3])
  # Report gene matrix dimensions:
  message(cat("\n Original number of genes: ", dim(all.data)[1],
              "\n Original number of cells: ", dim(all.data)[2] - 1,
              "\n Filtered set number of genes:", dim(gene.matrix)[1],
              "\n Filtered set number of cells:", dim(gene.matrix)[2] - 1))
}

# If opt$rds is used, save it as an RDS file, otherwise, default is save as tsv
if (!opt$rds) {
  # Save this overall gene matrix to a tsv file
  readr::write_tsv(gene.matrix, file.path(opt$output, paste0("counts", opt$label,
                                                             ".tsv")))
} else {
  saveRDS(gene.matrix, file.path(opt$output, paste0("counts", opt$label,
                                                    ".RDS")))
}

# Save sample key to a tsv file
readr::write_tsv(data.frame(samples), file.path(opt$output,
                                                paste0("sample_key",
                                                       opt$label, ".tsv")))
