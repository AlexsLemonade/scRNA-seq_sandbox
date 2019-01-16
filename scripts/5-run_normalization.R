# C.Savonen, CCDL for ALSF
# 2019
#
# Purpose: Normalizing count data from single rna-seq
#
# Adapted from Automated Single-cell Analysis Pipeline (ASAP) github url:
# https://github.com/DeplanckeLab/ASAP/blob/master/R_Python/normalization.R
#
# Options:
# '-d' : Path to gene matrix in tab delimited format, gene x sample with gene
#        info as the first column
# '-a' : Normalization method to use. To run all methods use "all", to run more
#        than one method, separate the names by a space. eg -a tmm log
#        Options:
#           'scale'- scale and center each sample's data using base::scale
#           'log'- just log transform alone
#           'voom'- quantile normalization from the DEseq2::voom function
#           'deseq2'- Use DEseq2 size factors
#           'vsd'- Use DEseq2 varianceStabilizingTransformation function
#           'tmm'- trimmed mean M values by edgeR
#           'all'- run all the above methods
# '-o' : Directory where you would like the output to go. Default is current
#        directory
# '-l' : Optional label for output files"
#
# Command line example:
#
# Rscript scripts/5-run_normalization.R \
# -d data/darmanis_counts.tsv \
# -a log \
# -o normalized_data \
# -l darmanis
#
# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(opt_str = c("-d", "--data"), type = "character", default = "none",
              help = "Path to gene matrix in tab delimited format, gene x sample
              with gene info as the first column",
              metavar = "character"),
  make_option(opt_str = c("-a", "--algorithm"), type = "character",
              default = "none", help = "Normalization method to use. Options:
              'scale', log', 'voom', 'deseq2', 'vsd', 'rld' or 'tmm' ",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character",
              default = getwd(), help = "Directory where you would like the
              output to go", metavar = "character"),
  make_option(opt_str = c("-l", "--label"), type = "character",
              default = "", help = "Optional label for output files",
              metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Stop if no input data matrix is specified
if (opt$data == "none") {
    stop("Error: no specified input gene matrix file. Use option -d to specify
    the gene matrix tsv data you would like to normalize")
}

# Create the output for results folder if it does not exist
if (!dir.exists(opt$output)) {
    message(paste("Output folder:", opt$output, "does not exist, creating one"))
    dir.create(opt$output)
}

#---------------------------Set up algorithms to run---------------------------#
# Check that the dimension reduction option given is supported
all.algorithms <- c( 'scale', 'log', 'voom','tmm', 'deseq2', 'vsd', 'rld')

# If "all" is chosen, run all the normalization methods
if (opt$algorithm == "all") {
    opt$algorithm <- all.algorithms
} else {
    opt$algorithm <- unlist(strsplit(opt$algorithm, " "))
}

# Check that the algorithms chosen are in the supported list
if (any(is.na(match(opt$algorithm, all.algorithms)))) {
    stop("That is not a normalization method supported by this script.
    Check for typos. Acceptable options:'scale', 'log', 'voom','tmm',
    'deseq2', 'vsd', 'rld'")
}
#----------------------------------Load data-----------------------------------#
# Read in a tsv file of data
dataset <- readr::read_tsv(opt$data)

# Separate genes from the numeric data
genes <- dataset$gene
dataset <- dataset[, -1]

#-----------------------Run each algorithm in the list-------------------------#
# For each algorithm, run through this loop
for (algorithm in opt$algorithm) {
    
    # Print out progress message:
    cat("\n \n Running:", algorithm, "...")
    
    # Run normalization algorithms
    if (algorithm  == "scale"){
        data.out <- as.data.frame(scale(dataset))
        title <- "Scaled Expression"
        
    } else  if (algorithm == "log") {
        data.out <- sign(dataset) * log2(1 + abs(dataset))
        title <- "Log2 Expression"
        
    } else  if (algorithm == "voom") {
        data.out <- as.data.frame(limma::voom(counts = dataset,
        normalize.method = "quantile",
        plot = FALSE)$E)
        title <- "Log2 Expression (Voom)"
        
    } else if (algorithm == "tmm") {
        data.dge <- edgeR::DGEList(counts = dataset)
        data.dge <- edgeR::calcNormFactors(data.dge)
        data.out <- as.data.frame(edgeR::cpm(data.dge, normalized.lib.sizes = TRUE,
        log = TRUE))
        title <- "Log2 Expression (TMM / edgeR)"
        
    } else if (algorithm == "deseq") {
        data.colData <- data.frame(row.names = colnames(dataset))
        data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
        design = ~1)
        data.dds <- DESeq2::estimateSizeFactors(data.dds)
        data.out <- as.data.frame(DESeq2::counts(data.dds, normalized = TRUE))
        data.out <- sign(data.out) * log2(1 + abs(data.out))
        title <- "Log2 Expression (DEseq2)"
        
    } else if (algorithm == "vsd") {
        data.colData <- data.frame(row.names = colnames(dataset))
        data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
        design = ~1)
        data.out <- SummarizedExperiment::assay(
        DESeq2::varianceStabilizingTransformation(data.dds,
        blind = TRUE))
        title <- "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
        
    }
    # Calculate percent zeroes
    perc.zero <- round(length(which(data.out == 0)) / length(as.matrix(data.out)), 3)
    
    # Build output file name:
    output.file <- paste0(opt$output, "/", algorithm, "_", opt$label, ".tsv")
    
    # Print out summary:
    cat("\n normalization method:", title,
    "\n number of genes:", nrow(data.out),
    "\n number of cells:", ncol(data.out),
    "\n percent zeroes:", perc.zero,
    "\n\n results file:", output.file)
    
    # Save normalized data to a tsv file
    data.out <- data.frame("genes" = genes, data.out)
    
    # Save normalized data to a tsv file
    readr::write_tsv(data.out, output.file)
}
