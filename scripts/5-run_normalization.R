# C.Savonen, CCDL for ALSF
# 2019
# Adapted from ASAP github 

# Purpose: Normalizing count data from single rna-seq

# Options:
# '-d' : Path to gene matrix in tab delimited format, gene x sample with gene 
#        info as the first column
# '-a' : Normalization method to use. 
#        Options: 
#           'scale'-
#           'log'-
#           'voom'-
#           'deseq2'-
#           'vsd'-
#           'rld'-
#           'tmm'-
# '-o' : Directory where you would like the output to go
# '-l' : Optional label for output files"
# 
# 
# Command line example:
#
# Rscript scripts/normalization_methods.R \
# -d data/darmanis_counts.tsv \
# -a log \
# -o normalized_data \
# -l "darmanis"
# 
# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
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

### Default Parameters
opt$data <- "darmanis_data/normalized_darmanis/darmanis_counts.tab"
opt$output <- "darmanis_data/normalized_darmanis"
opt$algorithm <- "vsd"
opt$label <- "darmanis"

# Check that the dimension reduction option given is supported
if (!(opt$algorithm %in% c( 'scale', 'log', 'voom', 'deseq2', 'vsd', 'rld', 'tmm'))){
  stop("That is not a normalization method supported by this script. 
       Check for typos. Acceptable options:'scale', 'log', 'voom', 'deseq2', 
       'vsd', 'rld' or 'tmm'")
}

# Create the output for results folder if it does not exist
if (!dir.exists(opt$output)) {
  message(paste("Output folder:", opt$output, "does not exist, creating one"))
  dir.create(opt$output)
}

#----------------------------------Load data-----------------------------------#
# Read in a tsv file of data
dataset <- readr::read_tsv(opt$data)

# Separate genes from the numeric data
genes <- dataset$gene
dataset <- dataset[, -1]

# Run normalization algorithms
if (opt$algorithm == "scale"){ # default []
  data.out <- as.data.frame(scale(dataset))
  title <- "Scaled Expression"  
  
} else  if (opt$algorithm == "log") { # default []
  data.out <- sign(dataset) * log2(1 + abs(dataset)) 
  title <- "Log2 Expression"
  
} else  if (opt$algorithm == "voom") { 
  data.out <- as.data.frame(limma::voom(counts = dataset, 
                                   normalize.method = "quantile",
                                   plot = FALSE)$E)
  title <- "Log2 Expression (Voom)"
  
} else if (opt$algorithm == "tmm") { 
  data.dge <- edgeR::DGEList(counts = dataset)
  data.dge <- edgeR::calcNormFactors(data.dge)
  data.out <- as.data.frame(edgeR::cpm(data.dge, normalized.lib.sizes = TRUE,
                                       log = TRUE))
  title <- "Log2 Expression (TMM / edgeR)"
  
} else if (opt$algorithm == "deseq") { 
  data.colData <- data.frame(row.names = colnames(dataset))
  data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
                                             design = ~1)
  data.dds <- DESeq2::estimateSizeFactors(data.dds)
  data.out <- as.data.frame(DESeq2::counts(data.dds, normalized = TRUE))
  data.out <- sign(data.out) * log2(1 + abs(data.out))
  title <- "Log2 Expression (DEseq2)"
  
} else if (opt$algorithm == "vsd") { 
  data.colData <- data.frame(row.names = colnames(dataset))
  data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
                                             design = ~1)
  data.out <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(data.dds, 
                                                              blind = TRUE))
  title <- "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
  
} else if (opt$algorithm == "rld") { 
  data.colData <- data.frame(row.names = colnames(dataset))
  data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
                                             design = ~1)
  data.out <- SummarizedExperiment::assay(DESeq2::rlogTransformation(data.dds, blind = TRUE))
  title <- "Log2 Expression (rlogTransformation - DEseq2)"
}

# Calculate percent zeroes
perc.zero <- round(length(which(data.out == 0)) / length(as.matrix(data.out)), 3)

# Print out summary: 
cat("normalization method:", title, 
    "\n number of genes:", nrow(data.out), 
    "\n number of cells:", ncol(data.out),
    "\n percent zeroes:", perc.zero)
    
# Save normalized data to a tsv file
data.out$genes <- genes

# Save normalized data to a tsv file
readr::write_tsv(data.out, 
                 paste0(opt$output, "/", opt$algorithm, "_", opt$label, ".tsv"))
