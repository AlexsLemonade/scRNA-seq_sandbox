# C.Savonen, CCDL for ALSF
# 2019
# Adapted from ASAP github 

# Purpose: Using Seurat for post-processing of scRNA-seq gene matrix .tsv 
#          In gene x sample format. 

# Options:
# "-d" - Directory of where dataset(s) to analyze exist
# "-m" - file of metadata information to analyze. Each column in this metadata 
#        will be analyzed. Input should be a path to a tsv file 
# 
# Command line example:
#
# Rscript scripts/4a-seurat_normalize.R \
# -d data/salmon_quants \
# -o data \
# -l "patel"
# 
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Import functions for this analysis
source(file.path("scripts", "util", "clustering_statistics_functions.R"))

# ggplot2 library
library(optparse)
require(jsonlite)

# These are the necessary packages
packages <- c("optparse", "scLVM", "DESeq2", "limma", "edgeR")

# Check if these packages are installed and install them if they aren't
lapply(packages, function(package) {
  if (!(package %in% installed.packages())) {
    install.packages(package)
  }
})
#--------------------------------Set up options--------------------------------#
option_list <- list(
  make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
              help = "Path to gene matrix in tab delimited format, gene x sample
              with gene info as the first column",
              metavar = "character"),
  make_option(opt_str = c("-b", "--batch"), type = "character",
              default = "none", help = "Path to metadata file that contains 
              only the batch info that you wish to normalize by.",
              metavar = "character"),
  make_option(opt_str = c("-a", "--algorithm"), type = "character",
              default = "none", help = "Normalization method to use. Options: 
              'log', 'voom', 'DESeq2', 'TMM', 'scVLM' ", metavar = "character"),
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
opt$data <- args[1]
opt$output <- args[2]
opt$algorithm <- args[3]
opt$batch <- args[4]

### Load data
data.parsed <- readr::read_tsv(opt$data)

### Run Normalization algorithms
if (opt$algorithm == "scale"){ # default []
  data.out <- as.data.frame(scale(data.parsed))
  title <- "Scaled Expression"  
  
} else  if (opt$algorithm == "log") { # default []
  data.out <- sign(data.parsed) * log2(1 + abs(data.parsed)) # Take into account normalized data with negative values
  title <- "Log2 Expression"
  
} else  if (opt$algorithm == "voom") { # Default []
  require(limma)
  data.out <- as.data.frame(voom(counts = data.parsed, 
                                 normalize.method = "quantile",
                                 plot = FALSE)$E)
  title <- "Log2 Expression (Voom)"
} else if (opt$algorithm == "tmm") { # default []
  require(edgeR)
  data.dge <- DGEList(counts = data.parsed)
  data.dge <- calcNormFactors(data.dge)
  data.out <- as.data.frame(cpm(data.dge, normalized.lib.sizes=TRUE, log = TRUE))
  title <- "Log2 Expression (TMM / edgeR)"
  
} else if (opt$algorithm == "deseq") { # default []
  require(DESeq2)
  data.colData <- data.frame(row.names = colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
  data.dds <- estimateSizeFactors(data.dds)
  data.out <- as.data.frame(counts(data.dds, normalized=TRUE))
  data.out <- sign(data.out) * log2(1 + abs(data.out))
  title <- "Log2 Expression (DEseq2)"
  
} else if (opt$algorithm == "vsd") { # default []
  require(DESeq2)
  data.colData <- data.frame(row.names = colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
  data.out <- assay(varianceStabilizingTransformation(data.dds, blind = TRUE))
  title <- "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
  
} else if (opt$algorithm == "rld") { 
  require(DESeq2)
  data.colData <- data.frame(row.names = colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
  data.out <- assay(rlogTransformation(data.dds, blind = TRUE))
  title <- "Log2 Expression (rlogTransformation - DEseq2)"
}


# Print out summary: 
cat("normalization method:", title, 
    "\n number of genes:", nrow(data.out), 
    "\n number of cells:", ncol(data.out),
    "\n number of zeroes:", length(which(data.out == 0)),
    "\n warnings:", data.warnings)
    
# Save normalized data to a tsv file
data.out$Genes = rownames(data.out)
readr::write_tsv(data.out[,c("Genes", data.out.cols)], 
                 paste0(opt$output, "/", opt$algorithm, "_", opt$label, ".tsv"),
                 sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)