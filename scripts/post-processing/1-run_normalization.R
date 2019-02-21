# C.Savonen, CCDL for ALSF
# 2019
#
# Purpose: Normalizing count data from single rna-seq
#
# Adapted from Automated Single-cell Analysis Pipeline (ASAP) github url:
# https://github.com/DeplanckeLab/ASAP/blob/88ca3716d09197f7d178761015ea4f6a2783dd5e/R_Python/normalization.R
#
# Options:
# '-d' : Path to gene matrix in tab delimited format, gene x sample with gene
#        info as the first column. With gene info column labeled as 'gene' or
#        'genes' (not case sensitive)",
# '-a' : Normalization method to use. To run all methods use "all", to run more
#        than one method, separate the names by a space. eg -a tmm log
#        If scran is one of the methods chosen, it will be ran first because 
#        samples with negative factor sizes will need to be removed. To make it
#        parallel with the other methods, those same samples will be removed 
#        from the dataset before other normalization methods are run. 
#        Options:
#           'scale'- scale and center each sample's data using base::scale
#           'log'- just log transform alone
#           'voom'- quantile normalization from the DEseq2::voom function
#           'deseq2'- Use DEseq2 size factors
#           'vsd'- Use DEseq2 varianceStabilizingTransformation function
#           'tmm'- trimmed mean M values by edgeR
#           'scran' - pooled CPM values
#           'all'- run all the above methods
# '-o' : Directory where you would like the output to go. Default is current
#        directory. New folder will be created if the specified folder doesn't
#        exist.
# '-l' : Optional label for output files"
# '-n' : Option to override scran negative factor warning, remove the problematic
#        genes and normalize. Default is to not continue with scran normalization 
#        if 100 or more genes need to be removed.  
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
              with gene info column labeled as 'gene' or 'genes' (not case sensitive)",
              metavar = "character"),
  make_option(opt_str = c("-a", "--algorithm"), type = "character",
              default = "none", help = "Normalization method to use. Options:
              'scale', log', 'voom', 'deseq2', 'vsd', , 'tmm', 'scran' or use 
              'all' to run all algorithms. If multiple methods are wanted, 
              separate with a space eg: 'log voom'",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character",
              default = getwd(), help = "Directory where you would like the
              output to go", metavar = "character"),
  make_option(opt_str = c("-l", "--label"), type = "character",
              default = "", help = "Optional label for output files",
              metavar = "character"),
  make_option(opt_str = c("-n", "--negative"), action = "store_true",
              default = FALSE, help = "Option to override scran negative factor 
              warning, and remove the problematic samples from dataset anyway
              and then normalize")
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

# Append an underscore to the label if it is given
if (opt$label != "") {
  opt$label <- paste0("_", opt$label)
}

#---------------------------Set up algorithms to run---------------------------#
# List all the supported algorithms
all.algorithms <- c('scran', 'scale', 'log', 'voom','tmm', 'deseq2', 'vsd')

# If "all" is chosen, put all the normalization methods in list
if (opt$algorithm == "all") {
    opt$algorithm <- all.algorithms
} else {
    # Separate if multiple algorithms have been requested
    opt$algorithm <- unlist(strsplit(opt$algorithm, " "))
    
    # If scran is in the list, make sure it is the first one to be run
    if ('scran' %in% opt$algorithm) {
      opt$algorithm <- c("scran", opt$algorithm[which(opt$algorithm != "scran")])
    }
}

# Check that the algorithms chosen are in the supported list
if (any(is.na(match(opt$algorithm, all.algorithms)))) {
    stop("No supported normalization method was given. Check for typos.
         Use option -a to choose and algorithm. Acceptable options:'scale', 'log',
         'voom','tmm','deseq2', 'vsd', 'scran', or 'all")
}
#----------------------------------Load data-----------------------------------#
# Read in a tsv file of data
dataset <- readr::read_tsv(opt$data, guess_max = 10000)

# Find which column has gene info
gene.col <- grep("gene", colnames(dataset), ignore.case = TRUE)

# Create the output for results folder if it does not exist
if (length(gene.col) < 1) {
  stop("Error: Cannot find a column with gene info. (Looking for 'gene' as
       column name")
}

# Separate genes from the numeric data
genes <- dataset[, gene.col]

# Make dataset the data only
dataset <- dataset[, -gene.col]

# Turn NAs into zeroes if NAs exist
if (any(is.na(dataset))) {
  num.nas <- length(is.na(dataset))
  message(paste("Warning: there are", num.nas, 
                "NA's in your dataset. Turning them into zeroes."))
  dataset <- dplyr::mutate_all(dataset, zoo::na.aggregate)
}

#-----------------------Run each algorithm in the list-------------------------#
# For each algorithm, run through this loop
for (algorithm in opt$algorithm) {
  
  # Print out progress message:
  cat("\n \n Running:", algorithm, "...")
  
  # Run normalization algorithms
  if (algorithm == "scran") {
    
    # scater wants the data to be rounded
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = round(as.matrix(dataset))))
    
    # Make some clusters
    sce <- suppressWarnings(scran::computeSumFactors(sce))
    neg.fact <- which(sce@int_colData@listData$size_factor < 0)
    
    # Dealing with negative size factors:
    if (length(neg.fact) < 100) {
      message(paste(length(neg.fact), "negative size factors from
                    scran::computeSizeFactors calculations were found.
                    These samples will be removed from normalization output."))
      
      # Remove these samples from final output as well as the original dataset
      sce <- sce[, -neg.fact]
      dataset <- dataset[, -neg.fact]
      
    } else {
      if (!opt$negative) {
        warning(paste("Stopping scran normalization. There are",
                      length(neg.fact),
                      "negative size factors from scran::computeSizeFactors 
                       calculations. Arbitrary cutoff is set at 100. Use -n to 
                       override, remove these samples, and normalize anyway."))
      } else {
        message(paste("Option -n is being used. ", length(neg.fact),
                      " samples are being removed because they have negative size 
                       factors due to heavy zero inflation."))
        
        # Remove these samples from final output as well as the original dataset
        sce <- sce[, -neg.fact]
        dataset <- dataset[, -neg.fact]
      }
    }
    if (length(neg.fact) > 0) {
      # Write copy of counts to file 
      readr::write_tsv(dataset, file.path(opt$output, "matching_counts.tsv"))
    }
    
    # Normalize the data
    sce <- scater::normalize(sce)
    
    # Convert to edgeR so we can extract the data as a matrix more easily
    data.dge <- scran::convertTo(sce, type = "edgeR")
    data.out <- as.data.frame(edgeR::cpm(data.dge, normalized.lib.sizes = TRUE,
                                         log = TRUE))
    title <- "Cluster 'scran' Normalized CPMs"
    
  } else if (algorithm  == "scale") {
    data.out <- as.data.frame(scale(dataset))
    title <- "Scaled Expression"
    
  } else  if (algorithm == "log") {
    data.out <- sign(dataset) * log2(1 + abs(dataset))
    title <- "Log2 Expression"
    
  } else  if (algorithm == "voom") {
    # Quantile normalization using limma
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
                                               design = ~ 1)
    data.dds <- DESeq2::estimateSizeFactors(data.dds)
    data.out <- as.data.frame(DESeq2::counts(data.dds, normalized = TRUE)) 
    data.out <- sign(data.out) * log2(1 + abs(data.out)) %>% d
    title <- "Log2 Expression (DEseq2)"
    
  } else if (algorithm == "vsd") {
    data.colData <- data.frame(row.names = colnames(dataset))
    data.dds <- DESeq2::DESeqDataSetFromMatrix(dataset, colData = data.colData,
                                               design = ~ 1)
    
    # This takes more time than the other methods:
    data.out <- SummarizedExperiment::assay(
      DESeq2::varianceStabilizingTransformation(data.dds, blind = TRUE))
    title <- "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
  }
  
  # Calculate percent zeroes
  perc.zero <- round(length(which(data.out == 0)) / length(as.matrix(data.out)), 3)
  
  # Build output file name:
  output.file <- file.path(opt$output, paste0(algorithm, opt$label, ".tsv"))
  
  # Print out summary:
  cat("\n normalization method:", title,
      "\n number of genes:", nrow(data.out),
      "\n number of cells:", ncol(data.out),
      "\n percent zeroes:", perc.zero,
      "\n\n results file:", output.file)
  
  # Tack on the gene column
  data.out <- data.frame(genes, data.out)
  
  # Save normalized data to a tsv file
  readr::write_tsv(data.out, output.file)
}

