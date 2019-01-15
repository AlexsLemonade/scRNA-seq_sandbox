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
data.parsed <- readr::read_tsv(opt$data, sep = "\t", header = TRUE, row.names = 1,
                          colClasses = "character", check.names = FALSE,
                          stringsAsFactors = FALSE)

### Functions
error.json <- function(displayed) {
  stats <- list()
  stats$displayed_error = displayed
  write(toJSON(stats, method="C", auto_unbox=T), file = paste0(opt$output,"/output.json"), append =F)
  stop(displayed)
}

### Run Normalization algorithms
if (opt$algorithm == "scale"){ # default []
  data.out <- as.data.frame(scale(data.parsed))
  ycol <- "Scaled Expression"  
} else  if (opt$algorithm == "log"){ # default []
  data.out <- sign(data.parsed) * log2(1 + abs(data.parsed)) # Take into account normalized data with negative values
  ycol <- "Log2 Expression"
} else  if (opt$algorithm == "voom"){ # Default []
  require(limma)
  data.out <- as.data.frame(voom(counts = data.parsed, normalize.method = "quantile", plot = FALSE)$E)
  ycol <- "Log2 Expression (Voom)"
} else if (opt$algorithm == "tmm"){ # default []
  require(edgeR)
  data.dge <- DGEList(counts = data.parsed)
  data.dge <- calcNormFactors(data.dge)
  data.out <- as.data.frame(cpm(data.dge, normalized.lib.sizes=TRUE, log = TRUE))
  ycol <- "Log2 Expression (TMM / edgeR)"
} else if (opt$algorithm == "deseq"){ # default []
    require(DESeq2)
    data.colData <- data.frame(row.names=colnames(data.parsed))
    data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
    data.dds <- estimateSizeFactors(data.dds)
    data.out <- as.data.frame(counts(data.dds, normalized=TRUE))
    data.out <- sign(data.out) * log2(1 + abs(data.out))
    ycol <- "Log2 Expression (DEseq2)"
} else if (opt$algorithm == "vsd"){ # default []
  require(DESeq2)
  data.colData <- data.frame(row.names = colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
  data.out <- assay(varianceStabilizingTransformation(data.dds, blind = TRUE))
  ycol <- "Log2 Expression (varianceStabilizingTransformation - DEseq2)"
} else if (opt$algorithm == "rld"){ # default []
  require(DESeq2)
  data.colData = data.frame(row.names = colnames(data.parsed))
  data.dds <- DESeqDataSetFromMatrix(data.parsed, colData = data.colData, design = ~1)
  data.out <- assay(rlogTransformation(data.dds, blind = TRUE))
  ycol = "Log2 Expression (rlogTransformation - DEseq2)"
} else if (opt$algorithm == "scLVM"){ # default []
  require(scLVM)
  require(DESeq2)
  fit.model <- args[5]
  if(is.na(fit.model) || fit.model == "") fit.model = 'log'
  keep.most.variable.genes <- T
  if(is.na(args[6]) || args[6] == "" || args[6] == "null" || args[6] == "false" || args[6] == "0") keep.most.variable.genes = F
  ercc.file <- args[7]
  data.ercc = NULL
  if(!is.null(ercc.file) & !is.na(ercc.file) & ercc.file != "" & ercc.file != "null"){
    data.ercc = read.table(ercc.file, sep="\t", header=T, row.names=1, colClasses=c(Genes="character"), check.names=F, stringsAsFactors=F)
    if(all(sort(colnames(data.parsed)) %in% sort(colnames(data.ercc)))) {
      print("ERCC file is correct. Process with normalization.")
    } else error.json("ERCC file is NOT correct.")
  }
  if(!is.null(data.ercc)){
    print("Normalize read counts based on ERCC.")
    sfERCC <- NULL
    tryCatch({
      sfERCC <<- estimateSizeFactorsForMatrix(data.ercc[,colnames(data.parsed)])
    }, error = function(err) {
      if(grepl("every gene contains at least one zero", err$message)) error.json("scLVM error: At least one Gene should not contain 0 reads.")
      error.json(err$message)
    })
    data.ercc.sf <- t( t(data.ercc[,colnames(data.parsed)]) / sfERCC )
    data.parsed.sf <- t( t(data.parsed) / sfERCC )
    png(paste0(opt$output,"/tech.noise.fit.png"), width=1000, height=600, type="cairo")
    data.tech.noise <- fitTechnicalNoise(data.parsed.sf,nCountsERCC=data.ercc.sf, fit_type = 'counts')
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="tech.noise.fit.png", description="scLVM fit of technical noise"))
    if(keep.most.variable.genes){
      png(paste0(opt$output,"/variable.genes.png"), width=500, height=600, type="cairo")
      data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, method = "fdr", threshold = 0.1, fit_type="counts",sfEndo=estimateSizeFactorsForMatrix(data.parsed), sfERCC=sfERCC)
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
      data.out <- as.data.frame(data.parsed.sf[data.variable.genes,])
    } else {
      data.out <- as.data.frame(data.parsed.sf)
    }
    data.out <- log2(1 + data.out)
    ycol = "Log2 ERCC Normalized Counts (scLVM [Brennecke et al.] + DEseq2)"
  } else {
    print("Normalize read counts without ERCC.")
    sfCounts <- NULL
    tryCatch({
      sfCounts <<- estimateSizeFactorsForMatrix(data.parsed)
    }, error = function(err) {
      if(grepl("every gene contains at least one zero", err$message)) error.json("scLVM error: At least one Gene should not contain 0 reads.")
      error.json(err$message)
    })
    data.parsed.sf <- t( t(data.parsed) / sfCounts )
    png(paste0(opt$output,"/tech.noise.fit.png"), width=500, height=600, type="cairo")
    data.tech.noise <- fitTechnicalNoise(data.parsed.sf,use_ERCC = F, fit_type = fit.model)
    dev.off()
    data.plots <- rbind(data.plots, data.frame(name="tech.noise.fit.png", description="scLVM fit of technical noise"))
    if(keep.most.variable.genes){
      png(paste0(opt$output,"/variable.genes.png"), width=500, height=600, type="cairo")
      data.variable.genes <- getVariableGenes(data.parsed.sf, data.tech.noise$fit, fit_type = fit.model, threshold = 0.1)
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="variable.genes.png", description="scLVM fit of technical noise + variable genes"))
      data.out <- as.data.frame(data.parsed.sf[data.variable.genes,])
    } else {
      data.out <- as.data.frame(data.parsed.sf)
    }
    data.out <- log2(1 + data.out)
    if(fit.model == 'log') ycol = "Log2 Normalized Counts (scLVM [log-linear fit] + DEseq2)"
    if(fit.model == 'logvar') ycol = "Log2 Normalized Counts (scLVM [2nd order polynomial regression (loess) fit] + DEseq2)"
  }
} else error.json("This normalization method is not implemented")

png(paste0(opt$output,"/boxplot.norm.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
boxplot(data.out, outline=F, las=2, ylab=ycol)
dev.off()

#BATCH EFFECT CORRECTION
if(!is.null(opt$batch) & !is.na(opt$batch) & opt$batch != "" & opt$batch != "null"){
  require(sva)
  print("Batch effect correction requested.")
  print("Reading batch file...")
  data.batch = read.table(opt$batch, sep="\t", header=F, row.names=1, check.names=F, stringsAsFactors=F)
  if(all(sort(colnames(data.out)) %in% sort(rownames(data.batch)))) {# If all names are in the batch file (handle filtered cells)
    print("Batch file is correct. Process with ComBat.")
    data.batch = factor(data.batch[colnames(data.out),])
    # Test that no batch has only one element
    if(length(which(table(data.batch) == 1)) != 0) {
      data.warnings <- c(data.warnings, paste0("It's not possible to perform batch effect correction with batches having only one sample [", paste(names(table(data.batch)[table(data.batch) == 1]), collapse = ", "),"]. Batch effect correction was skipped."))
    } else if(length(unique(data.batch)) == 1){
      data.warnings <- c(data.warnings, paste0("It's not possible to perform batch effect correction with only one batch across all samples [", as.character(unique(data.batch)), "]. Batch effect correction was skipped."))
    } else {
      # Test variance in every batch
      for(batch in unique(data.batch)) {
        no.variance.genes = which(apply(data.out[,which(data.batch==batch)], 1, var) == 0)
        if(length(no.variance.genes) != 0){
          data.warnings <- c(data.warnings, paste0(length(no.variance.genes), " gene(s) have no variance (var = 0) in batch ", batch,". ComBat (batch effect correction method) cannot be ran with zero-variance genes. They were removed from the dataset."))
          data.out = data.out[-no.variance.genes,]
        }
      }
      # Run Combat
      png(paste0(opt$output,"/batch.correction.output.png"), width=1000, height=1000, type="cairo")
      data.combat.2 = ComBat(dat=data.out, batch=data.batch, par.prior=T, prior.plots=T)
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="batch.correction.output.png", description="Output of ComBat for batch effect correction model."))
      data.out = data.combat.2
      png(paste0(opt$output,"/boxplot.norm.after.batch.correction.png"), width=(560 + min(10000, ncol(data.out) * 6.5)), height=600, type="cairo")
      boxplot(data.out, outline=F, las=2, ylab=paste0(ycol," + batch effect correction [ComBat]"))
      dev.off()
      data.plots <- rbind(data.plots, data.frame(name="boxplot.norm.after.batch.correction.png", description="Distribution of the expression of genes in each sample after normalization & batch effect correction."))
    }
  } else error.json("The batch file is not correct.")
}

stats <- list()
stats$nber_genes = nrow(data.out)
stats$nber_cells = ncol(data.out)
stats$nber_zeros = length(which(data.out == 0))
stats$list_plots = data.plots
stats$warnings = as.list(data.warnings)
write(toJSON(stats, method="C", auto_unbox=T), file = paste0(opt$output,"/output.json"), append=F)

data.out.cols = colnames(data.out)
data.out$Genes = rownames(data.out)
write.table(data.out[,c("Genes",data.out.cols)], file=paste0(opt$output,"/output.tab") , sep="\t", row.names = F, quote=F, append = F)