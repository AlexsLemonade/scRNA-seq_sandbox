# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: Using Seurat for post-processing of scRNA-seq Salmon processed data
#
# Options:
# "-d" - Directory of where individual samples' salmon folders are located.
# "-o" - Directory of where the output gene matrix RDS file should go.(Optional)
# "-l" - Optional label to add to output files. Generally necessary if processing
#        multiple datasets in the same pipeline.
# 
# Command line example:
#
# Rscript scripts/4a-seurat_normalize.R \
# -d data/salmon_quants \
# -o data \
# -l "patel"

#-------------------------- Get necessary packages-----------------------------#
# Adapted from the Seurat pbmc tutorial found here: https://satijalab.org/seurat/pbmc3k_tutorial.html)
# Get Seurat installed if it isn't already
if (!("Seurat" %in% installed.packages())) {
    install.packages("remotes")
    remotes::install_github("UCSF-TI/fake-hdf5r")
    install.packages('Seurat')
}
# Attach library
library(Seurat)
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
                help = "Path to counts gene matrix in tsv format, gene x sample",
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

#--------------------Import Salmon/tximport gene matrix----------------------#
# Import gene expression matrix data
tx.counts <- readr::write_tsv(opt$data)

# Mitochondrial genes by their ensembl IDs
# This list is from here :http://jdblischak.github.io/singleCellSeq/analysis/qc-filter-ipsc.html
# Identify mitochondrial genes by their ensembl IDs
mtgene <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
"ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
"ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
"ENSG00000198763","ENSG00000228253", "ENSG00000198938",
"ENSG00000198840")

# Get the gene symbol names for these genes
mtgene <- rownames(tx.counts)[match(mtgene, tx.counts$ensembl)]

# Add regex symbols for exact matches
mtgene <- paste0("^", mtgene, "$")

# Set up data in Seurat format
seurat <- CreateSeuratObject(raw.data = tx.counts, min.cells = 3,
                             min.genes = 200, project = "darmanis_gbm")

# Get those mito genes
mito.genes <- grep(pattern = paste0(mtgene, collapse = "|"),
                   x = rownames(seurat@data), value = TRUE)

# Calculate the percentage of mitochonrial transcripts to total transcripts
percent.mito <- Matrix::colSums(seurat@raw.data[mito.genes, ])/
Matrix::colSums(seurat@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
seurat <- AddMetaData(object = seurat, metadata = percent.mito,
                      col.name = "percent_mito")

# Filter out cells
seurat <- FilterCells(object = seurat, subset.names = c("nGene", "percent_mito"),
                      low.thresholds = c(200, -Inf),
                      high.thresholds = c(2500, 0.05))

# Normalizing the data
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000)

# Save this seurat object to an RDS file
readr::write_tsv(seurat, file = file.path(opt$output,
                                 paste0(opt$label, "seurat_normalize.tsv"))
