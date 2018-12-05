# C.Savonen
# CCDL for ALSF 2018
# 
# Using Seraut for post-processing of scRNA-seq Salmon processed data

# Get Seraut installed if it isn't already
if (!('Seraut' %in% installed.packages())) {
  remotes::install_github("UCSF-TI/fake-hdf5r")
  install.packages('Seurat')
}
# Attach library
library(Seurat)

# Import gene expression matrix data
salmon.data <- readRDS(file.path("data", "salmon.data.RDS"))

# Mitochondrial genes by their ensembl IDs
# This list is from http://jdblischak.github.io/singleCellSeq/analysis/qc-filter-ipsc.html

mtgene <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888", 
            "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
            "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
            "ENSG00000198763","ENSG00000228253", "ENSG00000198938", 
            "ENSG00000198840")

# Get the gene symbol names for these genes add regex for exact matches
mtgene <- salmon.data$gene[match(mtgene, salmon.data$ensembl)]
mtgene <- paste0("^", mtgene, "$")

# Get gene factor vector to sort by
genes <- as.factor(salmon.data$gene)

# Change data to a matrix with only one gene entry per 
salmon.data <- as.list(salmon.data[,-c(1:2)])
salmon.data <- do.call("cbind", 
              lapply(salmon.data, function(sample) tapply(sample, genes, sum))
              )

# Batch effect testing
geo.meta <- read.csv(file.path("data", "meta_data.csv"))
geo.meta <- geo.meta[match(colnames(salmon.data),geo.meta$geo_accession),]

# Make our data into a Seraut Object
salmon.seraut <- CreateSeuratObject(raw.data = salmon.data, min.cells = 3,
                                    min.genes = 200, project = "darmanis_gbm")

# Get those mito genes
mito.genes <- grep(pattern = paste0(mtgene, collapse = "|"), x = rownames(salmon.seraut@data), value = TRUE)

percent.mito <- Matrix::colSums(salmon.seraut@raw.data[mito.genes, ])/
  Matrix::colSums(salmon.seraut@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
salmon.seraut <- AddMetaData(object = salmon.seraut, metadata = percent.mito,
                             col.name = "percent_mito")

# Make a violin plot
png(file.path("results", "violin_plot.png"), width = 1000, height = 500)
VlnPlot(object = salmon.seraut, features.plot = c("nGene", "nUMI", "percent_mito"),
        nCol = 3)
dev.off()

# Examine relationships of nUMI, percent_mito and nGene
png(file.path("results", "gene_nUMI_plots.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))
GenePlot(object = salmon.seraut, gene1 = "nUMI", gene2 = "percent_mito")
GenePlot(object = salmon.seraut, gene1 = "nUMI", gene2 = "nGene")
dev.off()

# Filter out cells
salmon.seraut <- FilterCells(object = salmon.seraut, 
                             subset.names = c("nGene", "percent_mito"), 
                             low.thresholds = c(200, -Inf), 
                             high.thresholds = c(2500, 0.05))
# Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data.
#  By default, we employ a global-scaling normalization method “LogNormalize” 
# that normalizes the gene expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
salmon.seraut <- NormalizeData(object = salmon.seraut, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# Make a variance to average expression plot
png(file.path("results", "dispersion_avg_exp.png"), width = 1000, height = 700)
salmon.seraut <- FindVariableGenes(object = salmon.seraut,
                                   mean.function = ExpMean,
                                   dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.0125,
                                   x.high.cutoff = 3,
                                   y.cutoff = 0.5)
dev.off()

# Save this seraut object to an RDS file
saveRDS(salmon.seraut, file = file.path("data", "salmon.seraut"))
