# Second round of evaluating pre-processing 
# 
# tSNE
# Try ASAP cuz why not
# Try SCONE cuz why not
# Try Seraut
if (!("tsne" %in% installed.packages())) {
  install.packages("tsne")
}

# Import gene expression matrix data
salmon.data <- readRDS(file.path("data", "salmon.data.RDS"))
author.data <- readRDS(file.path("data", "author.data.RDS"))

# Batch effect testing
geo.meta <- read.csv(file.path("data", "meta_data.csv"))
geo.meta <- geo.meta[match(colnames(salmon.data),geo.meta$geo_accession),]

# Do PCA of samples
pca <- prcomp(t(salmon.data[,-c(1:2)]))

# Do tSNE of samples
tsne.dat <- tsne::tsne(t(salmon.data[,-c(1:2)]))

# Make PCA plot function that allows us 
pc.plot <- function(dat, var){
  colz <- colors(distinct = TRUE)[runif(length(levels(var)), min = 1,
                                        max = length(colors(distinct=TRUE)))]
  plot(dat, pch = 21, bg = colz[var]);
  legend(x = "topleft", legend = levels(var), fill = colz, cex = 0.8)
}

# Plate ID PCA
png(file.path( "results", paste0("Plate_ID_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$plate.id.ch1[-c(1:2)]))
dev.off()
summary(lm(pca$x[,1]~as.factor(geo.meta$plate.id.ch1[-c(1:2)])))

# Patient ID PCA
png(file.path( "results", paste0("Patient_ID_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$patient.id.ch1[-c(1:2)]))
dev.off()
summary(lm(pca$x[,1]~as.factor(geo.meta$patient.id.ch1[-c(1:2)])))

# Tissue Type PCA
png(file.path( "results", paste0("Tissue_Type_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$characteristics_ch1.3[-c(1:2)]))
dev.off()
summary(lm(pca$x[,1]~as.factor(geo.meta$characteristics_ch1.3[-c(1:2)])))

# tSNE Plots
png(file.path( "results", paste0("tSNE_Cluster_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$tsne.cluster.ch1[-c(1:2)]))
dev.off()

png(file.path( "results", paste0("Patient_ID_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$patient.id.ch1[-c(1:2)]))
dev.off()

png(file.path( "results", paste0("Author_tSNE_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$tsne.cluster.ch1[-c(1:2)]))
dev.off()

png(file.path( "results", paste0("Plate_ID_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$plate.id.ch1[-c(1:2)]))
dev.off()

png(file.path( "results", paste0("Tissue_Type_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$characteristics_ch1.3[-c(1:2)]))
dev.off()

#------------------------------Prepping data for ASAP--------------------------#
# Format the data so ASAP will like it. 
salmon.data <- apply(salmon.data[,-c(1:2)], 2, round)

# Write these data to tab delimted text files for use in ASAP
write.table(salmon.data, file.path("data", "salmon.data.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(author.data, file.path("data", "author.data.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write batch files
geo.meta <- GEOquery::getGEO("GSE84465", destdir = "data")
geo.meta <- data.frame(geo.meta[[1]]@phenoData@data)
geo.meta <- geo.meta[match(colnames(salmon.data)[-c(1:2)], geo.meta$geo_accession), ]

# Write batch info file for ASAP use
write.table(cbind(geo.meta$geo_accession, geo.meta$plate.id.ch1),
            file.path("data", "batch.info.txt"), row.names = FALSE, 
            col.names = FALSE, sep = "\t", quote = FALSE)

#----------------------------Using SCONE for analyses--------------------------#
# Installing SCONE
source("https://bioconductor.org/biocLite.R")
biocLite("scone")

#-----------------------------Using Seraut for analyses------------------------#
# Get Seraut
install.packages('Seurat')
library(Seurat)

# Make our data into a Seraut Object
author.seraut <- CreateSeuratObject(raw.data = author.data, min.cells = 3,
                                    min.genes = 200, project = "darmanis_gbm")

# Get those mito genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = author.seraut@data),
                   value = TRUE)
percent.mito <- Matrix::colSums(author.seraut@raw.data[mito.genes, ])/
  Matrix::colSums(author.seraut@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
author.seraut <- AddMetaData(object = author.seraut, metadata = percent.mito,
                             col.name = "percent.mito")
VlnPlot(object = author.seraut, features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)

# Filter out cells
author.seraut <- FilterCells(object = author.seraut, subset.names = c("nGene", "percent.mito"), 
                             low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
# Normalizing the data
# After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.

author.seraut <- NormalizeData(object = author.seraut, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

