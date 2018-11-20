# C. Savonen 
# CCDL ALSF 2018 
# Purpose: Compare author processed data, Salmon processed, and HISAT2 processed
# 
#-------------------------------Import Data------------------------------------#
library(GEOquery)
geo.meta <- getGEO("GSE84465", destdir = "data")
gunzip("data/GSE84465_series_matrix.txt.gz")

if (!file.exists("SRA.to.GSM.RDS")){
# Get geo metadata
samples <- geo.meta[[1]]@phenoData@data$geo_accession

# Get the damn SRA id from the webpage
sra <- lapply(samples, function(x) {
    sample.page <- readLines(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", x))
    sample.sra <- grep("SRA", sample.page)
    sra.id <- unlist(sample.page[sample.sra[2]+1])
    start.id <- regexpr("term=", sra.id)[1] + 5
    stop.id <- regexpr("\">", sra.id)[1] - 1
    return(substr(sra.id, start.id, stop.id))
    })

sra.gsm <- as.list(samples)
names(sra.gsm) <- unlist(sra)
saveRDS(sra.gsm, file = "SRA.to.GSM.RDS")

} else {
sra.gsm <- readRDS("SRA.to.GSM.RDS")
}

# Get author normalized data
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84465&format=file&file=GSE84465%5FGBM%5FAll%5Fdata%2Ecsv%2Egz", destfile = "data/GSE84465.csv.gz")
gunzip("data/GSE84465.csv.gz")

# Import all the different datasets
author.data <- read.csv("data/GSE84465.csv", stringsAsFactors = FALSE)
colnames(author.data) <- samples

# Load in file that contains both SRA and SRX IDs
sra.srx <- read.csv("SRA.files.csv")

# Load in Salmon quantification
salmon.data <- readRDS(file.path("results", "salmon.data.RDS"))

# Obtain GSM ids
gsm.ids <- dplyr::recode(sra.srx$experiment[match(
  colnames(salmon.data)[-1], sra.srx$run)], !!!sra.gsm)

# Make these the column names 
colnames(salmon.data)[-1] <- as.character(gsm.ids)

# Load in hisat data
hisat.data <- readRDS(file.path("results", "hisat.data.RDS"))

# Obtain GSM ids
gsm.ids <- dplyr::recode(sra.srx$experiment[match(
  colnames(hisat.data)[-1], sra.srx$run)], !!!sra.gsm)

# Make these the colnames
colnames(hisat.data)[-1] <- as.character(gsm.ids)

#------------------------Compare hisat and author data-------------------------#

compareData <- function (data = data, org.data = author.data, label = "label") {

# Get genes 
genes <- data$gene

# Only look at the data columns that are samples
data <- data[, grep("GSM", colnames(data))]

# Get indices of matching samples and genes
sample.ind <- match(colnames(data),
                    colnames(org.data))
gene.ind <- match(genes,
                  rownames(org.data))

# Combine the data into one dataset
combined <- cbind(data[!is.na(gene.ind),], 
org.data[gene.ind[!is.na(gene.ind)], sample.ind])

process <- as.factor(c(rep(label, ncol(data)), rep("original", ncol(org.data))))
pca <- prcomp(t(combined))

# This function will plot the distribution of the PC scores looks like and label samples by the different variables' different factor levels
pc.plot <- function(dat,var){
  colz <- colors(distinct=TRUE)[runif(length(levels(var)),min=1,max=length(colors(distinct=TRUE)))]
  plot(dat,pch=21,bg=colz[var]);
  legend(x="topleft", legend = levels(var),fill=colz,cex=.8)
}

# Print out a PCA plot
jpeg(file.path("..", "results", paste0("PCA", label, ".jpeg")))
pc.plot(pca$x, process)
dev.off()

Rsqs <- c()
  for (sample in 1:length(sample.ind)) {
    # Correlate the two datasets
    Rsq <- cor(data[!is.na(gene.ind),][,sample],
               
               org.data[gene.ind[!is.na(gene.ind)], sample.ind][,sample])

    # Keep all the Rsq's
    Rsqs <- c(Rsqs, Rsq)

    # Save scatterplot to jpeg
    jpeg(file.path("..", "results", paste0("Sample", sample, label, ".jpeg")))
    plot(data[!is.na(gene.ind),][,sample],
     org.data[gene.ind[!is.na(gene.ind)], sample.ind][,sample],
     xlab = paste0("HISAT"),
     ylab = paste0("Original"),
     main = paste0("R = ", round(Rsq, 3)))
    dev.off()
    }

}
# Run this with the hisat data
compareData(hisat.data, author.data, "HISAT")

# Run this with the salmon data
compareData(salmon.data, author.data, "Salmon")

