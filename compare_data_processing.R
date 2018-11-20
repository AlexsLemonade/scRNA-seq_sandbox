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
# Correlate hisat and author processed data
sample.ind <- match(colnames(hisat.data)[-c(1:2)],
                    colnames(author.data))
gene.ind <- match(hisat.data$gene,
                  rownames(author.data))

Rsqs <- c()

for (i in 1:length(sample.ind)) {
# Correlate the two datasets
Rsq <- cor(hisat.data[!is.na(gene.ind),-c(1:2)][,i], author.data[gene.ind[!is.na(gene.ind)], sample.ind][,i])

# Keep all the Rsq's
Rsqs <- c(Rsqs, Rsq)

# Save scatterplot to jpeg
jpeg(paste0("../results/Sample", i , ".jpeg"))
plot(hisat.data[!is.na(gene.ind),-c(1:2)][,i],
     author.data[gene.ind[!is.na(gene.ind)], sample.ind][,i],
     xlab = paste0("HISAT"),
     ylab = paste0("Original"),
     main = paste0("R = ", round(Rsq, 3)))
dev.off()
}

# Correlate salmon and author processed data
sample.ind <- match(dimnames(salmon.data)[[2]], dimnames(author.data)[[2]])
corr(salmon.data, author.data[gene.id, sample.ind])



