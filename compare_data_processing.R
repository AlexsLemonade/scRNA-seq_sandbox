# C. Savonen 
# CCDL ALSF 2018 
# Purpose: Compare author processed data, Salmon processed, and HISAT2 processed
# 
library(GEOquery)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#-----------------------------Import Author GEO Data---------------------------#
# Get geo metadata
geo.meta <- getGEO("GSE84465", destdir = "data")

# Get author normalized data
if (!(file.exists(file.path("data", "GSE84465.csv")))) {
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84465&format=file&file=GSE84465%5FGBM%5FAll%5Fdata%2Ecsv%2Egz",
              destfile = file.path("data", "GSE84465.csv.gz"))
gunzip("data/GSE84465.csv.gz")
}

# Import all the different datasets
author.data <- read.table(file.path("data", "GSE84465.csv"), sep = " ",
                          stringsAsFactors = FALSE) %>% 
                          tibble::rownames_to_column("gene")

#------------------------------Create ID conversion key------------------------#
# Get the GSM and SRX ids from the GEO metadata
id.key <- data.frame(gsm.ids = geo.meta[[1]]@phenoData@data$geo_accession,
                     srx.ids = stringr::word(geo.meta[[1]]@phenoData@data$relation.1,
                                             start = 2, sep = "term="))
# Get SRR to SRX conversion key
sra.files <- read.csv(file.path("data", "SRA.files.csv"))[, -1]

# Merge both keys into a single id.key dataframe
id.key <- merge(id.key, sra.files, by.x = "srx.ids", by.y = "experiment")

# Save this dataframe for later
saveRDS(id.key, file = file.path("data", "sample_id_key.RDS"))

# Make our conversion key 
sra.gsm <- as.list(as.character(id.key$gsm.ids))
names(sra.gsm) <- as.character(id.key$run)

#----------------------------- Set up Salmon data------------------------------#
salmon.data <- readRDS(file.path("results", "salmon.data.RDS"))

# Obtain GSM ids using conversion key and make these the column names 
colnames(salmon.data) <- dplyr::recode(colnames(salmon.data), !!!sra.gsm)

#---------------------------Set up HISAT2 mapped data--------------------------#
hisat.data <- readRDS(file.path("results", "hisat.data.RDS"))

# Obtain GSM ids using conversion key and make these the column names 
colnames(hisat.data) <- dplyr::recode(colnames(hisat.data), !!!sra.gsm)

#------------------------Compare hisat and author data-------------------------#

compareData <- function (data1 = data1, data2 = data2, label = "label") {
  # Uses correlation and PCA to compare two datasets
  #
  # Args:
  #   data1: a dataset to compare to data2, must have colnames that are "GSM" 
  #          and a column labeled "gene" that contains gene ids to match to data2
  #   data2: a dataset to compare to data1, must have colnames that are "GSM" 
  #          and a column labeled "gene" that contains gene ids to match to data2
  #   label: what label you would like the output to have
  #
  # Returns:
  #   PCA of data combined
  #   Gene quantification correlations of samples that are in both datasets 

# Only keep samples that are in both datasets
xx <- data1 %>% filter(as.name(column) == colnames(data2))
xx <- data2 %>% filter(as.name(column) == colnames(data1))

# Change column names so when we merge data1 has distinct column names
colnames(data1) <- paste0(colnames(data1), label)

# Join datasets into one
combined <- inner_join(data1, data2, by = "gene")

# Make labels for which data is processed which way
pca <- prcomp(t(combined))

# This function will plot the distribution of the PC scores looks like and label samples by the different variables' different factor levels
pc.plot <- function(dat, var){
  colz <- colors(distinct = TRUE)[runif(length(levels(var)), min = 1,
                                        max = length(colors(distinct=TRUE)))]
  plot(dat,pch = 21, bg = colz[var]);
  legend(x = "topleft", legend = levels(var), fill = colz, cex = 0.8)
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
     xlab = paste0(label),
     ylab = paste0("Original"),
     main = paste0("R = ", round(Rsq, 3)))
    dev.off()
    }

}

# Run this with the hisat data
compareData(hisat.data, author.data, "HISAT")

# Run this with the salmon data
compareData(salmon.data, author.data, "Salmon")

