# C. Savonen 
# CCDL ALSF 2018 
# Purpose: Compare author processed data, Salmon processed, and HISAT2 processed
# 
# Magrittr pipe
`%>%` <- dplyr::`%>%`
setwd("../..")

#--------------------------- Create ID conversion key--------------------------#
if (!file.exists(file.path("data", "sample_id_key.RDS"))) {
  # Get geo metadata
  geo.meta <- GEOquery::getGEO("GSE84465", destdir = "data")
  
  # Get the GSM and SRX ids from the GEO metadata
  id.key <- data.frame(gsm.ids = geo.meta[[1]]@phenoData@data$geo_accession,
                     plate.ids = geo.meta[[1]]@phenoData@data$description.1,
                     srx.ids = stringr::word(geo.meta[[1]]@phenoData@data$relation.1,
                                             start = 2, sep = "term="))
  # Import SRA file names
  sra.files <- read.csv(file.path("data", "SRA.files.csv"))[, -1]

  # Merge both keys into a single id.key dataframe
  id.key <- merge(id.key, sra.files, by.x = "srx.ids", by.y = "experiment")

  # Save this dataframe for later
  saveRDS(id.key, file = file.path("data", "sample_id_key.RDS"))
  
  # Save the rest of the metadata as a csv
  geo.meta <- data.frame(geo.meta[[1]]@phenoData@data)
  write.csv(geo.meta, file = file.path("data", "meta_data.csv"))
  
  rm(geo.meta)
  } else {
  id.key <- readRDS(file.path("data", "sample_id_key.RDS"))
  }
#----------------------------- Set up author data------------------------------#

if (!file.exists(file.path("data", "author.data.RDS"))) {
  # Get author normalized data if not downloaded yet
  if (!(file.exists(file.path("data", "GSE84465.csv")))) {
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84465&format=file&file=GSE84465%5FGBM%5FAll%5Fdata%2Ecsv%2Egz",
              destfile = file.path("data", "GSE84465.csv.gz"))
        gunzip(file.path("data", "GSE84465.csv.gz"))
  }

  # Import all the different datasets
  author.data <- read.table(file.path("data", "GSE84465.csv"), sep = " ",
                          stringsAsFactors = FALSE)

  # Make plate id to GSM id key
  convert.key <- as.list(id.key$gsm.ids)
  names(convert.key) <- as.character(id.key$plate.ids)

  # Obtain GSM ids using conversion key and make these the column names 
  colnames(author.data) <- dplyr::recode(gsub("^X", "", colnames(author.data)),
                                       !!!convert.key)

  # Add genes as their own column
  author.data <- author.data %>% tibble::rownames_to_column("gene")

  # Save cleaned data
  saveRDS(author.data, file = file.path("data", "author.data.RDS"))
  } else {
  author.data <- readRDS(file.path("data", "author.data.RDS"))
  }

#----------------------Change sample column names to GSM-----------------------#
# Make our SRA conversion key 
convert.key <- as.list(as.character(id.key$gsm.ids))
names(convert.key) <- as.character(id.key$run)

# Import datasets from their RDS files
salmon.data <- readRDS(file.path("data", "salmon.data.RDS"))
hisat.data <- readRDS(file.path("data", "hisat.data.RDS"))

# Obtain GSM ids using conversion key and make these the column names 
colnames(salmon.data) <- dplyr::recode(colnames(salmon.data), !!!convert.key)
colnames(hisat.data) <- dplyr::recode(colnames(hisat.data), !!!convert.key)

# Save objects with the GEO accession IDs
saveRDS(salmon.data, file.path("data", "salmon.data.RDS"))
saveRDS(hisat.data, file.path("data", "hisat.data.RDS"))

#------------------------------Compare datasets--------------------------------#
#
# This function will plot the distribution of the PC scores looks like and label samples by the different variables' different factor levels
pc.plot <- function(dat, var){
  colz <- colors(distinct = TRUE)[runif(length(levels(var)), min = 1,
                                        max = length(colors(distinct=TRUE)))]
  plot(dat,pch = 21, bg = colz[var]);
  legend(x = "topleft", legend = levels(var), fill = colz, cex = 0.8)
}

compareData <- function (data1 = data1, data2 = data2, label1 = "1", label2 = "2") {
  # Uses correlation and PCA to compare two datasets
  #
  # Args:
  #   data1: a dataset to compare to data2, must have colnames that are "GSM" 
  #          and a column labeled "gene" that contains gene ids to match to data2
  #   data2: a dataset to compare to data1, must have colnames that are "GSM" 
  #          and a column labeled "gene" that contains gene ids to match to data2
  #   label1: character string with the label you would like the output to have for data1
  #   label2: character string with the label you would like the output to have for data2
  #
  # Returns:
  #   PCA of data combined
  #   Gene quantification correlations of samples that are in both datasets 

  # Only keep columns that match the datasets  
  data2 <- data2 %>% dplyr::select(dplyr::matches(paste0(colnames(data1), collapse = "|")))  
  data1 <- data1 %>% dplyr::select(dplyr::matches(paste0(colnames(data2), collapse = "|")))
  
  # Join datasets into one
  combined <- dplyr::inner_join(data1, data2, by = "gene")

  # Get a vector of sample names 
  samples <- gsub("\\.x$",  "", grep("\\.x$", colnames(combined), value = TRUE))

  # Add respective column names
  colnames(combined) <- gsub("\\.x$", label1, colnames(combined))
  colnames(combined) <- gsub("\\.y$", label2, colnames(combined))

  # Make a vector with the dataset labels
  datasets <- rep(label2, ncol(combined) - 1)
  datasets[grepl(label1, colnames(combined))[-1]] <- label1
  datasets <- factor(datasets)

  # Make labels for which data is processed which way
  pca <- prcomp(t(combined[, -1]))

  # Print out a PCA plot
  png(file.path("results", paste0("PCA_", label1, "_vs_", label2, ".png")))
  pc.plot(pca$x, datasets)
  dev.off()
  
  corr.file <- file.path("results", paste0(label1, "_vs_", label2, "_sample_correlations"))
  # Create separate folder for sample correlations 
  if (!dir.exists(corr.file)) {
    dir.create(corr.file)
  }
  
  sample.corrs <- vapply(samples, function(sample) {
    # Get sample data
    sample.data <- combined %>% dplyr::select(dplyr::contains(sample))
    
    # Correlate the two datasets
    sample.corr <- cor(sample.data[ , 1], sample.data[,2])
    
    # Save scatterplot to jpeg
    png(file.path(corr.file, paste0(sample, "_", label1,"_vs_", label2, ".png")))
    plot(sample.data,
        xlab = paste0(label1),
        ylab = paste0(label2),
        main = paste0("R = ", round(sample.corr, 3)))
    dev.off()
    return(sample.corr)  
    }, FUN.VALUE = 1)

  # Plot sample correlations on a histogram so we can get an overall picture
  png(file.path( "results", paste0("hist_sample_corrs_", label1,"_vs_", label2, ".png")))
  hist(sample.corrs, xlab = "", main = paste(label1, "vs", label2, "Sample Correlations"), breaks = 20)
  dev.off()
}

# Run this with the hisat data
compareData(data1 = hisat.data, data2 = author.data, label1 = "HISAT", label2 = "author_data")

# Run this with the salmon data
compareData(data1 = salmon.data, data2 =  author.data, label1 = "Salmon", label2 = "author_data")
