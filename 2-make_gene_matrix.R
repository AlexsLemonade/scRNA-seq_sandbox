# C. Savonen 
# CCDL for ALSF 
#
# Purpose: Make gene matrix and do quality testing 
# Compare Salmon vs traditional mapping alignment 

# Use Subread's R package to analyze the bam files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsubread", version = "3.8")

library(Rsubread)

# Set directory
setwd(file.path("data", "aligned_reads"))
bam.files <- grep("\\.bam$", dir(), value = TRUE)

# Get the counts of all the gene features for each sample
# If you are using a non common genome, can use ext.ann argument to give your own 
# GTF or other file. 

hisat <- lapply(bam.files, function(x) {
  featureCounts(x, annot.inbuilt = 'hg38', isPairedEnd = TRUE) 
  })

# Make a matrix of the data
hisat.data <- do.call("cbind", lapply(hisat, function(x){
                                x$counts}))
# Carry over ID names
dimnames(hisat.data)[[2]] <- gsub("\\.bam", "", bam.files)

# Make an annotation matrix
hisat.annot <- hisat[[1]]$annotation

# Get the proportion of mapped reads
hisat.prop.assigned <- vapply(hisat, function(x){
                        x$stat[1,2]/sum(x$stat[1:2,2])},
                        FUN.VALUE = 1)

#------------------------------Import Salmon reads----------------------------#
# Check out the salmon files
setwd(file.path("..", "salmon_quants"))

# Get the names of all the folders
salmon.folders <- dir() 

# Read the data into a list
salmon <- lapply(salmon.folders, function(x) read.table(file.path(x, "quant.sf"),
                                                      header = TRUE))

# Make a matrix of the data
salmon.data <- do.call("cbind", lapply(salmon, function(x){ x$NumReads }))

# Carry over ID names
dimnames(salmon.data)[[2]] <- gsub("_quant", "", salmon.folders)

# Get Salmon ensembl gene annotation IDs
salmon.annot <- strsplit(as.character(salmon[[1]]$Name), "\\|")
salmon.annot <- lapply(salmon.annot, function(x) grep("ENSG", x, value = TRUE))

# Make into a make annotation into a dataframe
salmon.annot <- data.frame('ensembl' = unlist(salmon.annot), salmon[[1]][,2:3],
                           stringsAsFactors = FALSE)

#-------------------------Compare Salmon with HISAT2---------------------------#
library(org.Hs.eg.db)

# Make sure samples are in the same order:
all(match(dimnames(salmon.data)[[2]], dimnames(hisat.data)[[2]]) == 1:ncol(hisat.data))

hisat.data <- data.frame('ensembl' = mapIds(org.Hs.eg.db, keys = rownames(hisat.data), 
                                                column = "ENSEMBL", keytype = "ENTREZID"),
                          hisat.data)


annot.df <- transform(merge(symbol.df, df, by = "row.names"),
                      row.names = Row.names, Row.names = NULL)






