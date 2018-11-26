# C. Savonen 
# CCDL for ALSF 
#
# Purpose: Make gene matrix and do quality testing 

# Attach necessary libraries
library(Rsubread)
library(org.Hs.eg.db)
library(optparse)

#------------------------Get options using optparse----------------------------#
option_list <- list( 
  make_option(opt_str = c("-d", "--dir"), type = "character", default = NULL, 
              help = "Directory of the fastqc reports",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character", 
              default = getwd(), 
              help = "Directory where results should be placed",
              metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))


# Set directory
setwd(file.path("data", "aligned_reads"))
bam.files <- grep("\\.bam$", dir(), value = TRUE)

#-------------------------Get counts for HISAT Data----------------------------#
# Get the counts of all the gene features for each sample
# If you are using a non-common genome, can use ext.ann argument to give your own 
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

# Get the ensembl gene ids that correspond to the entrez gene ids
hisat.data <- data.frame('ensembl' = mapIds(org.Hs.eg.db, keys = rownames(hisat.data), 
                                            column = "ENSEMBL", keytype = "ENTREZID"),
                         'gene' = mapIds(org.Hs.eg.db, keys = rownames(hisat.data), 
                                            column = "SYMBOL", keytype = "ENTREZID"),
                         hisat.data, stringsAsFactors = FALSE)

# Save HISAT to an RDS file
saveRDS(hisat.data, file = file.path("..", "..", "results", "hisat.data.RDS"))

#--------------------------Proportion of mapped reads--------------------------#
# Get the proportion of mapped reads
hisat.prop.assigned <- vapply(hisat, function(x){
                        x$stat[1,2]/sum(x$stat[1:2,2])},
                        FUN.VALUE = 1)

# Make a histogram of this information
png("proportion_reads_assigned_hist.png")
hist(hisat.prop.assigned, main = "Proportion of Assigned Reads", breaks = 20)
dev.off()

setwd("../")
#------------------------------Import Salmon reads----------------------------#
# Check out the salmon files
setwd(file.path("data", "salmon_quants"))

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
salmon.annot <- gsub("\\.[0-9]+$", "", unlist(salmon.annot))

# Make into a make annotation into a dataframe
salmon.data <- data.frame('ensembl' = salmon.annot,
                          'gene' = mapIds(org.Hs.eg.db, keys = salmon.annot, 
                                          column = "SYMBOL", keytype = "ENSEMBL"),
                           salmon.data,
                           stringsAsFactors = FALSE)

# Save to an RDS file
saveRDS(salmon.data, file = file.path("..", "..", "results", "salmon.data.RDS"))
