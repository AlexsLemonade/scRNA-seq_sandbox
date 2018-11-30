# C. Savonen 
# CCDL for ALSF 
#
# Purpose: Make gene matrix and do quality testing 

# Install these libraries if they are not installed
if (!("org.Hs.eg.db" %in% installed.packages())) { 
  BiocManager::install("org.Hs.eg.db", suppressUpdates = FALSE)
  }
if (!("rjson" %in% installed.packages())) {
  install.packages("rjson", suppressUpdates = FALSE)
}

# Attach needed libraries
library(Rsubread)
library(org.Hs.eg.db)
library(optparse)

#------------------------Get options using optparse----------------------------#
#option_list <- list( 
  #make_option(opt_str = c("-d", "--dir"), type = "character", default = NULL, 
              #help = "Directory of the bam files",
              #metavar = "character"),
  #make_option(opt_str = c("-o", "--output"), type = "character", 
              #default = getwd(), 
              #help = "Directory where results should be placed",
              #metavar = "character"))

#opt <- parse_args(OptionParser(option_list = option_list))


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
hisat.data <- do.call("cbind", lapply(hisat, function(x) x$counts))

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
saveRDS(hisat.data, file = file.path("..", "hisat.data.RDS"))

#----------------------Proportion of HISAT2 mapped reads-----------------------#
# Get the proportion of mapped reads
hisat.prop.assigned <- vapply(hisat, function(x){
                        x$stat[1,2]/sum(x$stat[1:2,2])},
                        FUN.VALUE = 1)

# Make a histogram of this information
png(file.path("..", "..", "results", "hisat2_prop_reads_mapped_hist.png"))
hist(hisat.prop.assigned, xlab = "", main = "HISAT2 Proportion of Mapped Reads", breaks = 20)
dev.off()

#------------------------------Import Salmon reads-----------------------------#
# Check out the salmon files
setwd(file.path("..", "salmon_quants"))

# Get the names of all the folders
salmon.folders <- dir() 

# Read the data into a list
salmon <- lapply(salmon.folders, function(x) read.table(file.path(x, "quant.sf"),
                                                      header = TRUE))

# Get Salmon ensembl gene annotation IDs
salmon.annot <- strsplit(as.character(salmon[[1]]$Name), "\\|")
salmon.annot <- vapply(salmon.annot, function(x) {
                        x <- grep("ENST", x, value = TRUE)
                             gsub("\\.[0-9]+$", "", x)
                          }, FUN.VALUE = "character")

# Turn transcript ids to gene ids
salmon.annot <- mapIds(org.Hs.eg.db, keys = salmon.annot, column = "ENSEMBL", 
                       keytype = "ENSEMBLTRANS")

# Make a matrix of the data
salmon.data <- do.call("cbind", lapply(salmon, function(x) {
                          tapply(x$NumReads, salmon.annot, sum)
                          }))

# Carry over ID names
dimnames(salmon.data)[[2]] <- gsub("_quant", "", salmon.folders)

# Make into a make annotation into a dataframe (removes data without a gene)
salmon.data <- data.frame('ensembl' = rownames(salmon.data),
                          'gene' = mapIds(org.Hs.eg.db, keys = rownames(salmon.data), 
                                          column = "SYMBOL", keytype = "ENSEMBL"),
                           salmon.data,
                           stringsAsFactors = FALSE)

# Save to an RDS file
saveRDS(salmon.data, file = file.path("..", "salmon.data.RDS"))

#---------------------Salmon proportion of mapped reads------------------------#
# Get the proportion of mapped reads
salmon.prop.assigned <- vapply(salmon.folders, function(x) {
  rjson::fromJSON(file = file.path(x, "aux_info", "meta_info.json"))$percent_mapped/100
  }, FUN.VALUE = 1)

# Make a histogram of this information
png(file.path("..", "..", "results", "salmon_prop_reads_mapped_hist.ensembl.png"))
hist(salmon.prop.assigned, xlab = "", main = "Salmon Proportion of Mapped Reads", breaks = 20)
dev.off()
