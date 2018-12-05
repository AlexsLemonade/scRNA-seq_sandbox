# C. Savonen 
# CCDL for ALSF 
#
# Purpose: Make gene matrix and do quality testing 
#-------------------------- Get necessary packages-----------------------------#
if (!("org.Hs.eg.db" %in% installed.packages())) { 
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db", suppressUpdates = TRUE)
}
if (!("rjson" %in% installed.packages())) {
  install.packages("rjson", suppressUpdates = FALSE)
}

# Attach needed libraries
library(org.Hs.eg.db)
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list( 
  make_option(opt_str = c("-d", "--dir"), type = "character", default = getwd(),
              help = "Directory where salmon quantification folders are located",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character", default = getwd(),
              help = "Directory where you would like the output to go",
              metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))

#------------------------------Import Salmon reads-----------------------------#
# Get the names of all the folders
salmon.folders <- dir(opt$dir) 

# Read the data into a list
salmon <- lapply(salmon.folders, function(x) read.table(file.path(opt$dir, x, "quant.sf"),
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
dimnames(salmon.data)[[2]] <- salmon.folders
e
# Make into a make annotation into a dataframe (removes data without a gene)
salmon.data <- data.frame('ensembl' = rownames(salmon.data),
                          'gene' = mapIds(org.Hs.eg.db, keys = rownames(salmon.data), 
                                          column = "SYMBOL", keytype = "ENSEMBL"),
                           salmon.data,
                           stringsAsFactors = FALSE)

# Save to an RDS file
saveRDS(salmon.data, file = file.path(opt$output, "salmon.data.RDS"))

#---------------------Salmon proportion of mapped reads------------------------#
# Get the proportion of mapped reads
salmon.prop.assigned <- vapply(salmon.folders, function(x) {
  rjson::fromJSON(file = file.path(opt$dir, x, "aux_info", "meta_info.json"))$percent_mapped/100
  }, FUN.VALUE = 1)

# Make a histogram of this information
png(file.path(opt$output, "salmon_prop_reads_mapped_hist.ensembl.png"))
hist(salmon.prop.assigned, xlab = "", main = "Salmon Proportion of Mapped Reads", breaks = 20)
dev.off()
