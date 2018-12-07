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

#--------------------------- Create ID conversion key--------------------------#
if (!file.exists(file.path("sample_id_key.RDS"))) {
    # Get geo metadata
    geo.meta <- GEOquery::getGEO("GSE84465", destdir = "data")
    
    # Get the GSM and SRX ids from the GEO metadata
    id.key <- data.frame(gsm.ids = geo.meta[[1]]@phenoData@data$geo_accession,
    plate.ids = geo.meta[[1]]@phenoData@data$description.1,
    srx.ids = stringr::word(geo.meta[[1]]@phenoData@data$relation.1,
    start = 2, sep = "term="))
    # Import SRA file names
    sra.files <- read.csv(file.path("SRA.files.csv"))[, -1]
    
    # Merge both keys into a single id.key dataframe
    id.key <- merge(id.key, sra.files, by.x = "srx.ids", by.y = "experiment")
    
    # Save this dataframe for later
    saveRDS(id.key, file = file.path("sample_id_key.RDS"))
    
    # Save the rest of the metadata as a csv
    geo.meta <- data.frame(geo.meta[[1]]@phenoData@data)
    write.csv(geo.meta, file = file.path("meta_data.csv"))
    
    rm(geo.meta)
} else {
    id.key <- readRDS(file.path("sample_id_key.RDS"))
}

#----------------------Change sample column names to GSM-----------------------#
# Make our SRA conversion key
convert.key <- as.list(as.character(id.key$gsm.ids))
names(convert.key) <- as.character(id.key$run)

# Import datasets from their RDS files
salmon.data <- readRDS(file.path("salmon.data.RDS"))

# Obtain GSM ids using conversion key and make these the column names
colnames(salmon.data) <- dplyr::recode(colnames(salmon.data), !!!convert.key)

# Save objects with the GEO accession IDs
saveRDS(salmon.data, file.path("salmon.data.RDS"))
