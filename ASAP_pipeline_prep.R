# C.Savonen
# CCDL for ALSF 2018
# 
# Prepping your data for use in ASAP online pipeline

# Import gene expression matrix data
salmon.data <- readRDS(file.path("data", "salmon.data.RDS"))

# Batch effect testing
geo.meta <- read.csv(file.path("data", "meta_data.csv"))
geo.meta <- geo.meta[match(colnames(salmon.data),geo.meta$geo_accession),]

#------------------------------Prepping data for ASAP--------------------------#
# Format the data so ASAP will like it. 
salmon.data <- apply(salmon.data[,-c(1:2)], 2, round)

# Write these data to tab delimted text files for use in ASAP
write.table(salmon.data, file.path("data", "salmon.data.txt"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write batch files
geo.meta <- GEOquery::getGEO("GSE84465", destdir = "data")
geo.meta <- data.frame(geo.meta[[1]]@phenoData@data)
geo.meta <- geo.meta[match(colnames(salmon.data)[-c(1:2)], geo.meta$geo_accession), ]

# Write batch info file for ASAP use
write.table(cbind(geo.meta$geo_accession, geo.meta$plate.id.ch1),
            file.path("data", "batch.info.txt"), row.names = FALSE, 
            col.names = FALSE, sep = "\t", quote = FALSE)
