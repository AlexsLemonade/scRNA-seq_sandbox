# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: After salmon has been run successfully on your samples. Assemble the
# Individual samples quantification data into a matrix. Also use the ensembl IDs
# to obtain gene symbols. Then save this matrix to an tsv file.

# Options:
# "-d" - Directory of where individual samples' salmon folders are located.(Optional)
# "-o" - Directory of where the output gene matrix RDS file should go.(Optional)
# "-g" - GEO ID for the dataset so the metadata can be downloaded.
# "-m" - Percent mapped reads (reported as a decimal) cutoff for filtering 
#        samples. Default is 0.5.
# "-l" - Optional label to add to output files. Generally necessary if processing
#        multiple datasets in the same pipeline.
# 
# Command line example:
#
# Rscript scripts/3-make_gene_matrix.R \
# -d data/salmon_quants \
# -o data \
# -m 0.5 \
# -l "patel"

#-------------------------- Get necessary packages-----------------------------#
if (!("org.Hs.eg.db" %in% installed.packages())) { 
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db", suppressUpdates = TRUE)
}
if (!("tximport" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("tximport", suppressUpdates = FALSE)
}
if (!("rjson" %in% installed.packages())) {
  install.packages("rjson", suppressUpdates = FALSE)
}

# Attach needed libraries
library(org.Hs.eg.db)
library(optparse)
library(tximport)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

#--------------------------------Set up options--------------------------------#
option_list <- list(
    make_option(opt_str = c("-d", "--dir"), type = "character", default = getwd(),
                help = "Directory where salmon quantification folders are located",
                metavar = "character"),
    make_option(opt_str = c("-o", "--output"), type = "character",
                default = getwd(), help = "Directory where you would like the
                output to go", metavar = "character"),
    make_option(opt_str = c("-m", "--mapped"), type = "numeric", 
                default = "0.5", help = "Cutoff for what percent mapped_reads samples
                should have. Any samples with less than the cutoff will be removed.",
                metavar = "numeric"),
    make_option(opt_str = c("-l", "--label"), type = "character",
                default = "", help = "Optional label for output files",
                metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

opt$data <- "patel_data/salmon_quants"
opt$output <- "results"
opt$label <- "patel_data"
#------------------------------Import Salmon reads-----------------------------#
# Get quant files
quant.files <- list.files(opt$dir, recursive = TRUE, full.names = TRUE,
                          pattern = "quant.sf")

# Carry sample names over
sample.names <- stringr::word(quant.files, 3, sep = "/")

# Get transcript IDs
transcripts <- read.table(quant.files[[1]], header = TRUE, 
                          colClasses = c("character", rep("NULL", 4)))

# Get rid of transcript version numbers because mapIDs will not recognize them
transcripts.wout.ver <- gsub("\\.[0-9]*$", "", transcripts$Name)

# Turn transcript ids to gene ids
tx.gene.key <- data.frame('transcript' = transcripts$Name,
                          'gene' = mapIds(org.Hs.eg.db, keys = transcripts.wout.ver,
                                         column = "ENSEMBL", keytype = "ENSEMBLTRANS"),
                          stringsAsFactors = FALSE)

# Remove transcripts without genes
tx.gene.key <- tx.gene.key[!is.na(tx.gene.key$gene),]

# Do the thing
tx.import <- tximport::tximport(quant.files, type = "salmon", tx2gene = tx.gene.key,
                                countsFromAbundance = "no")

# Make a counts dataframe
tx.counts <- data.frame('gene' = rownames(tx.import),
                        tx.counts$counts, stringsAsFactors = FALSE)

# Make a tpms dataframe
tx.tpm <- data.frame('gene' = rownames(tx.import),
                     tx.import, stringsAsFactors = FALSE)

#---------------------Salmon proportion of mapped reads------------------------#
# Get the proportion of mapped reads by reading the meta files
salmon.prop.assigned <- vapply(sample.names, function(x) {
  rjson::fromJSON(file = file.path(opt$dir, x, "aux_info",
                                   "meta_info.json"))$percent_mapped/100
  }, FUN.VALUE = 1)

# Make a histogram of this information
png(file.path(opt$output, paste0(opt$label, "_prop_reads_mapped_hist.png")))
hist(salmon.prop.assigned, xlab = "", main = "Proportion of Mapped Reads",
     breaks = 20)
abline(v = opt$mapped, col = "red")
dev.off()

# Filter out samples with too low of mapped reads
tx.counts <- tx.counts[, which(salmon.prop.assigned > opt$mapped)]

# Print report of how many samples are left
cat("Number of samples that have percent mapped reads greater than cutoff:",
    length(which(salmon.prop.assigned > opt$mapped)))

# Save both tpms and counts to tsv files
readr::write_tsv(tx.counts, file.path(opt$output,
                                      paste0(opt$label, "_counts.tsv")))
readr::write_tsv(tx.tpm, file.path(opt$output,
                                      paste0(opt$label, "_tpms.tsv")))

