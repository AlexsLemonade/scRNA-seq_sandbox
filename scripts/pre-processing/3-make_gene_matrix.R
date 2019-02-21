# CCDL for ALSF 2018
# C. Savonen
#
# Purpose: After salmon has been run successfully on your samples, this script
# will use tximport to quantify by transcripts and put all your samples together
# in one gene matrix tsv file.

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
# -g GSE86445 \
# -m 0.5 \
# -l "patel"

#-------------------------- Get necessary packages-----------------------------#
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
    make_option(opt_str = c("-g", "--geo"), type = "character", default = getwd(),
                help = "GEO ID of dataset for downloading the metadata",
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

# Add an underscore if label is specified
if (!is.null(opt$label)){
  opt$label <- paste0(opt$label, "_")
}

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
                          'gene' = mapIds(org.Hs.eg.db, 
                                          keys = transcripts.wout.ver,
                                          column = "ENSEMBL", 
                                          keytype = "ENSEMBLTRANS"),
                          stringsAsFactors = FALSE)

# Remove transcripts without genes
tx.gene.key <- tx.gene.key[!is.na(tx.gene.key$gene),]

# Do the thing
tx.counts <- tximport::tximport(quant.files, type = "salmon", 
                                tx2gene = tx.gene.key,
                                countsFromAbundance = "no")

# Save to RDS file temporarily
saveRDS(tx.counts, "tximport_obj.RDS")
# tx.counts <- readRDS("tximport_obj.RDS")

# Make as a dataframe
tx.counts <- data.frame(tx.counts$counts, stringsAsFactors = FALSE)

# Bring the sample names with
colnames(tx.counts) <- sample.names

#---------------------Salmon proportion of mapped reads------------------------#
# Get the proportion of mapped reads by reading the meta files
salmon.prop.assigned <- vapply(sample.names, function(x) {
  rjson::fromJSON(file = file.path(opt$dir, x, "aux_info",
                                   "meta_info.json"))$percent_mapped/100
  }, FUN.VALUE = 1)

# Make a histogram of this information
png(file.path(opt$output, paste0(opt$label, "prop_reads_mapped_hist.png")))
hist(salmon.prop.assigned, xlab = "", main = "Proportion of Mapped Reads",
     breaks = 20)
dev.off()

# Filter out samples with too low of mapped reads
tx.counts <- tx.counts[, which(salmon.prop.assigned > opt$mapped)]

# Make a gene column so read_tsv will have the info
tx.counts <- tx.counts %>% tibble::rownames_to_column("gene")

# Save to tsv file
readr::write_tsv(tx.counts, file.path(opt$output,
                                      paste0(opt$label, "counts.tsv")))

