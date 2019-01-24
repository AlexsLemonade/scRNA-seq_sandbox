# C. Savonen 
# CCDL for ALSF 
# 2019

# Purpose: Make a human transcript-gene key list 

# Two ways to do it: 
##########################################################################
# Method 1) Get the genes and transcripts from AnnotationDbi R packages

# Get ensembl genes
genes <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENSEMBL")

# Retrieve transcripts for those
gene_2_tx <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                  keys = genes, 
                                  column = "ENSEMBLTRANS", 
                                  keytype = "ENSEMBL", 
                                  multiVals="CharacterList")

# Make into dataframe:
gene_2_tx <- reshape2::melt(gene_2_tx@listData)

# Get rid of transcripts without genes
gene_2_tx <- gene_2_tx[!is.na(gene_2_tx$value), ]

# Write to tsv file
readr::write_tsv(gene_2_tx, "genes_2_tx.tsv", col_names = FALSE)

##########################################################################
# Method 2) Get the genes and transcripts directly from the transcriptome file

# Get the ensembl transcriptome file
download.file(
"ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz")

# Unzip the file
R.utils::gunzip("ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz")

# Read in ensembl transcriptome
genome <- readLines("ref_files/Homo_sapiens.GRCh38.cdna.all.fa")
genome <- genome[grep(">", genome)]

# Extract transcripts
transcripts <- stringr::word(genome, sep = " ", 1)
transcripts <- gsub(">", "", transcripts)

# Extract genes
genes <- stringr::word(genome, sep = "gene:", 2)
genes <- stringr::word(genes, sep = " ", 1)
# genes <- gsub("\\.[0-9]*$", "", genes)

# Write this to a tsv file
readr::write_tsv(data.frame(transcripts, genes), "genes_2_tx.tsv",
                 col_names = FALSE)
