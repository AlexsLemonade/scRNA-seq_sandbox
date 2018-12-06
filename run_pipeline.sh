#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
mkdir data/aligned_reads
mkdir data/fastqc_reports
mkdir data/salmon_quants
mkdir results

#-------------------------- Index Genome for Salmon ---------------------------#
# Get the human transcriptome
curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-o data/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Index the human transcriptome
salmon --threads=16 --no-version-check index \
-t data/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i data/human_index \
-k 23

#--------------------------- Download fastq data-------------------------------#
# Note: because running all ~3800 samples takes quite a bit of time, use -n option 
# to set what number of randomly selected samples you would like  
Rscript scripts/0-get_sample_download_list.R \
-i SRP079058 \
-d data/salmon_quants \
-q data/SRAmetadb.sqlite \
-r

# Change directory to the data
cd data

# Download each sample and run Salmon on it. Then remove the sample to save room. 
for line in `cat ../files.to.download.txt`
do
  # Download forward and reverse fastq files
  Rscript ../scripts/1-download_sra.R -s $line -q SRAmetadb.sqlite -d raw_data
  # Run sequence quality control with FASTQC
  /FastQC/fastqc raw_data/* --outdir fastqc_reports
  # For each fastq file pair, do QC and then salmon
  for f in `ls raw_data/*_1.fastq.gz | sed 's/_1.fastq.gz//' `
  do
    "Processing sample ${f}"
    # Run Salmon
    salmon quant -i human_index -l A \
    -1 ${f}_1.fastq.gz \
    -2 ${f}_2.fastq.gz \
    -p 8 -o salmon_quants/$line \
    --gcBias --seqBias --biasSpeedSamp 5
  done
  rm raw_data/*
done

# Obtain summary report of fastqc:
Rscript scripts/2-get_fastqc_reports.R -d data/fastqc_reports -o results

# Make a gene matrix out of the Salmon quantification data
Rscript scripts/3-make_gene_matrix.R -d data/salmon_quants -o results
