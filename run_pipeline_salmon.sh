#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
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
Rscript 0-download_fastq_data.R \
-i SRP079058 \
-d data/raw_data 

cd data

for line in `cat ../files.2.download.txt`
do
  # Download forward and reverse fastq files
  Rscript ../download_sra.R -s $line -q SRAmetadb.sqlite -d raw_data
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
Rscript 1-get_fastqc_reports.R -d data/fastqc_reports -o results
