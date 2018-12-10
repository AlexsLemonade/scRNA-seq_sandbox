#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018
# Purpose: running the pre-processing steps for single cell data.

#-------------------------------Make directories-------------------------------#
mkdir anoop_data
mkdir anoop_data/raw_data
mkdir anoop_data/aligned_reads
mkdir anoop_data/fastqc_reports
mkdir anoop_data/salmon_quants
mkdir results
mkdir ref_files

#-------------------------- Index Genome for Salmon ---------------------------#
# Will check for genome index first before running
if [ ! -e ref_files/human_index ]; then
# Get the human transcriptome
curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-o ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Index the human transcriptome
salmon --threads=16 --no-version-check index \
-t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i ref_files/human_index \
-k 23
fi

#--------------------------- Download fastq anoop_data-------------------------------#
# Note: because running all ~3800 samples takes quite a bit of time, use -n option 
# to set what number of randomly selected samples you would like  
Rscript scripts/0-get_sample_download_list.R \
-i SRP042161 \
-d anoop_data/salmon_quants \
-q ref_files/SRAmetadb.sqlite \
-r

# Change directory to the data
cd anoop_data

# Download each sample and run Salmon on it. Then remove the sample to save room. 
for line in `cat ../files.to.download.txt`
do
  # Download forward and reverse fastq files
  Rscript ../scripts/1-download_sra.R \
  -s $line \
  -q ../ref_files/SRAmetadb.sqlite \
  -d raw_data

  # Run sequence quality control with FASTQC
  /FastQC/fastqc raw_data/* --outdir fastqc_reports

  # For each fastq file pair, do QC and then salmon
  for f in `ls raw_data/*_1.fastq.gz | sed 's/_1.fastq.gz//' `
  do
    echo "Processing sample ${f}"

    # Run Salmon
    salmon quant -i ../ref_files/human_index -l A \
    -1 ${f}_1.fastq.gz \
    -2 ${f}_2.fastq.gz \
    -p 8 -o salmon_quants/$line \
    --gcBias --seqBias --biasSpeedSamp 5
  done
  rm raw_data/*
done

# Obtain summary report of fastqc:
Rscript scripts/2-get_fastqc_reports.R \
-d anoop_data/fastqc_reports \
-o results

# Make a gene matrix out of the Salmon quantification data
Rscript scripts/3-make_gene_matrix.R \
-d anoop_data/salmon_quants \
-g GSE57872
-o anoop_data
