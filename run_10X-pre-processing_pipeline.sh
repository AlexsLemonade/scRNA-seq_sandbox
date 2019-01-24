#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: running the pre-processing steps for single cell RNA-seq data.
# 
# Change your directory name and desired labels here.
dir=pbmc_data
label=pbmc_1k_v2

#----------------------------Make directories---------------------------------#
mkdir ${dir}
mkdir ref_files
mkdir results
mkdir ${dir}/alevin_output
mkdir ${dir}/normalized_${label}

#---------------------------Download the dataset-------------------------------#
cd ${dir}
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_fastqs.tar
tar -xvf ${label}.tar

#--------------Will run setup only if it hasn't been ran before----------------#
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

#--------------------------- Quantify samples-------------------------------#
# For each fastq file pair run salmon/alevin for quantfication
cd ${label}
  
for f in `ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' `
  do
  echo "Processing sample ${f}"

  # Run Salmon with Alevin
  salmon alevin -l ISR  \
    -i ../../ref_files/human_index \
    -1 ${f}_R1_001.fastq.gz \
    -2 ${f}_R2_001.fastq.gz\
    --chromium  \
    -p 10 \
    -o ../alevin_output/${f} \
    --tgMap ../../ref_files/genes_2_tx.tsv \
    --dumpCsvCounts \
    --dumpFeatures
  done

#-------------------Make gene matrix from alevin output------------------------#
# Back out of the current directory
cd ../..

# Run this script, which will create alevinQC output for you
Rscript scripts/10x-pre-processing/1-make_gene_matrix-alevin.R \
  -d ${dir}/alevin_output \
  -o ${dir}/normalized_${label}/counts_${label}.tsv \
  -q ../results \
  -l ${label}
