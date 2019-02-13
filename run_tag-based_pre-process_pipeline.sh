#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: running the pre-processing steps for tag-based single cell RNA-seq data.
# 
# Change your directory name and desired labels here.
dir=tabula_muris_data
label=tabula_muris

#----------------------------Make directories---------------------------------#
mkdir ${dir}
mkdir ref_files
mkdir results
mkdir ${dir}/alevin_output
mkdir ${dir}/normalized_${label}
mkdir ${dir}/raw_data

#---------------------------Download the dataset-------------------------------#
cd ${dir}

# If downloading from a url, use this: 
#curl -O http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar
#tar -xvf pbmc_10k_v3_fastqs.tar

# If downloading from NCBI's SRR run selector, use this: 
# Rscript scripts/pre-processing/0-get_sample_download_list.R \
#  -i ${SRP} \
#  -o results \
#  -d ${dir}/salmon_quants \
#  -q ref_files/SRAmetadb.sqlite

#------------------------------------------------------------------------------#
# Download each sample in the supplied list.
for line in `cat SRR_Acc_List_TabMur.txt`
  do
  fastq-dump --split-3 $line
  done
  
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
cd pbmc_10k_v3_fastqs
  
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
Rscript scripts/tag-based-pre-processing/1-make_gene_matrix-alevin.R \
  -d ${dir}/alevin_output \
  -o ${dir}/normalized_${label} \
  -q results \
  -l ${label}
