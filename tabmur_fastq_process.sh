#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running pre-processing steps for tabula muris
# 
#----------------------------Make directories----------------------------------#
dir=tab_mur_data
label=tab_mur

mkdir tab_mur_data
mkdir tab_mur_data/tab_mur_bam
mkdir tab_mur_data/tab_mur_fastqs
mkdir tab_mur_data/alevin_output

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

#------------------------Download bam files from aws---------------------------#
cd tab_mur_data

# Make a list of the files we want
aws --no-sign-request s3 ls s3://czbiohub-tabula-muris/10x_bam_files/  file.list.txt

# Download each of them (NOTE: These are large files ~20GB each)
for line in $(cat file.list.txt | awk ‘{print $4}’)
  do aws --no-sign-request s3 cp s3://czbiohub-tabula-muris/10x_bam_files/$line tab_mur_bam
  done

#------------------------Convert bam files to fastq files----------------------#
cd tab_mur_bam

# Run bamtofastq on all the files in tab_mur_bam folder
for file in `ls *bam`
  do
  cellranger bamtofastq $file ../tab_mur_fastqs
  done

#-------------------------Get rid of extra folders-----------------------------#
cd ../tab_mur_fastqs

# Get rid of those extra folders that bamtofastq makes  
for folder in `ls`
  do
  mv `ls $folder`/* .
  rmdir `ls $folder`
  done
  
#-------------------------------Quantify samples-------------------------------#
# For each fastq file pair run salmon/alevin for quantfication
# Change to the individual folder and then run Alevin for all the files there
for folder in `ls`
  do
  cd ${folder}
  for f in `ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' `
    do
    echo "Processing sample ${f}"
    # Run Salmon with Alevin
    salmon alevin -l ISR  \
      --no-version-check \
      -i ../../../ref_files/human_index \
      -1 ${f}_R1_00*.fastq.gz \
      -2 ${f}_R2_00*.fastq.gz \
      --chromium  \
      -p 10 \
      -o ../../alevin_output/${f} \
      --tgMap ../../../ref_files/genes_2_tx.tsv \
      --dumpCsvCounts \
      --dumpFeatures
    done
  cd ..
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
