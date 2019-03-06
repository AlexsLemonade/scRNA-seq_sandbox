#!/bin/bash
# C. Savonen
# CCDL for ALSF 2019

# Purpose: running pre-processing steps for Tabula Muris data
#
#----------------------------Make directories----------------------------------#
mkdir tab_mur_data
mkdir tab_mur_data/tab_mur_bam
mkdir tab_mur_data/tab_mur_fastqs
mkdir tab_mur_data/alevin_output

# Used the same transcriptome index as has been used before.
# Made using these steps:
# curl ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
#   -o ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz
#
# salmon --threads=16 --no-version-check index \
#   -t ref_files/Homo_sapiens.GRCh38.cdna.all.fa.gz \
#   -i ref_files/human_index \
#   -k 23

#------------------------Download bam files from aws---------------------------#
cd tab_mur_data

# Make a list of the files we want
aws --no-sign-request s3 ls s3://czbiohub-tabula-muris/10x_bam_files/10X_P4_ > file.list.txt

# Select only the bam files, and not the bam.bai
grep -w "10X_P4_[0-7]_possorted_genome.bam$" file.list.txt > file.list.txt

# Download each of them (NOTE: These are large files ~20GB each)
for file in $(file.list.txt | awk '{print $4}')
  do
  echo "Downloading bam file: $file"
  aws --no-sign-request s3 cp s3://czbiohub-tabula-muris/10x_bam_files/$file tab_mur_bam
  echo "Converting $file to fastq file"
  cellranger bamtofastq tab_mur_bam/$file tab_mur_fastqs/$file
  echo "Removing $file"
  rm tab_mur_bam/$file
  done

#-------------------------------Quantify samples-------------------------------#
# Change to where the data are:
cd tab_mur_fastqs

# Set up a file to hold the times in seconds for file processing:
echo "File time (secs)" > ../file.run.duration.txt

# For each fastq file sets run Salmon/Alevin for quantfication
for folder in `ls`
  do
  cd $folder/`ls $folder`/
  for file in `ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' `
    do
    echo "Processing sample ${folder}/${file}"
    start=`date +%s`
    echo "Start time : $start"
    salmon alevin -l ISR  \
      --no-version-check \
      -i ../../../human_index \
      -1 ${file}_R1_00*.fastq.gz \
      -2 ${file}_R2_00*.fastq.gz \
      --chromium  \
      -p 10 \
      -o ../../alevin_output/${folder}/${file} \
      --tgMap ../../../genes_2_tx.tsv \
      --dumpCsvCounts \
      --dumpFeatures
    done
    end=`date +%s`
  cd ../..
  echo "${file} $((end-start)) seconds" >> ../file.run.duration.txt
  done

#-------------------Make gene matrix from alevin output------------------------#
# Back out of the current directory
cd ../..

# Run this script, which will create alevinQC output for you (if you designate
# an output folder for it with the -q option )
Rscript scripts/tag-based-pre-processing/1-make_gene_matrix-alevin.R \
  -d tab_mur_data/alevin_output \
  -o tab_mur_data/normalized_tab_mur \
  -l tab_mur \
  -r \
  -f
