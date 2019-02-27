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
# Make a list of the files we want
aws --no-sign-request s3 ls s3://czbiohub-tabula-muris/10x_bam_files/10X_P4_ > file.list.txt
# grep -w "10X_P4_[5-7]_possorted_genome.bam$" file.list.txt

# Download each of them (NOTE: These are large files ~20GB each)
for line in $(grep -w "10X_P4_[6-7]_possorted_genome.bam$" file.list.txt | awk '{print $4}')
  do 
  echo "Downloading bam file: $line"
  aws --no-sign-request s3 cp s3://czbiohub-tabula-muris/10x_bam_files/$line tab_mur_bam
  echo "Converting $line to fastq file"
  cellranger bamtofastq tab_mur_bam/$line tab_mur_fastqs/$line
  echo "Removing $line"
  rm tab_mur_bam/$line
  done
  
#-------------------------------Quantify samples-------------------------------#
# Change to where the data are:
cd tab_mur_fastqs

# Set up a file to hold the times in seconds for file processing: 
echo "File time(secs)" > ../file.run.duration.txt

# For each fastq file sets run Salmon/Alevin for quantfication
for folder in `ls`
  do
  cd $folder/`ls $folder`/
  for f in `ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' `
    do
    echo "Processing sample ${folder}/${f}"
    start=`date +%s`
    echo "Start time : $start"
    salmon alevin -l ISR  \
      --no-version-check \
      -i ../../../human_index \
      -1 ${f}_R1_00*.fastq.gz \
      -2 ${f}_R2_00*.fastq.gz \
      --chromium  \
      -p 10 \
      -o ../../alevin_output/${folder}/${f} \
      --tgMap ../../../genes_2_tx.tsv \
      --dumpCsvCounts \
      --dumpFeatures
    done
    end=`date +%s`
  cd ../..
  echo "${f} $((end-start)) seconds" >> ../file.run.duration.txt
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
