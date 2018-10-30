#!/bin/bash

# Set directory
cd home/rstudio/kitematic/scRNA-seq_workflow

# Run STAR to align to genome 
STAR --runThreadN 4 \
--genomeDir indices/STAR \
--readFilesIn data/fastqc_trimmed/* \
--outFileNamePrefix data/STAR/ \

# Prep the human genome reference sequences
mkdir -p data/GRCh38/sequence
cd data/GRCh38/sequence/
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r77.all.fa
cd ../../../

# Prep human genome annotations
mkdir -p data/GRCh38/annotation
cd data/GRCh38/annotation/
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip Homo_sapiens.GRCh38.77.gtf.gz
cd ../../../
mkdir -p data/GRCh38/star_indices_overhang100

bin/star --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir data/GRCh38/star_indices_overhang100/ \
--readFilesIn fastqc_trimmed/* \
--outFileNamePrefix STAR/ \
--sjdbGTFfile data/GRCh38/annotation/Homo_sapiens.GRCh38.77.gtf \
--sjdbOverhang 100
