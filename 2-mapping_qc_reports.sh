#!/bin/bash
# Use RSeQC to check mapping quality

# Get human gene model
wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/
wget http://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/GRCh38_rRNA.bed.gz/download

bam_stat.py  -i aligned_reads/*
