#!/bin/bash
# Use RSeQC to check mapping quality
mkdir results/map_qc

# Get genome bed file
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz'
gunzip hg38_RefSeq.bed.gz

# Get rRNA bedfile to use as quality control
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_rRNA.bed.gz'
gunzip hg38_rRNA.bed.gz

# Get stats on aligned reads
python /usr/local/bin/bam_stat.py -i ../aligned_reads/* -o ../../results/map_qc/

# Get stats on gene body coverage
python /usr/local/bin/geneBody_coverage.py \
-i aligned_reads \
-r hg38_RefSeq.bed \
-o ../results/map_qc

# For each file, get FPKM counts and split into rRNA and non rRNA
cd aligned_reads
for f in `ls -1 * | sed 's/.bam//' `
do
python /usr/local/bin/FPKM_count.py \
-i aligned_reads/${f}.bam \
-r ../hg38_RefSeq.bed \
-o ../../results/map_qc/${f}
python /usr/local/bin/split_bam.py \
-i ${f}.bam \
-r ../hg38_rRNA.bed \
-o ../../results/map_qc/${f}
done


