#!/bin/bash
# Use RSeQC to check mapping quality
mkdir results/map_qc

# Get genome bed file and rRNA bedfile
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz'
gunzip hg38_RefSeq.bed.gz
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_rRNA.bed.gz'
gunzip hg38_rRNA.bed.gz

cd aligned_reads
for f in `ls -1 * | sed 's/.bam//' `
do
python /usr/local/bin/geneBody_coverage.py \
-i aligned_reads/${f}.bam \
-r ../hg38_RefSeq.bed \
-o ../../results/map_qc/${f}_geneBody.txt
done
python /usr/local/bin/bam_stat.py \
-i ${f}.bam -r ../hg38_RefSeq.bed \
-o ../../results/map_qc/${f}bam_stat.txt
python /usr/local/bin/split_bam.py \
-i ${f}.bam \
-r ../hg38_rRNA.bed \
-o ../../results/map_qc/${f}split_bam.txt
done


