#!/bin/bash
#--------------- Align trimmed reads to genome using HISAT2--------------------#
# Get premade human genome index:
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz -P ../
tar data/grch38.tar.gz

# Go to trimmed reads directory
cd ../fastqc_trimmed

# Do genome alignment for all sample pairs
for f in `ls -1 *_1_val_1.fq.gz | sed 's/_1_val_1.fq.gz//' `
do
/hisat2-2.1.0/hisat2 -x ../grch38/genome \
-1 ${f}_1_val_1.fq.gz \
-2 ${f}_2_val_2.fq.gz \
-S ../aligned_reads/${f}.bam
done
