#!/bin/bash

# Get hisat2 for genome alignment
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip

# Get premade human genome index:
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz

# Build the genome
/hisat2-2.1.0/hisat2-build data/grch38.tar.gz ht25

# Using the genome index file for human, align a single sample to the genome
/hisat2-2.1.0/hisat2 -x ht25 \
-1 data/fastqc_trimmed/SRR3934448_1.fastq.gz \
-2 data/fastqc_trimmed/SRR3934448_2.fastq.gz

# This *should* Run all the samples for us
for f in `ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//' `
do
echo /hisat2-2.1.0/hisat2 -x ht25 -1 ${f}_1.fastq.gz -2 ${f}_2.fastq.gz -S ${f}.bam
done



