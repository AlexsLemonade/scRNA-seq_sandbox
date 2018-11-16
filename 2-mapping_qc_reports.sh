#!/bin/bash
# Purpose: Use RSeQC to check mapping quality

# Get genome bed file
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz'
gunzip hg38_RefSeq.bed.gz

# Get rRNA bedfile to use as quality control
wget 'https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_rRNA.bed.gz'
gunzip hg38_rRNA.bed.gz

mv data/aligned_reads/*.sorted.bam data/sorted_bams/*.sorted.bam
#ls aligned_reads/*.sorted.bam.bai | rename  's/.sorted.bam.bai/.bam.bai/'

cd aligned_reads
# For each file, get FPKM counts
for f in `ls *.bam | sed 's/.bam//' `
do
python /usr/local/bin/FPKM_count.py \
-i ${f}.bam \
-r ../hg38_RefSeq.bed \
-o ../fpkms/${f}c
done

# Get stats on aligned reads
python /usr/local/bin/geneBody_coverage.py \
-r hg38_RefSeq.bed \
-i aligned_reads/*.bam \
-o ../results/map_qc/


for f in `ls -1 * | sed 's/.bam//' `
do
python /usr/local/bin/split_bam.py \
-i ${f}.bam \
-r ../hg38_rRNA.bed \
-o ../../results/map_qc/${f}
done


# Install cmake so we can build Salmon tools later
RUN wget http://www.cmake.org/files/v2.8/cmake-2.8.3.tar.gz && \
tar xzf cmake-2.8.3.tar.gz && \
cd cmake-2.8.3 && \
./configure --prefix=/opt/cmake && \
make && \
make install

# Need cmake for Salmon to install
RUN cd / && \
wget https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.tar.gz && \
tar -zxvf cmake-3.11.4-Linux-x86_64.tar.gz && \
cd cmake-3.11.4-Linux-x86_64 && \
mkdir build
cd build
~/cmake-3.11.4-Linux-x86_64/bin/cmake ..

RUN cd / && \
wget https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.tar.gz && \
tar xzf cmake-3.11.4-Linux-x86_64.tar.gz && \
cd cmake-3.11.4-Linux-x86_64 && \
./configure --prefix=/opt/cmake && \
make && \
make install
