#!/bin/bash
mkdir data/fastqc_trimmed

#--------------Paired trimming of the adapters from these data:----------------#
# Go to the raw_data directory 
cd data/raw_data

# Do adapter trimming for all sample file pairs
for f in `ls *_1.fastq.gz | sed 's/_1.fastq.gz//'`
do
/TrimGalore-0.4.5/trim_galore \
--nextera \
--paired ${f}_1.fastq.gz ${f}_2.fastq.gz \
-o ../fastqc_trimmed
done
