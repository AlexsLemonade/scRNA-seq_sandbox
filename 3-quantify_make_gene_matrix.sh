#!/bin/bash

# Use aligned reads to get quantification with Salmon
./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant

# This *should* Run all the samples for us
for f in `ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//' `
do
echo ./bin/salmon quant -t transcripts.fa -l -1 ${f}_1.fastq.gz -2 ${f}_2.fastq.gz -S ${f}.bam
done
