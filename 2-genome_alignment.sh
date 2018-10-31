#!/bin/bash

# Get hisat2 for genome alignment
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip

# Get human genome. For other species, find the link at ftp://ftp.ncbi.nlm.nih.gov/genomes/
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p12_top_level.gff3.gz

# Get premade human genome index:
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz

hisat2-2.1.0/hisat2-build ref_GRCh38.p12_top_level.gff3.gz

hisat2 -x data/grch38.tar.gz -1 data/raw_data/SRR3934448_1.fastq.gz -2 data/raw_data?SRR3934448_1.fastq.gz
