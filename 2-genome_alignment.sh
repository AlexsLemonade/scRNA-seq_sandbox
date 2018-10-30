#!/bin/bash

# Get hisat2 for genome alignment
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip

# Get human genome. For other species, find the link at ftp://ftp.ncbi.nlm.nih.gov/genomes/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p12_top_level.gff3.gz


hisat2 -x ref_GRCh38.p12_top_level.gff3.gz -1 sample_1.fq.gz -2 sample_2.fq.gz
