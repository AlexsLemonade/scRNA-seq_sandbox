#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir data
mkdir data/raw_data
mkdir data/aligned_reads
mkdir data/fastqc_reports
mkdir data/fastqc_trimmed

# Note: because running all ~3800 samples takes quite a bit of time, this 
# example run is set to only run 10 samples. When you want to run the full set, 
# delete the -n argument.  

# Download fastq data
Rscript 0-download_fastq_data.R -i SRP079058 -n 10 -d data/raw_data

# Run fastqc and get reports
/bin/FastQC/fastqc raw_data/* --outdir data/fastqc_reports
Rscript 1-get_fastqc_reports.R -d data/fastqc_reports

# Trim the adapters from these data:
/bin/TrimGalore-0.4.5/trim_galore --nextera -o data/fastqc_trimmed data/raw_data/*

# Run genome alignmnent with STAR - This is step is not set up yet
bash 2-genome_alignment.sh
