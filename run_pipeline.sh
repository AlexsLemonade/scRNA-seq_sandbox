#!/bin/bash
# Note, change the directory to where you want these things to appear before running this code.

# Make directories
mkdir STAR
mkdir fastqc_reports
mkdir fastqc_trimmed

# Note: because running all ~3800 samples takes quite a bit of time, this 
# example run is set to only run 10 samples. When you want to run the full set, 
# delete the -n argument.  

# Download fastq data
Rscript 0-download_fastq_data.R -d SRP079058 -n 10

# Run fastqc and get reports
/bin/FastQC/fastqc raw_data/* --outdir fastqc_reports
Rscript 1-get_fastqc_reports.R -d ./fastqc_reports

# Trim the adapters from these data:
/bin/TrimGalore-0.4.5/trim_galore --nextera -o ./fastqc_trimmed ./raw_data/*

# Run genome alignmnent with STAR - This is step is not set up yet
# bash 2-genome_alignment.sh
