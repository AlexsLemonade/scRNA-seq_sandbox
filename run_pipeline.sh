#!/bin/bash

# Cannot run this script before 0-set_up.sh is ran, but should only need to run
# that scriptonce

# Set directory
cd home/rstudio/kitematic/scRNA-seq_workflow

# Make directories
mkdir STAR
mkdir fastqc_reports
mkdir fastqc_trimmed

# Download fastq data
Rscript 1-download_fastq_data.R -d SRP079058 -n 100

# Run fastqc and get reports
FASTQC/fastqc raw_data/* --outdir fastqc_reports
Rscript 2-get_fastqc_reports.R

# Trim the adapters from these data:
/bin/TrimGalore-0.4.5/trim_galore --nextera -o fastqc_trimmed raw_data/*

# Run genome alignmnent with STAR
bash 3-genome_alignment.sh
