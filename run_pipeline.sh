#!/bin/bash

# Cannot run this script before 0-set_up.sh is ran, but should only need to run
# that script once

# Set directory
cd home/rstudio/kitematic/scRNA-seq_workflow

# Make directories
mkdir fastqc_data
mkdir fastqc_trimmed_results

# Download fastq data
Rscript 1-download_fastq_data.R

# Run fastqc and get reports
FASTQC/fastqc raw_data/* --outdir fastqc_data
Rscript 2-get_fastqc_reports.R

# Trim the adapters from these data:
TrimGalore-0.4.5/trim_galore --nextera -o fastqc_trimmed_results fastqc_data/*.zip
