#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running the post-processing steps for single cell RNA-seq data.

# Change your directory name, GEO ID, and SRP here. Then run the script.
dir=darmanis_data
GSE=GSE84465
SRP=SRP079058
label=darmanis

#-------------------------------Run normalization------------------------------#
Rscript scripts/5-run_normalization.R \
  -d ${dir}/normalized_${label}/counts_${label}.tsv \
  -a all \
  -o ${dir}/normalized_${label} \
  -l ${label}

#------------------------------Dimension reduction-----------------------------#
Rscript scripts/6-dim_reduction_analysis.R \
  -d ${dir}/normalized_${label} \
  -m ${dir}/metadata.tsv \
  -r pca \
  -l ${label} \
  -o results/pca_${label} 
  
#------------------------------Clustering analysis-----------------------------#
Rscript scripts/7-cluster_analysis.R \
  -d results/pca_${label} 
  -m ${dir}/metadata.tsv \
  -l ${label} \
  -o results/pca_${label} 