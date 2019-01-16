#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running the post-processing steps for single cell RNA-seq data.

# Change your directory name, and desired label here. Then run the script.
dir=darmanis_data
label=darmanis

#---------------------------Prep the data with this Rmd------------------------#
# Need to set up this Rmd for each dataset
Rscript -e "rmarkdown::render('darmanis_data_prep.Rmd')"

#-------------------------------Run normalization------------------------------#
Rscript scripts/post-processing/5-run_normalization.R \
  -d ${dir}/normalized_${label}/counts_${label}.tsv \
  -a all \
  -o ${dir}/normalized_${label} \
  -l ${label}

#------------------------------Dimension reduction-----------------------------#
Rscript scripts/post-processing/6-dim_reduction_analysis.R \
  -d ${dir}/normalized_${label} \
  -m ${dir}/metadata.tsv \
  -r pca \
  -l ${label} \
  -o results/pca_${label} 
  
#------------------------------Clustering analysis-----------------------------#
Rscript scripts/post-processing/7-cluster_analysis.R \
  -d results/pca_${label} \
  -m ${dir}/metadata.tsv \
  -l ${label} \
  -o results/pca_results 