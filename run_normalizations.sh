#!/bin/bash
# C. Savonen
# CCDL for ALSF 2018

# Purpose: running the normalization on gene matrices previously made
# by the run_pipeline.sh scripts

# The three normalization methods: seurat, scnorm, and RUVnormalize

# Process the patel data in the three ways
Rscript scripts/4a-seurat_normalize.R \
-d patel_data/patel.counts.tsv \
-o normalized_data \
-l "patel"

Rscript scripts/4b-scnorm_normalize.R \
-d patel_data/patel.counts.tsv \
-o normalized_data \
-l "patel"

# Process the darmanis data in the three ways
Rscript scripts/4a-seurat_normalize.R \
-d darmanis_data/darmanis.counts.tsv \
-o normalized_data \
-l "darmanis"

Rscript scripts/4b-scnorm_normalize.R \
-d darmanis_data/darmanis.counts.tsv \
-o normalized_data \
-l "darmanis"
