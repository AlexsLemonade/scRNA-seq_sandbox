
# Import libraries
import os as os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rcParams

# Declare directory path
input_dir = 'Users/candacesavonen/Desktop/GitRepos/scRNA-seq_workflow/tab_mur_data'

# Declare input file paths
data_file = os.path.join(input_dir, 'filtered_counts_tab_mur.tsv')
metadata_file = os.path.join(input_dir, 'filtered_metadata_tab_mur.tsv')

# Read in the gene expression data
sc_data = sc.read(data_file, delimiter='\t')

# Read in the metadata file
metadata = pd.read_csv(metadata_file, delimiter="\t")

# Set the tissue information as the Anndata observation catregory
sc_data.obs['tissue'] = metadata['tissue']  

# Apply pre-processing
sc.pp.recipe_zheng17(sc_data)

# Compute PCA
sc.tl.pca(sc_data, svd_solver='arpack')

# Find nearest neighbors
sc.pp.neighbors(sc_data, n_neighbors=4, n_pcs=20)

# Draw the graph
sc.tl.draw_graph(sc_data)


sc_data.write(filename)
sc_data.write_csvs(filename)

