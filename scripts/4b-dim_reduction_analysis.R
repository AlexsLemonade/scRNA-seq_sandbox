# C.Savonen, CCDL for ALSF
# 2019
# 
# Purpose: Using Seurat for post-processing of scRNA-seq gene matrix .tsv 
#          In gene x sample format. 
#
# Options:
# "-d" - Directory of where dataset(s) to analyze exist
# "-m" - file of metadata information to analyze. Each column in this metadata 
#        will be analyzed. Input should be a path to a tsv file 
# 
# Command line example:
#
# Rscript scripts/4a-seurat_normalize.R \
# -d data/salmon_quants \
# -o data \
# -l "patel"
# 
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Import functions for this analysis
source(file.path("scripts", "util", "clustering_statistics_functions.R"))

# ggplot2 library
library(ggplot2)
library(optparse)

# These are the necessary packages
packages <- c("Rtsne", "NMI", "caret")

# Check if these packages are installed and install them if they aren't
lapply(packages, function(package) {
  if (!(package %in% installed.packages())) {
    install.packages(package)
  }
})
#--------------------------------Set up options--------------------------------#
option_list <- list(
  make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
              help = "Path to gene matrix in tab delimited format, gene x sample
              with gene info as the first column",
              metavar = "character"),
  make_option(opt_str = c("-m", "--metadata"), type = "character",
              default = "none", help = "Path to metadata file",
              metavar = "character"),
  make_option(opt_str = c("-r", "--reduce"), type = "character",
              default = "none", help = "Dimension reduction technique to use.
              Options are: 'pca', 'tsne', or 'umap' If none is given, full 
              datasets will be used.", metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character",
              default = getwd(), help = "Directory where you would like the
              output to go", metavar = "character"),
  make_option(opt_str = c("-l", "--label"), type = "character",
              default = "", help = "Optional label for output files",
              metavar = "character")
  )

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

#--------------------------------Set up options--------------------------------#

opt$metadata <- file.path("darmanis_data", "meta_data.csv")
opt$data <- file.path("darmanis_data", "normalized_darmanis")
opt$label <- "darmanis"
opt$output <- "results"
opt$reduce <- "tsne"

# Check that the dimension reduction option given is supported
if (!(opt$reduce %in% c("tnse", "pca", "umap", "none"))){
  stop("That is not a dimension reduction technique supported by this script. 
       Check for typos. Acceptable options:'tsne', 'pca', or 'umap'")
}

# Create the output for results folder if it does not exist
if (!dir.exists(opt$output)) {
  message(paste("Output folder:", opt$output, "does not exist, creating one"))
  dir.create(opt$output)
}

#---------------------------------Read in data---------------------------------#
# Get the file names of all the normalized files
dataset.files <- dir(opt$data, full.names = TRUE)

# Read in each of the normalization files
datasets <- lapply(dataset.files, readr::read_tsv)

# Keep all the gene lists
genes <- lapply(datasets, function(x) x[, 1])

# Get names of datasets
dataset.names <- dir(opt$data)

# Remove file extension name and make it the name in the list
names(datasets) <- gsub("\\..*$", "", dataset.names)

# Read in the data
meta <- readr::read_csv(opt$metadata)

#-------------------------------Dimension Reduction----------------------------#
if (opt$reduce != "none") {
  # Run dimension reduction on each dataset and extract x and y coordinates
  dim.red.data <- lapply(datasets, function(dataset) {
  
    if (opt$reduce == "tsne") {
      # Run tsne
      dim.red <- Rtsne::Rtsne(t(dataset[,-1]), check_duplicates = FALSE)
  
      # Only keep the dimension coordinates
      dim.red <- dim.red$Y
    }
    if (opt$reduce == "pca") {
      # Run tsne
      dim.red <- prcomp(t(dataset[,-1]), scale. = TRUE)
      
      # Only keep the scores for first two PCs
      dim.red <- dim.red$x[, 1:2]
    }
    if (opt$reduce == "umap") {
      # Run tsne
      dim.red <- umap::umap(t(dataset[,-1]))
      
      # Only keep the dimension coordinates
      dim.red <- dim.red$Y
    }
    # Keep the sample names
    names(dim.red) <- samples
    
    # Save these dimensions to a tsv file
    set.name <- names(datasets)[parent.frame()$i[]]
  })
  
}

# Plot with metadata labels
metadata.plots <- lapply(dim.red, function(dataset) {
  
# Get data normalizaion name
set.name <- names(dim.red)[parent.frame()$i[]]

# Make plots for cell type and plate batch 
DimPlot(dataset, metadata, xlabel = paste(opt$reduce, "dim 1"),
        ylabel = paste(opt$reduce, "dim 2"), name = set.name)
})

# Extract the legend
legend <- cowplot::get_legend(metadata.plots[[1]])

# Surpressing the legend in the files
cell.type.plots <- lapply(cell.type.plots, function(a.plot) {
a.plot + theme(legend.position = 'none')
})

# Put all plots and legend together
main.plot <- cowplot::plot_grid(
cowplot::plot_grid(plotlist = metadata.plots, ncol = 2, align='hv'),
cowplot::plot_grid(NULL, legend, ncol = 2),
rel_widths=c(1, 0.2))

# Save to png
ggplot2::ggsave(plot = main.plot, 
                filename = file.path("results", paste0(opt$label, "tsne.png")))

# Print out plot here
main.plot

#Run cell type analyses on all the dim.red.data datasets
# Get knn and kmeans results for all dim.red.data's of all datasets
cell.type.results <- lapply(dim.red.data, function(dataset) {
  # Get clustering results of the data
  knn.results <- KnnEval(dataset, metadata = cell.types)
  kmeans.results <- KmeansEval(dataset, metadata = cell.types)
  
  # Return data frame of combined results
  data.frame(knn.results, kmeans.results)
})


# Run plate batch analyses on all the tsne datasets
# Get knn and kmeans results for all tsne's of all datasets
batch.results <- lapply(tsne, function(dataset) {
  # Get clustering results of the data  
  knn.results <- KnnEval(dataset, metadata = plate.batch)
  kmeans.results <- KmeansEval(dataset,
                               metadata = plate.batch)
  
  # Return data frame of combined results  
  data.frame(knn.results, kmeans.results)
})


# Plot all the cell type statistics on a boxplot
# Melt the  list into one data.frame
cell.type.df <- reshape::melt(cell.type.results)

# Make the plot
celltype.plot <- ggplot(data = cell.type.df, aes(x = L1, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge()) + 
  xlab("Normalization method") +
  facet_wrap(~variable)

# Save plots to png
ggplot2::ggsave(plot = celltype.plot,
                filename = file.path("results", "arnon_cell_type_plots.png"), 
                width = 10)

# Print plot
celltype.plot


Plot all the batch statistics on a boxplot
# Melt the  list into one data.frame
batch.df <- reshape::melt(batch.results)

# Make the plot
batch.plot <- ggplot(data = batch.df, aes(x = L1, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge()) + 
  xlab("Normalization method") +
  facet_wrap(~variable)

# Save plots to png
ggplot2::ggsave(plot = batch.plot,
                filename = file.path("results", "arnon_batch_plots.png"), 
                width = 10)

# Print plot
batch.plot
