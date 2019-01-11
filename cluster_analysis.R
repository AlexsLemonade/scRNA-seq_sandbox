# C.Savonen, CCDL for ALSF
# 2019
# 
### Purpose: compare normalization methods

# Import functions for this analysis
source(file.path("scripts", "util", "clustering_statistics_functions.R"))

# ggplot2 library
library(ggplot2)

# These are the necessary packages
packages <- c("Rtsne", "NMI", "caret")

# Check if these packages are installed and install them if they aren't
lapply(packages, function(package) {
  if (!(package %in% installed.packages())) {
    install.packages(package)
  }
})

# Get the file names of all the normalized files
normalized.files <- dir(file.path("arnon_data", "normalized_arnon"))

# Read in each of the normalization files
normalized.data <- lapply(normalized.files, function(file) {
  readr::read_tsv(file.path("arnon_data",
                            "normalized_arnon",
                            file))
})

# Keep all the gene lists
genes <- lapply(normalized.data, function(x) x[,1])

# Keep the names with it.
names(normalized.data) <- gsub("\\.tab|arnon_", "", normalized.files)

# Read in the data
meta <- readr::read_csv(file.path("arnon_data", "meta_data.csv"))

# Get sample names from columns
samples <- colnames(normalized.data[[1]])[-1]

# Keep metadata only for the samples we have
meta <- meta[match(samples, meta$title), ]

# Cell types info as it's own vector
cell.types <- as.factor(meta$cell.type.ch1)

# Plate batch info as a factor
plate.batch <- as.factor(meta$cohort.ch1)

# Run tsne on each dataset and extract
tsne <- lapply(normalized.data, function(dataset) {
  # Run tsne
  tsne.res <- Rtsne::Rtsne(t(dataset[,-1]), check_duplicates = FALSE)
  
  # Only keep the dimension coordinates
  tsne.res <- tsne.res$Y
  
  # Keep the sample names
  names(tsne.res) <- samples
})

# Plot with cell type labels
cell.type.plots <- lapply(tsne, function(dataset) {
# Get data normalizaion name
set.name <- names(tsne)[parent.frame()$i[]]

# Make plots for cell type and plate batch 
DimPlot(dataset, cell.types, xlabel = "tsne dim 1",
ylabel = "tsne dim 2", name = set.name)
})

# Extract the legend
legend <- cowplot::get_legend(cell.type.plots[[1]])

# Surpressing the legend in the files
cell.type.plots <- lapply(cell.type.plots, function(a.plot) {
a.plot + theme(legend.position = 'none')
})

# Put all plots and legend together
main.plot <- cowplot::plot_grid(
cowplot::plot_grid(plotlist = cell.type.plots, ncol = 2, align='hv'),
cowplot::plot_grid(NULL, legend, ncol = 2),
rel_widths=c(1, 0.2))

# Save to png
ggplot2::ggsave(plot = main.plot, 
filename = file.path("results", "arnon_cell_type_tsne.png"))

# Print out plot here
main.plot

# Plot batches
batch.plots <- lapply(tsne, function(dataset) {
# Get data normalizaion name
set.name <- names(tsne)[parent.frame()$i[]]

# Make plots for cell type and plate batch 
DimPlot(dataset, plate.batch, xlabel = "tsne dim 1",
ylabel = "tsne dim 2", name = set.name)
})

# Extract the legend
legend <- cowplot::get_legend(cell.type.plots[[1]])

# Surpressing the legend in the files
batch.plots <- lapply(batch.plots, function(a.plot) {
a.plot + theme(legend.position = 'none')
})

# Put all plots and legend together
main.plot <- cowplot::plot_grid(
cowplot::plot_grid(plotlist = batch.plots, ncol = 2, align='hv'),
cowplot::plot_grid(NULL, legend, ncol = 2),
rel_widths=c(1, 0.2))

# Save to png
ggplot2::ggsave(plot = main.plot, 
filename = file.path("results", "arnon_batch_tsne.png"))

# Print out plot here
main.plot

#Run cell type analyses on all the tsne datasets
# Get knn and kmeans results for all tsne's of all datasets
cell.type.results <- lapply(tsne, function(dataset) {
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
