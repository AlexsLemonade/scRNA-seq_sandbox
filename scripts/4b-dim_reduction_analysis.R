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
              default = "none", help = "Path to metadata file that contains 
              only the metadata variables you wish to test and label by.",
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

opt$metadata <- file.path("darmanis_data", "GSE84465_meta.tsv")
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

#-----------------------Plot with metadata variable labels---------------------#
# Read in the metadata
meta <- readr::read_tsv(opt$metadata)

# Obtain variable names from metadata import
variables <- names(metadata)

for (variable in ncol(metadata)) {
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
  main.plot <- cowplot::plot_grid(cowplot::plot_grid(plotlist = metadata.plots,
                                                    ncol = 2, align='hv'),
                                  cowplot::plot_grid(NULL, legend, ncol = 2),
                                  rel_widths=c(1, 0.2))

  # Save to png
  ggplot2::ggsave(plot = main.plot, 
                  filename = file.path("results", paste0(opt$label, opt$reduce,
                                                         ".png")))

  # Print out plot here
  main.plot
}