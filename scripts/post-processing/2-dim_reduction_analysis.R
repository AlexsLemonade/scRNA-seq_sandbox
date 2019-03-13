# C.Savonen, CCDL for ALSF
# 2019
# 
# Purpose: Perform dimension reduction analyses of scRNA-seq gene matrix .tsv 
#          
# "-d" : Path to gene matrix in tab delimited format, gene x sample with gene 
#        info as the first column with gene info column labeled with the word 
#        'gene' in it
# "-m" : Path to metadata file that contains only the metadata variables you 
#        wish to test and label by. First column must contain the sample names
#        that are in the datasets
# "-r" : Dimension reduction technique to use. Options are: 'pca', 'tsne', or 
#        'umap'. Default is pca."
# "-o" : Directory where you would like the output to go. Default is current 
#        directory
# "-l" : Optional label for output files. 

# Command line example:
#
# Rscript scripts/6-dim_reduction_analysis.R \
# -d data/gene_matrix.tsv \
# -m metadata.tsv \
# -r pca \
# -o results \
# -l darmanis
# 
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Import functions for this analysis
source(file.path("scripts", "util", "clustering_statistics_functions.R"))

# ggplot2 library
library(ggplot2)
library(optparse)

#--------------------------------Set up options--------------------------------#
option_list <- list(
  make_option(opt_str = c("-d", "--data"), type = "character", default = getwd(),
              help = "Path to gene matrix in tab delimited format, gene x sample
              with gene info column labeled with the word 'gene' in it",
              metavar = "character"),
  make_option(opt_str = c("-m", "--metadata"), type = "character",
              default = "none", help = "Path to metadata file that contains 
              only the metadata variables you wish to test and label by. First
              column must contain the sample names that are in the datasets",
              metavar = "character"),
  make_option(opt_str = c("-r", "--reduce"), type = "character",
              default = "pca", help = "Dimension reduction technique to use.
              Options are: 'pca', 'tsne', or 'umap'. Default is pca.",
              metavar = "character"),
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
# Check that the dimension reduction option given is supported
if (!(opt$reduce %in% c("tnse", "pca", "umap"))){
  stop("That is not a dimension reduction technique supported by this script. 
       Check for typos. Acceptable options: 'tsne', 'pca', or 'umap'")
}

# Create the output for results folder if it does not exist
if (!dir.exists(opt$output)) {
  message(paste("Output folder:", opt$output, "does not exist, creating one"))
  dir.create(opt$output)
}

# Append an underscore to the label if it is given
if (opt$label != "") {
  opt$label <- paste0("_", opt$label)
}

# Add this warning to prevent the whole thing running and then not working when 
# it gets to the metadata step
if (!file.exists(opt$metadata)) {
  warning("Metadata file not found. Check the path given for option -m.")
}
#---------------------------------Read in data---------------------------------
# Get the file names of all the normalized files
dataset.files <- dir(opt$data, full.names = TRUE)

# Read in each of the normalization files
datasets <- lapply(dataset.files, readr::read_tsv, guess_max = 10000)

# Get names of datasets
dataset.names <- dir(opt$data)

# Remove file extension name and make it the name in the list
names(datasets) <- gsub("\\..*$", "", dataset.names)

# Read in the metadata
meta <- readr::read_tsv(opt$metadata)
colnames(meta)[1] <- "sample_id"

# Filter the metadata if needed
meta <- meta %>% dplyr::filter(sample_id %in% colnames(datasets[[1]]))

#-------------------------------Dimension Reduction----------------------------#
# Run dimension reduction on each gene expression dataset and extract x and y 
# coordinates
# Although a typical gene expression matrix consists of rows being genes and
# columns being samples, these methods generally assume columns to be the different
# features, so the data are transposed in order to obtain values for each sample. 
dim.red.data <- lapply(datasets, function(dataset) {
  
  # Find if there's a gene column so we can get rid of it
  gene.col <- grep("gene", colnames(dataset), ignore.case = TRUE)
  
  if (length(gene.col > 0 )) {
    # Make dataset the data only
    dataset <- dataset[, -gene.col]
  }
  
  # Check that the metadata and data are reasonably compatible
  if (ncol(dataset) != nrow(meta)) {
    stop("Metadata length doesn't match dataset's number of columns.
         Make sure the metadata that you are providing matches the order 
         and length of the dataset")
  }
  
  # Extract sample names
  samples <- colnames(dataset)
  
  # Run the dimension reduction technique
  if (opt$reduce == "tsne") {
    # Run tsne
    dim.red <- Rtsne::Rtsne(t(dataset), check_duplicates = FALSE)
  
    # Only keep the dimension coordinates
    dim.red <- data.frame(dim.red$Y)
  }
  if (opt$reduce == "pca") {
    # Run pca
    dim.red <- prcomp(t(dataset))
      
    # Only keep the scores for first two PCs
    dim.red <- data.frame(dim.red$x[, 1:2])
  }
  if (opt$reduce == "umap") {
    # Run umap
    dim.red <- umap::umap(t(dataset))
      
    # Only keep the dimension coordinates
    dim.red <- data.frame(dim.red$Y)
  }
  # Keep the sample names
  rownames(dim.red) <- samples
    
  # Extract this dataset's name 
  set.name <- names(datasets)[parent.frame()$i[]]
    
  # Save these dimensions to a tsv file with their dataset name
  readr::write_tsv(dim.red, 
                   file.path(opt$output, paste0(opt$reduce, "_", set.name,
                                                ".tsv")))
})
# saveRDS(dim.red.data, file.path(opt$output,
#                                 paste0(opt$reduce, "_", set.name,".RDS"))
                                 
#-----------------------Plot with metadata variable labels---------------------#
# Format metadata for easier use
meta <- meta %>% 
  dplyr::mutate_all(as.factor) %>% 
  as.list()

# Obtain variable names from metadata import
variable.names <- gsub(".ch1", "", names(meta))
variable.names <- gsub("\\.", "_", variable.names)

for (variable in 1:length(meta)) {
  # Start with legend as NULL
  legend <- NULL
  
  # Plot with metadata labels
  metadata.plots <- lapply(dim.red.data, function(dataset) {
  
    # Get data normalizaion name
    set.name <- names(dim.red.data)[parent.frame()$i[]]
    
    # Make plots for cell type and plate batch 
    DimPlot(dataset, meta[[variable]], xlabel = paste(opt$reduce, "dim 1"),
            ylabel = paste(opt$reduce, "dim 2"), 
            name = paste0(variable.names[variable], "_", set.name))
  })
  
  # Try to obtain legend if it exists
  try(legend <- cowplot::get_legend(metadata.plots[[1]]), silent = TRUE) 

  # Surpressing the legend in the files
  metadata.plots <- lapply(metadata.plots, function(a.plot) {
                           a.plot + theme(legend.position = 'none')
                          })

  # Put all plots and legend together
  main.plot <- cowplot::plot_grid(cowplot::plot_grid(plotlist = metadata.plots,
                                                     ncol = 2, align = 'hv'),
                                  cowplot::plot_grid(NULL, legend, ncol = 2),
                                  rel_widths = c(1, 0.2))

  # Save to png
  ggplot2::ggsave(plot = main.plot, 
                  file.path(opt$output, paste0(opt$reduce, opt$label, "_",
                                               variable.names[variable], ".png")))
}
