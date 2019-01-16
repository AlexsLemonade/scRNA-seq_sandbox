# C.Savonen, CCDL for ALSF
# 2019
# 
# Purpose: Analyze clustering of metadata scRNA-seq gene matrix .tsv using 
#          kmeans and KNN clustering
#          
# "-d" : Path to gene matrix in tab delimited format, gene x sample with gene 
#        info as the first column
# "-m" : Path to metadata file that contains only the metadata variables you 
#        wish to test and label by.
# "-o" : Directory where you would like the output to go. Default is current 
#        directory
# "-l" : Optional label for output files. 

# Command line example:
#
# Rscript scripts/4b-dim_reduction_analysis.R \
# -d data/gene_matrix.tsv \
# -m metadata.tsv \
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
              help = "Directory path to folder of dataset(s) of gene matrices in
              tsv format, gene x sample to test.",
              metavar = "character"),
  make_option(opt_str = c("-m", "--metadata"), type = "character",
              default = "none", help = "Path to metadata file that contains 
              only the metadata variables you wish to test and label by.",
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

#---------------------------------Read in data---------------------------------#
# Get the file names of all the normalized files
dataset.files <- dir(opt$data, full.names = TRUE)

# Read in each of the normalization files
datasets <- lapply(dataset.files, readr::read_tsv)

# Get names of datasets
dataset.names <- dir(opt$data)

# Remove file extension name and make it the name in the list
names(datasets) <- gsub("\\..*$", "", dataset.names)

#-------------------------------Set up metadata--------------------------------#
# Read in the metadata
meta <- as.list(readr::read_tsv(opt$metadata))

# Make all the metadata into factors
meta <- lapply(meta, as.factor)

# Obtain variable names from metadata import
variable.names <- gsub(".ch1", "", names(meta))
variable.names <- gsub("\\.", "_", variable.names)

#-------------------------------Run cluster analyses---------------------------#
# For each metadata variable, run the clustering and plot it
lapply(meta, function(meta.var) {
  
  # Extract the variable name 
  variable.name <- variable.names[parent.frame()$i[]]
  
  # Check how large the smallest group in the metadata is
  no.knn <- min(summary(meta.var)) < 4
  if (no.knn) {
    message("Groups are too irregular to run KNN clustering only kmeans will be
            done")
  }

  # Get clustering results for all datasets
  cluster.results <- lapply(datasets, function(dataset) {
                            # Run kmeans clustering 
                            kmeans.results <- KmeansEval(dataset, metadata = meta.var)
    
                            # Run KNN if groupings are large enough
                            if (no.knn) {
                              return(kmeans.results)
                            } else {
                              knn.results <- KnnEval(dataset, metadata = meta.var)
                              
                              # Return data frame of combined results
                              return(data.frame(knn.results, kmeans.results))
                            }
  })

  # Plot cluster statistics on a boxplot
  # Melt the  list into one data.frame
  cluster.results.df <- reshape::melt(cluster.results)

  # Make the plot
  cluster.stats.plot <- ggplot(data = cluster.results.df, 
                               aes(x = L1, y = value, fill = variable)) +
    geom_boxplot(position = position_dodge()) + 
    xlab("Normalization method") +
    facet_wrap(~variable)

  # Save plots to png
  ggplot2::ggsave(plot = cluster.stats.plot,
                  file.path(opt$output, paste0(opt$label, "_", variable.name,
                                               "_plots.png")), 
                  width = 10)
                
})
