# C.Savonen, CCDL for ALSF
# 2019
# 
# Purpose: Analyze clustering of metadata scRNA-seq gene matrix .tsv using 
#          kmeans and KNN clustering
#          
# "-d" : Path to gene matrix in tab delimited format, gene x sample with gene 
#        info labeled with the word "gene" as the column name. Default is to 
#        read all files in given directory with ".tsv" into the dataset. If you
#        prefer to explicitly list the files to be read in, separate file paths/names
#        with a space. eg "pca_dataset_1.tsv pca_dataset_2.tsv"
# "-m" : Path to metadata file that contains only the metadata variables you 
#        wish to test and label by. Should be in same order as the samples. 
#        Generally expected to be variables cleaned up from GEO metadata, but
#        other types of metadata in a tsv also work. This script will assume 
#        the metadata comes with column names.  
# "-o" : Directory where you would like the output to go. Default is current 
#        directory
# "-l" : Optional label for output files. 

# Command line example:
#
# Rscript scripts/7-cluster_analysis.R \
# -d pca_data \
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
              help = "Directory path to folder of dataset(s) of gene matrices 
              (probably dimension reduced) in tsv format, gene x sample to test.",
              metavar = "character"),
  make_option(opt_str = c("-m", "--metadata"), type = "character",
              default = "none", help = "Path to metadata file that contains 
              only the metadata variables you wish to test and label by. Typically
              expected to be metadata from GEO but other categorizations of data
              in a tsv format, in the same order as the samples work.",
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

# Create the output for results folder if it does not exist
if (!dir.exists(opt$output)) {
  message(paste("Output folder:", opt$output, "does not exist, creating one"))
  dir.create(opt$output)
}

# Append an underscore to the label if it is given
if (opt$label != "") {
  opt$label <- paste0("_", opt$label)
}

#---------------------------------Read in data---------------------------------#
# Obtain list of file names: 
if (dir.exists(opt$data)) { # If directory is given, read all tsv files
  # Get the file names of all the normalized files
  dataset.files <- dir(opt$data, full.names = TRUE)

  # Only read in tsv files
  dataset.files <- grep("\\.tsv$", dataset.files, value = TRUE)
  
} else { # Otherwise, make file list by separating by spaces
  dataset.files <- opt$data <- strsplit(opt$data, " ")
}

# Read in each of the normalization files
datasets <- lapply(dataset.files, readr::read_tsv)

# Obtain names of datasets
dataset.names <- stringr::word(dataset.files, sep = "/", -1)

# Remove file extension name and make it the name in the list
names(datasets) <- gsub("\\..*$", "", dataset.names)

#-------------------------------Set up metadata--------------------------------#
# Read in the metadata
# This script is generally built on these metadata being from GEO metadata
# files. 
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
    message(paste0("Groups for variable: '", variable.name, 
                   "' are too irregular to run KNN clustering only kmeans will 
                   be done"))
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
  cluster.results.df <- reshape2::melt(cluster.results)
  
  # Put more sensible colnames
  colnames(cluster.results.df) <- c("stat", "value", "dataset")
  
  # Save these data to a tsv 
  readr::write_tsv(cluster.results.df, 
                   file.path(opt$output, paste0("cluster_stats_", variable.name,
                                                opt$label, ".tsv")))
  
  # Make the plot
  cluster.stats.plot <- ggplot(data = cluster.results.df, 
                               aes(x = dataset, y = value, fill = stat)) +
    geom_boxplot(position = position_dodge()) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Normalization method") +
    facet_wrap(~stat)

  # Save plots to png
  ggplot2::ggsave(plot = cluster.stats.plot,
                  file.path(opt$output, paste0(variable.name, opt$label,
                                               "_cluster_stats_barplots.png")), 
                  width = 10)
                
})
