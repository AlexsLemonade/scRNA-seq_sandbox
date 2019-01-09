# C.Savonen, CCDL for ALSF
# 2018
# Make functions to use for comparing normalization methods using clustering 

ClusterPlot <- function(dat, metadata, name = "name"){
  # This function is to plot cluster data and label it by a given variable
  # Args:
  #  dat: a data.frame with two columns of data
  #  var: vector contains metadata labels
  #  name: a character string for the plot title and for to save as png
  # Returns:
  #   png and plot in Rmd of the provided data with labels of the metadata provided
  # convert celltype into numbers
  colz <- colors(distinct = TRUE)[runif(length(levels(metadata)), min = 1,
                                        max = length(colors(distinct = TRUE)))]
  plot(dat, pch = 21, bg = colz[metadata], main = name, ylab = "tsne dim 2",
       xlab = "tsne dim 1");
  legend(x = "bottomleft", legend = levels(metadata), fill = colz, cex = 0.6)
}

KmeansEval <- function(feature, metadata = metadata, iter = 10) {
  # This function is used to perform iterative k-means clustering based on projected 
  # features for single cell data and then evaluate the performance according to 
  # Normalized mutual information (NMI) and adjusted rand index (ARI)
  #
  # Args:
  #  feature: a data.frame contains clustering dimensions
  #  metadata: vector contains cell type or other metadata information
  #  iter: number of interation for k-means clustering
  # Returns:
  #   NMI and ARI results for each iteration, as a dataframe
  #   
  # Convert celltype into numbers
  metadata_num <- as.numeric(factor(metadata))
  sample_id <- seq(1:nrow(feature))
  
  # iterative k-means
  nmi_score_all <- c()
  ari_score_all <- c()
  all_cluster <- list()
  
  for(i in 1:iter){
    # set k equal to the number of celltypes in the dataset
    k <- length(unique(metadata))
    
    # perform k means clustering
    km <- kmeans(feature, k)
    
    # true clusters
    orignal_data <- data.frame(sample_id, metadata_num)
    
    # predicted clusters
    cl_data <- data.frame(sample_id, km$cluster)
    
    # calculate NMI and ARI score
    nmi_score <- NMI::NMI(orignal_data, cl_data)$value
    ari_score <- mclust::adjustedRandIndex(km$cluster, metadata_num)
    
    nmi_score_all <- c(nmi_score_all, nmi_score)
    ari_score_all <- c(ari_score_all, ari_score)
  }
  
  # Compile all results into a data.frame
  results <- data.frame(ari = ari_score_all, nmi = nmi_score_all)
  return(results)
}

KnnEval <- function(feature, metadata = metadata, iter = 10){
  # This function performs knn based evaluation 
  # Args:
  #  feature: a data.frame contains clustering dimensions
  #  metadata: vector contains cell type or other metadata information  
  #  iter: number of interation for knn clustering
  # Returns:
  #  list of accuracy scores for each iteration of KNN, as a dataframe
  
  # Make the data into a data.frame:
  feature <- data.frame("tsne" = feature, "metadata" = metadata)
  
  # Split observations into groups
  cv <- cvTools::cvFolds(nrow(feature), K = iter, R = 1)
  
  # Create empty objects to store the performance information for each iteration
  perf.eval <- list()
  confusion.matrix <- 0
  
  # Go through this iteration that k times
  for (i in 1:iter) {
    # Isolate samples for training the model
    train <- feature[cv$subsets[-which(cv$which == i)], ]
    
    # Isolate samples for testing the model
    test <- feature[cv$subsets[which(cv$which == i)], ]
    
    # Perform KNN model fitting
    knn.fit <- caret::train(metadata~. , data = train, method = "knn",
                            trControl = caret::trainControl(method = "cv", number = 3),
                            preProcess = c("center", "scale"),
                            tuneLength = 10)
    
    # Evaluate the model
    knn.pred <- predict(knn.fit, newdata = subset(test, select = -c(metadata)))
    perf.eval[[i]] <- round(cal_performance(knn.pred, test$metadata, 3), 2)
    
    # Make the results into a matrix
    matrix <- as.matrix(table(test$metadata, knn.pred, deparse.level = 0))
    confusion.matrix <- matrix + confusion.matrix
  }
  
  # Get mean performance of cross validation
  perf.eval <- dplyr::bind_rows(perf.eval)
  accuracy <- perf.eval$accuracy
  
  return(data.frame("knn" = accuracy))
}


PlotResults <- function(results.list, name = "results") {
  # This function makes a boxplot of the cluster statistic results 
  # Args:
  #  results.list: a list of dataframes which each column contains a different 
  #                statistic for 10 iterations
  #  name: name to use for the png to be saved and the plot title 
  # Returns: boxplots of the normalization method cluster statistics. Prints the
  #          plots and save the plot as a png
  
  # Get meta info
  meta <- names(unlist(results.list))
  
  # Transform list into a dataframe:
  ggplot.df <- data.frame("method" = stringr::word(meta, 1, sep = "\\."),
                          "test" = gsub("[0-9]*$", "", 
                                        stringr::word(meta, 2, sep = "\\.")),
                          "iter" = rep(1:nrow(results.list[[1]]),
                                       length(results.list)*ncol(results.list[[1]])),
                          "values" = unlist(results.list))
  
  # Make the plot
  plot <- ggplot(data = ggplot.df, aes(x = method, y = values, fill = test)) +
    geom_boxplot(position = position_dodge()) + 
    xlab("Normalization method") +
    ggtitle(name) +
    facet_wrap(~test)
  
  # Save plot to png
  ggsave(paste0(name, "_cluster_results.png"), width = 10) 
  
  # Print plot
  plot
}

