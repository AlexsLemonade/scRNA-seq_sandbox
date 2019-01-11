# C.Savonen, CCDL for ALSF
# 2018
# Make functions to use for comparing normalization methods using clustering 

# KmeansEval and KnnEval functions are adapted from Hu and Greene, 2018 CZI scripts
# which can be found here: 
# https://github.com/greenelab/CZI-Latent-Assessment/tree/master/single_cell_analysis

DimPlot <- function(feature, metadata, xlabel = "x dim", ylabel = "y dim",
                    name = ""){
  # This function is to plot data and label it by a given variable
  # Args:
  #  feature: a data.frame with two columns of data with coordinates to be plotted
  #       on x and y axis respectively. This function assumes the first column is
  #       x coordinates and the second column is y coordinates
  #  metadata: vector contains metadata labels
  #  xlabel: what to label the x dimension in the plot
  #  ylabel: what to label the y dimension in the plot
  #  name: what you would like the title of the plot to be
  # Returns:
  #  scatterplot with labels of the metadata provided
  library(ggplot2)
  
  # Make metadata a factor if it is not
  feature <- data.frame("x_data" = feature[, 1],
                        "y_data" = feature[, 2],
                        "metadata" = as.factor(metadata))
  
  # Make a scatter plot with the labels our metadata
  plot.data <- ggplot(feature, aes(x = x_data, y = y_data, group = metadata)) +
               geom_point(aes(colour = metadata)) + 
               colorblindr::scale_fill_OkabeIto() +
               xlab(xlabel) +
               ylab(ylabel) +
               ggtitle(name) 
  
  # Don't print legend if there are more than 10 groups, it's just too 
  # ridiculously big
  if (length(unique(metadata)) > 10) {
    plot.data <- plot.data + theme(legend.position = "none")
  }
  return(plot.data) 
}

KmeansEval <- function(feature, metadata = metadata, n.iter = 10) {
  # This function is used to perform iterative k-means clustering based on projected 
  # features for single cell data and then evaluate the performance according to 
  # Normalized mutual information (NMI) and adjusted rand index (ARI)
  #
  # Args:
  #  feature: a data.frame contains x and y dimensions as the first and second 
  #           columns
  #  metadata: vector contains cell type or other metadata information
  #  iter: number of iterations for k-means clustering
  # Returns:
  #   NMI and ARI results for each iteration, as a dataframe
  #   
  # Convert celltype into numbers
  metadata.num <- as.numeric(factor(metadata))
  sample.id <- seq(1:nrow(feature))
  
  # Objects to store the stats over the iterations
  nmi.score.all <- c()
  ari.score.all <- c()
  
  for(iter in 1:n.iter){
    # Set k equal to the number of celltypes in the dataset
    k <- length(unique(metadata))
    
    # Perform k means clustering
    km <- kmeans(feature, k)
    
    # True clusters
    orignal.data <- data.frame(sample.id, metadata.num)
    
    # Predicted clusters
    cl.data <- data.frame(sample.id, km$cluster)
    
    # Calculate NMI and ARI score
    nmi.score <- NMI::NMI(orignal.data, cl.data)$value
    ari.score <- mclust::adjustedRandIndex(km$cluster, metadata.num)
    
    nmi.score.all <- c(nmi.score.all, nmi.score)
    ari.score.all <- c(ari.score.all, ari.score)
  }
  
  # Compile all results into a data.frame
  results <- data.frame(ari = ari.score.all, nmi = nmi.score.all)
  return(results)
}

KnnEval <- function(feature, metadata = metadata, n.iter = 10){
  # This function performs knn based evaluation for the number of iterations and
  # returns kappa of each iteration
  # Args:
  #  feature: a data.frame contains x and y dimensions
  #  metadata: vector contains cell type or other metadata information  
  #  n.iter: number of iterations for knn clustering
  # Returns:
  #  data.frame with kappa scores for each iteration of KNN
  
  # Make the data into a data.frame:
  feature <- data.frame(feature, "metadata" = metadata)
  
  # Split observations into training and testing groups
  data.split <- caret::createDataPartition(feature$metadata, p = .5,
                                           times = n.iter)
  
  # Create empty objects to store the performance information for each iteration
  kappas <- c()
  
  # Go through this iteration that k times
  for (iter in 1:n.iter) {
    
    # Isolate samples for training the model,
    train <- feature[data.split[[iter]], ]
    
    # Isolate samples for testing the model
    test <-  feature[-data.split[[iter]], ]
    
    # Perform KNN model fitting
    knn.fit <- caret::train(metadata~. , data = train, method = "knn",
                            trControl = caret::trainControl(method = "cv", number = 3),
                            preProcess = c("center", "scale"),
                            tuneLength = 10)

    # Evaluate the model
    knn.pred <- predict(knn.fit, newdata = subset(test, select = -c(metadata)))
    
    # Make a confusion matrix
    confusion <- confusionMatrix(test$metadata, knn.pred)
  
    # Add this iteration's kappa to the list  
    kappas <- c(kappas, confusion$overall[2])
  }
  # Return the list of kappas for each iteration as a data frame
  return(data.frame("kappa" = kappas))
}

