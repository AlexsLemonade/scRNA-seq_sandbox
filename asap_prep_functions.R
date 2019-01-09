# C.Savonen, CCDL for ALSF
# 2018
# 
# These functions are used for prepping tximport counts data into a filtered 
# dataset that is ready for uploading to [ASAP online](https://asap.epfl.ch/)

AsapFilter <- function(data, genes, min_counts = 1, perc_genes = 0.01, num_genes = 100) {
  # This function is filters and makes data into ASAP format
  # Args:
  #  data: a gene expression data.frame that is gene x samples
  #  genes: a character vector containing the gene names (corresponding to the
  #         rows in data)
  #  min_counts: the cutoff for minimum number of counts for a gene to be 
  #              considered expressed in a particular sample
  #  perc_genes: The minimum percent of samples that need to express that gene 
  #              to keep the gene in the set
  #  num_genes: The minimum number of genes a particular sample must express to 
  #             be kept in the gene set
  # Returns:
  #   A filtered gene expression data.frame, that is able to be submitted to ASAP
  # Get rid of decimal points
  data <- apply(data[,-1], 2, round)
  rownames(data) <- genes

  # Identify which cells have at least 1 count
  boolean <- data >= min_counts

  # Filter genes that are expressed in 1% of cells
  gene.sum <- apply(boolean, 1, sum)
  data <- data[which(gene.sum > ncol(data)*perc_genes), ]

  # Identify which cells have at least 1 count
  boolean <- data >= min_counts

  # Filter samples that express at least 100 genes
  sample.sum <- apply(boolean, 2, sum)
  data <- data[, which(sample.sum > num_genes) ]

  # Need the genes to be its own column
  data <- data.frame("gene" = rownames(data), data)
  
  return(data)
}
