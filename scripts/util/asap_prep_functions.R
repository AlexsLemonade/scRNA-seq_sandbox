# C.Savonen, CCDL for ALSF
# 2018
# 
# These functions are used for prepping tximport counts data into a filtered 
# dataset that is ready for uploading to [ASAP online](https://asap.epfl.ch/)

AsapFilter <- function(data, min_counts = 1, perc_genes = 0.01, num_genes = 100) {
  # This function is filters and makes data into ASAP format and assumes gene info 
  # is the first column
  # Args:
  #  data: a gene expression data.frame that is gene x samples, with the first column
  #        containing the gene names
  #  min_counts: the cutoff for minimum number of counts for a gene to be 
  #              considered expressed in a particular sample
  #  perc_genes: The minimum percent of samples that need to express that gene 
  #              to keep the gene in the set
  #  num_genes: The minimum number of genes a particular sample must express to 
  #             be kept in the gene set
  # Returns:
  #   A filtered gene expression data.frame, that is able to be submitted to ASAP
  # Get rid of decimal points
  # 
  # Store the genes separately
  gene <- data[, 1]
  
  # Get rid of decimals (even if they are .000, ASAP doesn't like them)
  data <- apply(data[, -1], 2, round)

  # Filter genes that are expressed in 1% of cells
  gene.sum <- apply(data >= min_counts, 1, sum)
  data <- data[which(gene.sum > ncol(data)*perc_genes), ]

  # Filter samples that express at least 100 genes
  sample.sum <- apply(data >= min_counts, 2, sum)
  data <- data[, which(sample.sum > num_genes) ]

  # Need the genes to be its own column
  data <- data.frame(gene, data)
  
  return(data)
}
