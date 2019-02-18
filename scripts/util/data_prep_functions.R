# C.Savonen, CCDL for ALSF
# 2018
#
# These functions are used for prepping tximport counts data into a filtered
# dataset that is ready for uploading to [ASAP online](https://asap.epfl.ch/)

GeneMatrixFilter <- function(dataset, min_counts = 1, perc_genes = 0.01, num_genes = 100) {
    # This function is filters and makes dataset into ASAP format and assumes gene info
    # is the first column
    # Args:
    #  dataset: a gene expression data.frame that is gene x samples, with the first column
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
    gene <- dataset[, 1]
    
    # Get rid of decimals (even if they are .000, ASAP doesn't like them)
    dataset <- dataset[, -1] %>% apply(., 2, round)
    
    # Find genes that are expressed in 1% of cells
    gene.sum <- apply(dataset >= min_counts, 1, sum)
    perc.genes <- which(gene.sum > ncol(data)*perc_genes)
    
    # Only keep those genes
    dataset <- dataset[perc.genes, ]
    gene <- gene[perc.genes, ]
    
    # Filter samples that express at least 100 genes
    sample.sum <- apply(dataset >= min_counts, 2, sum)
    dataset <- dataset[, which(sample.sum > num_genes) ]
    
    # Need the genes to be its own column
    dataset <- data.frame(gene, dataset)
    
    return(dataset)
}
