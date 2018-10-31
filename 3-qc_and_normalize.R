# C. Savonen 
# CCDL ALSF 2018 
# Normalization of single-cell RNA-seq data
# 
# 
df <- CreateSeuratObject(raw.data = "", project = "glioblastoma")

AddMetaData(object = "meta", metadata = "meta", col.name = "meta")

# Make a violin plot
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent"), nCol = 3)

# plot the genes
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

data.filtered <- FilterCells(object = df, subset.names = c("nGene", "percent"), 
                             low.thresholds = c(200, -Inf),
                             high.thresholds = c(2500, 0.05))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)