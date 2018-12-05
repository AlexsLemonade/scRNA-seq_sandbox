# C.Savonen 
# CCDL for ALSF 2018
# 
# Run tSNE and PCA for data, look at where different sample types end up clustered

if (!("tsne" %in% installed.packages())) {
  install.packages("tsne")
}

# Import gene expression matrix data
salmon.data <- readRDS(file.path("data", "salmon.data.RDS"))
salmon.data <- salmon.data[ ,-c(1,2)]
salmon.data <- salmon.data[, which(salmon.prop.assigned > 0.50)]

# Get sample info for batch effect testing
geo.meta <- read.csv(file.path("data", "meta_data.csv"))
geo.meta <- geo.meta[match(colnames(salmon.data),geo.meta$geo_accession),]
geo.meta.levels <- apply(geo.meta, 2, function(x) length(levels(as.factor(x))))
geo.meta <- as.list(geo.meta[,-which(geo.meta.levels==nrow(geo.meta)|geo.meta.levels==1)])

# Do PCA of samples
pca <- prcomp(t(salmon.data))

# Do tSNE of samples
tsne.dat <- tsne::tsne(t(salmon.data))

# Do linear regression of PC scores with the various variables we have
results <- lapply(geo.meta, function(y) {
  -log10(summary(lm(pca$x[,1]~as.factor(y)))$coefficients[2,4])
  })
results <- unlist(results)[-grep("characteristics", names(geo.meta))]
names(results) <- gsub("\\.ch1", "", names(results))

# Make a plot of the regression results
results <- results[c(5, 9, 2:4, 7:8)]
png(file.path("results", "PC_score_regression.png"))
par(mar = c(10,3,3,3))
barplot(unlist(results), las = 2, main = "-log10(pval) for PCA Score Regression")
dev.off()

# Make PCA plot function that allows us 
pc.plot <- function(dat, var){
  colz <- colors(distinct = TRUE)[runif(length(levels(var)), min = 1,
                                        max = length(colors(distinct=TRUE)))]
  plot(dat, pch = 21, bg = colz[var]);
  legend(x = "topleft", legend = levels(var), fill = colz, cex = 0.8)
}

#----------------------------------PCA Plots-----------------------------------#
png(file.path( "results", paste0("Plate_ID_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$plate.id.ch1))
dev.off()

png(file.path( "results", paste0("Well_ID_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$characteristics_ch1.2))
dev.off()

png(file.path( "results", paste0("Patient_ID_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$patient.id.ch1))
dev.off()

png(file.path( "results", paste0("Tissue_Type_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$characteristics_ch1.3))
dev.off()

png(file.path( "results", paste0("Cell_Type_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$characteristics_ch1.6))
dev.off()

png(file.path( "results", paste0("tSNE_Cluster_PCA.png")))
pc.plot(pca$x, as.factor(geo.meta$tsne.cluster.ch1))
dev.off()

#----------------------------------tSNE Plots----------------------------------#
png(file.path( "results", paste0("Patient_ID_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$patient.id.ch1))
dev.off()

png(file.path( "results", paste0("Well_ID_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$characteristics_ch1.2))
dev.off()

png(file.path( "results", paste0("Author_tSNE_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$tsne.cluster.ch1))
dev.off()

png(file.path( "results", paste0("Plate_ID_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$plate.id.ch1))
dev.off()

png(file.path( "results", paste0("Cell_Type_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$characteristics_ch1.6))
dev.off()

png(file.path( "results", paste0("Tissue_Type_tSNE.png")))
pc.plot(tsne.dat, as.factor(geo.meta$characteristics_ch1.3))
dev.off()
