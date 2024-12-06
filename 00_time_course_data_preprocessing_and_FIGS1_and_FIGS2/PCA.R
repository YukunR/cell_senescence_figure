library(PCAtools)
library(corrplot)
library(pheatmap)
library(ggforce)
library(tidyverse)

PCA <- function(gene.dat, sample.info, output.dir) {
  # >>> data preparation >>>
  row.names(gene.dat) <- gene.dat$Accession
  gene.dat <- subset(gene.dat, select = -Accession)
  row.names(sample.info) <- sample.info$Sample
  sample.info <- subset(sample.info, select = -Sample)
  
  # >>> correlation >>>
  cor <- cor(gene.dat)
  corplot <- pheatmap(cor)
  pdf(paste0(output.dir, "corplot.pdf"))
  corplot
  dev.off()
  
  # dendrogram
  dist <- dist(t(gene.dat))
  hc <- hclust(dist)
  pdf(file = paste0(output.dir, "dendrogram.pdf"))
  plot(hc)
  dev.off()
  
  # scree plot
  pca.dat <- pca(gene.dat, metadata = sample.info)
  screeplot <- screeplot(pca.dat)
  ggsave(paste0(output.dir, "screeplot.pdf"), screeplot)
  
  # PCA
  pca <- biplot(pca.dat, colby = "Group") + 
    geom_mark_ellipse(aes(fill = sample.info$Group,
                          color = sample.info$Group), 
                      show.legend=F) +
    theme(legend.position = 'bottom') + 
    coord_equal()
  ggsave(paste0(output.dir, "pca.pdf"), pca, width = 10, height = 7)
  
  # loading
  loading <- plotloadings(pca.dat)
  ggsave(paste0(output.dir, "loading.pdf"), loading, width = 7, height = 5)
}