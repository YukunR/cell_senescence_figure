library(PCAtools)
library(ggplot2)
library(ggforce)

colSumNormalize <- function(df) {
  col_sums <- colSums(df)
  # 计算所有列总和的平均值
  average_col_sum <- mean(col_sums)
  # 按列归一化，使每列的和为列总和的平均值
  normalized_df <- sweep(df, 2, col_sums, FUN="/") * average_col_sum
}


gene.data.293t <- read.csv("./res/293t/Norm/gene_data_imputation.csv")
gene.data.lo2 <- read.csv("./res/lo2/gene_data_remove_batch_effect.csv")
gene.data.all <- merge(gene.data.293t, gene.data.lo2, by = "Accession", suffixes = c(".293t", ".lo2"))
row.names(gene.data.all) <- gene.data.all$Accession
gene.data.all <- subset(gene.data.all, select = -Accession)

gene.data.norm <- colSumNormalize(gene.data.all)


sample.info <- read.delim("./data/sample_info_all.txt")
gene.dat <- na.omit(gene.data.norm)


calculate_cv <- function(x) {
  sd(x) / mean(x) * 100
}
# 初始化一个逻辑向量，用于标记保留的行
keep_rows <- rep(TRUE, nrow(gene.dat))
# 遍历每个组，计算CV并更新保留的行
for (group in unique(sample.info$Group)) {
  # 构造该组内所有列的名称
  columns <- sample.info$Sample[sample.info$Group == group]
  
  # 计算该组的CV
  group_cv <- apply(gene.dat[, columns], 1, calculate_cv)
  
  # 更新保留的行
  keep_rows <- keep_rows & (group_cv < 40)
}
# 筛选数据
gene.dat <- gene.dat[keep_rows, ]


# remove batch effect
library(limma)
library(sva)
library(ggalt)
batch.info <- c(rep("batch1", 18), rep("batch2", 18))
gene.dat <- log10(gene.dat)
gene.dat <- removeBatchEffect(gene.dat, batch = batch.info)
gene.dat <- as.data.frame(10^gene.dat)
row.names(sample.info) <- sample.info$Sample
sample.info <- subset(sample.info, select = -Sample)


pca.dat <- pca(gene.dat, metadata = sample.info)
screeplot <- screeplot(pca.dat)


# PCA
color <- c("0 h" = "#fff99e", "3 h" = "#bfff7f", "6 h" = "#ffbf7f", 
           "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
pca.dat$metadata$shapeGroup <- c(rep("293t", 18), rep("lo2", 18))
pca.dat$metadata$colorGroup <- c(rep("0 h", 3), rep("3 h", 3), rep("6 h", 3), 
                                 rep("9 h", 3), rep("12 h", 3), rep("24 h", 3), 
                                 rep("0 h", 3), rep("3 h", 3), rep("6 h", 3), 
                                 rep("9 h", 3), rep("12 h", 3), rep("24 h", 3))
pca <- biplot(pca.dat, colby = "colorGroup", shape = "shapeGroup", encircle = T) + 
  # geom_mark_ellipse(aes(fill = pca.dat$metadata$colorGroup,
  #                       color = pca.dat$metadata$colorGroup),
  #                   show.legend=F) +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) + 
  theme(legend.position = 'bottom')
  
pca
ggsave("./res/pca.pdf", pca, width = 9, height = 7)

