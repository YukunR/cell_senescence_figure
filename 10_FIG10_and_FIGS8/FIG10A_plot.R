library(dplyr)


# WGCNA marker 标准：>0.8的相关性或最相关的正负相关模块的softConnectivity前5%蛋白
# mfuzz marker 标准：lo2 2/3/7 cluster和293t1 1/4 cluster蛋白且ANOVAp < 0.05的交集
# marker：WGCNA和mfuzz marker的交集


mfuzz.lo2 <- read.csv("./res/lo2/mFuzz/gene_dat_cluster.csv")
mfuzz.lo2 <- mfuzz.lo2 %>% filter(ANOVA_p < 0.05 & cluster %in% c(2, 3, 7))

mfuzz.293t <- read.csv("./res/293t/mFuzz/gene_dat_cluster.csv")
mfuzz.293t <- mfuzz.293t %>% filter(ANOVA_p < 0.05 & cluster %in% c(1, 4))

mfuzz.marker <- merge(data.frame(Accession = mfuzz.lo2$Accession), data.frame(Accession = mfuzz.293t$Accession))

wgcna.marker <- read.csv("./res/WGCNA/wgcna_marker.csv")
marker <- merge(wgcna.marker, mfuzz.marker, by="Accession")

write.csv(marker, file = "./res/marker.csv", row.names = F)


gene.dat.293t <- read.csv("./res/lo2/gene_data_remove_batch_effect.csv")
gene.dat.lo2 <- read.csv("./res/293t/gene_data_p.csv")[1: 19]
zscoreTransform <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


gene.data.z.293t <- t(apply(gene.dat.293t[-1], 1, zscoreTransform))
gene.dat.293t <- cbind(gene.dat.293t[1], gene.data.z.293t)
gene.data.z.lo2 <- t(apply(gene.dat.lo2[-1], 1, zscoreTransform))
gene.dat.lo2 <- cbind(gene.dat.lo2[1], gene.data.z.lo2)


gene.dat.293t <- gene.dat.293t %>% 
  filter(Accession %in% marker$Accession)
gene.dat.lo2 <- gene.dat.lo2 %>% 
  filter(Accession %in% marker$Accession)

gene.noms <- read.csv("./res/lo2/Norm/data_noms.csv")
gene.dat.293t <- merge(gene.noms[1: 3], gene.dat.293t, by = "Accession")
row.names(gene.dat.293t) <- gene.dat.293t$GeneName
gene.dat.293t <- gene.dat.293t[-c(1: 3)]

gene.dat.lo2 <- merge(gene.noms[1: 3], gene.dat.lo2, by = "Accession")
row.names(gene.dat.lo2) <- gene.dat.lo2$GeneName
gene.dat.lo2 <- gene.dat.lo2[-c(1: 3)]

library(stringr)

t <- c(0, 3, 6, 9, 12, 24)
for (x in t) {
  gene.dat.293t <- gene.dat.293t %>%
    mutate(!!paste0("S", x) := rowMeans(across(starts_with(paste0("S", x, ".")))))
  gene.dat.lo2 <- gene.dat.lo2 %>%
    mutate(!!paste0("S", x) := rowMeans(across(starts_with(paste0("S", x, ".")))))
}
gene.dat.293t <- gene.dat.293t[19: 24]
gene.dat.293t <- gene.dat.293t %>% select(rev(colnames(gene.dat.293t)))
gene.dat.lo2 <- gene.dat.lo2[19: 24] 
gene.dat.lo2 <- gene.dat.lo2 %>% select(rev(colnames(gene.dat.lo2)))
library(circlize)
library(ComplexHeatmap)
circos.clear()  # 确保清空之前的图形

circos.par(gap.after = c(60))
#设置颜色 
color = colorRamp2(c(-3, 0, 3), c("#adc6ad","white","#d69c9b"))
#画图
circos.heatmap(as.matrix(gene.dat.lo2),
               col=color,
               dend.side="inside", #聚类放在环形内侧，inside为显示在圆环内圈，outside为显示在圆环外圈               
               dend.track.height = 0.18,
               track.height=0.2,
               rownames.side="outside", #基因名放在环形外侧,控制矩阵行名的方向,注意与dend.side不能在同一侧，必须一内一外 
               rownames.col="black",#行颜色 
               rownames.cex=0.8, #行字体大小  
               rownames.font=0.8, #行字体粗细 
               cell.border ="white", #边缘颜色 
               show.sector.labels = T,  #显示分的区域的标签 
               cluster=F #cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
)

circos.heatmap(as.matrix(gene.dat.293t),
               col=color,
               dend.side="inside", #聚类放在环形内侧，inside为显示在圆环内圈，outside为显示在圆环外圈               
               dend.track.height = 0.18,
               track.height=0.2,
               rownames.side="none", #基因名放在环形外侧,控制矩阵行名的方向,注意与dend.side不能在同一侧，必须一内一外 
               rownames.col="black",#行颜色 
               rownames.cex=0.8, #行字体大小  
               rownames.font=0.8, #行字体粗细 
               cell.border ="white", #边缘颜色 
               show.sector.labels = T,  #显示分的区域的标签 
               cluster=F #cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
)

# 添加文本到外圈
for(i in 1: ncol(gene.dat)){
  circos.text(x = -0.5,
              y = 0-0.2*i,
              labels = rep(colnames(gene.dat), 2)[i],
              facing = "clockwise", niceFacing = TRUE,
              cex = 0.5, adj = c(0, 0))
}


#画图例
lg = Legend(title="",col_fun=color, direction = c("horizontal"))
circle_size = unit(0.1, "snpc")
draw(lg, x = circle_size)
circos.clear()
write.csv(gene.dat.293t, "./res/293t/293t_marker_heatmap.csv")
write.csv(gene.dat.lo2, "./res/lo2/lo2_marker_heatmap.csv")

library(ggplot2)
library(reshape2)
# ------------- boxplot ---------------------
gene.dat.293t$Accession <- row.names(gene.dat.293t)
gene.data.293t.melt <- gene.dat.293t %>%
  melt(id.vars = "Accession", 
       variable.name = "Sample", 
       value.name = "zscore") %>% 
  mutate(Group = factor(paste(gsub("S(\\d+)", "\\1", Sample), "h"), levels = paste(c(0, 3, 6, 9, 12, 24), "h", sep = " ")))

fill <- c("0 h" = "#fff99e", "3 h" = "#bfff7f", "6 h" = "#ffbf7f", 
          "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
color <- c("0 h" = "#ffe00e", "3 h" = "#bfda00", "6 h" = "#ffbf0f", 
           "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
p <- ggplot(gene.data.293t.melt, aes(Group, zscore, fill = Group)) +
  geom_boxplot(linewidth = 0.7, alpha = 0.8, outlier.shape = NA) +
  geom_line(aes(group = Accession), 
            color = "grey70", size = 0.5, linetype = "dashed", position = position_dodge2(0.3), 
            alpha = 0.35) +
  geom_point(aes(group = Group, color = Group), position = position_dodge2(0.3), size = 1, alpha = 0.9) + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color)
p <- p + theme_bw() + 
  theme(panel.grid = element_blank())
p
ggsave("./res/293t/293t_marker_boxplot.pdf", p, width = 6, height = 3)


gene.dat.lo2$Accession <- row.names(gene.dat.lo2)
gene.data.lo2.melt <- gene.dat.lo2 %>%
  melt(id.vars = "Accession", 
       variable.name = "Sample", 
       value.name = "zscore") %>% 
  mutate(Group = factor(paste(gsub("S(\\d+)", "\\1", Sample), "h"), levels = paste(c(0, 3, 6, 9, 12, 24), "h", sep = " ")))

fill <- c("0 h" = "#fff99e", "3 h" = "#bfff7f", "6 h" = "#ffbf7f", 
          "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
color <- c("0 h" = "#ffe00e", "3 h" = "#bfda00", "6 h" = "#ffbf0f", 
           "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
p <- ggplot(gene.data.lo2.melt, aes(Group, zscore, fill = Group)) +
  geom_boxplot(linewidth = 0.7, alpha = 0.8, outlier.shape = NA) +
  geom_line(aes(group = Accession), 
            color = "grey70", size = 0.5, linetype = "dashed", position = position_dodge2(0.3), 
            alpha = 0.35) +
  geom_point(aes(group = Group, color = Group), position = position_dodge2(0.3), size = 1, alpha = 0.9) + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color)
p <- p + theme_bw() + 
  theme(panel.grid = element_blank())
p
ggsave("./res/lo2/lo2_marker_boxplot.pdf", p, width = 6, height = 3)

