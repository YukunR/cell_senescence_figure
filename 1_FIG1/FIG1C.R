
# ===================
# ====== Mfuzz ======
# ===================

library(Mfuzz)
library(dplyr)
gene.dat.p <- read.csv("./res/293t/gene_data_p.csv")
gene.dat.time <- gene.dat.p %>% filter(ANOVA_p < 0.05) %>% select(1: 19)
group.names <- unique(sample.info$Group)
for (group.name in group.names) {
  columns <- names(gene.dat.time) %in% c(sample.info$Sample[sample.info$Group == group.name])
  gene.dat.time[[group.name]] <- rowMeans(gene.dat.time[columns])
}
gene.dat.time <- gene.dat.time %>% 
  select(all_of(c("Accession", group.names)))
row.names(gene.dat.time) <- gene.dat.time$Accession
gene.dat.time <- as.matrix(select(gene.dat.time, -Accession))

mfuzz.class <- new('ExpressionSet', exprs = gene.dat.time)
mfuzz.class <- standardise(mfuzz.class)

cluster.num <- 7
mfuzz.cluster <- mfuzz(mfuzz.class, c = cluster.num, m = mestimate(mfuzz.class))
library(RColorBrewer)
color <- colorRampPalette(rev(c("#df6561", "#f2d086", "#95be7e", "#57b1a9")))(1000)
pdf("./res/lo2/mFuzz/cluster.pdf", width = 4, height = 10)
mfuzz.plot2(mfuzz.class, cl = mfuzz.cluster, mfrow = c(7, 1), time.labels = colnames(gene.dat.time), 
            x11 = F, colo = "fancy", centre = T, Xwidth = 4, Xheight = 10)
dev.off()

png("./res/lo2/mFuzz/cluster.png", width = 4000, height = 2500, res = 300)
# mfuzz.plot(mfuzz.class, cl = mfuzz.cluster, mfrow = c(3, 3), time.labels = colnames(gene.dat.time), 
#            new.window = F, colo = color)
mfuzz.plot2(mfuzz.class, cl = mfuzz.cluster, mfrow = c(3, 3), time.labels = colnames(gene.dat.time), 
            x11 = F, colo = "fancy", centre = T)

dev.off()

gene.cluster <- data.frame(cluster = mfuzz.cluster$cluster)
gene.dat.time <- as.data.frame(gene.dat.time)
gene.dat.time$Accession <- row.names(gene.dat.time)
gene.cluster$Accession <- row.names(gene.cluster)
gene.dat.cluster <- merge(gene.dat.p, gene.cluster, by = "Accession")
write.csv(gene.dat.cluster, "./res/lo2/mFuzz/gene_dat_cluster.csv", row.names = F)

# ===================
# ===== heatmap =====
# ===================

library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(dendextend)
library(dplyr)
library(circlize)

z_score_transform <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

decorate_dendrogram <- function(dend, linewidth = 0.3) {
  dend %>% set("branches_lwd", linewidth) %>% set("branches_col", "#000000")  # 黑色线条
}


heatmap.data <- read.csv("./res/lo2/mFuzz/gene_dat_cluster.csv")
row.names(heatmap.data) <- heatmap.data$Accession
condition <- c(rep("0 h", 3), rep("3 h", 3), rep("6 h", 3), 
               rep("9 h", 3), rep("12 h", 3), rep("24 h", 3))
mycol <- c("0 h" = "#fff99e", "3 h" = "#bfff7f", "6 h" = "#ffbf7f", 
           "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
col_func = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

p <- NULL
for (c in 1: 7) {
  
  heatmap.data.tmp <- heatmap.data %>% filter(cluster == c)
  heatmap.data.tmp <- heatmap.data.tmp[2: 19]
  heatmap.data.tmp.z <- t(apply(heatmap.data.tmp, 1, z_score_transform))
  row_dendrogram <- as.dendrogram(hclust(dist(heatmap.data.tmp.z)))
  
  row_dendrogram <- decorate_dendrogram(row_dendrogram)
  
  
  # 使用ComplexHeatmap包绘制热图和聚类
  if (is.null(p)) {
    p <- Heatmap(heatmap.data.tmp.z,
                 cluster_rows = row_dendrogram,
                 col = col_func, 
                 cluster_columns = F,
                 show_row_names = F,
                 show_column_names = TRUE,
                 column_names_rot = 45, 
                 top_annotation = HeatmapAnnotation(df = data.frame(Condition = condition), col = list(Condition = mycol)),
                 right_annotation = rowAnnotation(foo = anno_barplot(rowMeans(heatmap.data.tmp.z), 
                                                                     numbers_rot = 45, 
                                                                     ylim = c(-2e-15, 2e-15), 
                                                                     axis = F, 
                                                                     labels_rot = 45), 
                                                  show_annotation_name = F),
                 heatmap_legend_param = list(title = "Z score", legend_position = "top"), 
                 column_names_gp = gpar(fontsize = 8), 
                 width = 5, height = 1.5, border = T)
    draw(p)
  } else {
    if (c == 7) {
      p <- p %v% Heatmap(heatmap.data.tmp.z,
                         cluster_rows = row_dendrogram,
                         col = col_func, 
                         cluster_columns = F,
                         show_row_names = F,
                         show_column_names = TRUE,
                         column_names_rot = 45,
                         right_annotation = rowAnnotation(foo = anno_barplot(rowMeans(heatmap.data.tmp.z), 
                                                                             numbers_rot = 45, 
                                                                             ylim = c(-2e-15, 2e-15), 
                                                                             axis = T, 
                                                                             labels_rot = 45), 
                                                          show_annotation_name = F), 
                         column_names_gp = gpar(fontsize = 8), 
                         width = 5, height = 1.5, show_heatmap_legend = F, border = T)
      draw(p)
    } else {
      p <- p %v% Heatmap(heatmap.data.tmp.z,
                         cluster_rows = row_dendrogram,
                         col = col_func, 
                         cluster_columns = F,
                         show_row_names = F,
                         show_column_names = TRUE,
                         column_names_rot = 45,
                         right_annotation = rowAnnotation(foo = anno_barplot(rowMeans(heatmap.data.tmp.z), 
                                                                             numbers_rot = 45, 
                                                                             ylim = c(-2e-15, 2e-15), 
                                                                             axis = F, 
                                                                             labels_rot = 45), 
                                                          show_annotation_name = F), 
                         column_names_gp = gpar(fontsize = 8), 
                         width = 5, height = 1.5, show_heatmap_legend = F, border = T)
    }
    
  }
}

pdf("./res/lo2/mFuzz/heatmap.pdf", width = 5, height = 7)
p
dev.off()


# ===================
# ======= GO ========
# ===================

gene.dat.cluster <- read.csv("./res/lo2/mFuzz/gene_dat_cluster.csv")
gene.noms <- read.csv("./res/lo2/Norm/data_noms.csv")
output.dir <- "./res/lo2/volcano/"

gene.dat.3h.0h <- gene.dat.cluster %>% select(c(Accession, p_3.h.0.h, FC_3.h_0.h, cluster))
colnames(gene.dat.3h.0h)[1: 3] <- c("Gene", "p", "FC")

gene.dat.6h.0h <- gene.dat.cluster %>% select(c(Accession, p_6.h.0.h, FC_6.h_0.h, cluster))
colnames(gene.dat.6h.0h)[1: 3] <- c("Gene", "p", "FC")

gene.dat.9h.0h <- gene.dat.cluster %>% select(c(Accession, p_9.h.0.h, FC_9.h_0.h, cluster))
colnames(gene.dat.9h.0h)[1: 3] <- c("Gene", "p", "FC")

gene.dat.12h.0h <- gene.dat.cluster %>% select(c(Accession, p_12.h.0.h, FC_12.h_0.h, cluster))
colnames(gene.dat.12h.0h)[1: 3] <- c("Gene", "p", "FC")

gene.dat.24h.0h <- gene.dat.cluster %>% select(c(Accession, p_24.h.0.h, FC_24.h_0.h, cluster))
colnames(gene.dat.24h.0h)[1: 3] <- c("Gene", "p", "FC")

source("./GO.R", local = T)
source("./KEGG.R", local = T)

go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
kegg.background <- read.delim("./data/pathfromKegg_hsa.txt")
kegg.background.sig <- read.delim("./data/pathfromKegg_signal_hsa.txt")

output.dir <- "./res/lo2/mfuzz/time_GO/"
group <- c("3 h", "6 h", "9 h", "12 h", "24 h")
go.dat <- list("3 h" = gene.dat.3h.0h, "6 h" = gene.dat.6h.0h, 
               "9 h" = gene.dat.9h.0h, "12 h" = gene.dat.12h.0h, 
               "24 h" = gene.dat.24h.0h)
for (c in 1: 7) {
  go.res <- list()
  for (i in group) {
    go.res[[i]] <- runGOAnalysis(data.frame(Gene = go.dat[[i]]$Gene[go.dat[[i]]$p < 0.05 & go.dat[[i]]$cluster == c]), go.background)
  }
  
  go.res.dataframe <- data.frame()
  for (i in group) {
    
    go.res[[i]]$Group <- i
    go.res.dataframe <- rbind(go.res.dataframe, go.res[[i]])
  }
  
  go.res.dataframe$GeneRatio <- as.numeric(gsub(
    "([0-9]*)/[0-9]*", 
    "\\1", 
    go.res.dataframe$GeneRatio)) / 
    as.numeric(gsub(
      "[0-9]*/([0-9]*)", 
      "\\1", 
      go.res.dataframe$GeneRatio))
  write.csv(go.res.dataframe, paste0(output.dir, "GO_result_cluster_", c, ".csv"), row.names = F)
  
  
  go.res.dataframe.bp <- go.res.dataframe[go.res.dataframe$Category == "BP", ]
  go.res.dataframe.mf <- go.res.dataframe[go.res.dataframe$Category == "MF", ]
  go.res.dataframe.cc <- go.res.dataframe[go.res.dataframe$Category == "CC", ]
  
  go.res.dataframe.bp.plot <- go.res.dataframe.bp %>% 
    group_by(Group) %>% 
    slice_head(n = 3) %>% 
    ungroup()
  
  go.res.dataframe.bp.plot <- merge(go.res.dataframe.bp, data.frame(ID = unique(go.res.dataframe.bp.plot$ID)), by = "ID")
  go.res.dataframe.bp.plot$Group <- factor(go.res.dataframe.bp.plot$Group, levels = c("0 h", "3 h", "6 h", "9 h", "12 h", "24 h"))
  go.res.dataframe.bp.plot <- go.res.dataframe.bp.plot[go.res.dataframe.bp.plot$pvalue < 0.05,]
  select.path <- unique(go.res.dataframe.bp.plot$ID)[1: 6]
  go.res.dataframe.bp.plot <- go.res.dataframe.bp.plot[go.res.dataframe.bp.plot$ID %in% select.path, ]
  
  group.go.plot.bp <- ggplot(go.res.dataframe.bp.plot, aes(Group, Description)) + 
    geom_point(aes(color = pvalue, size = GeneRatio)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "right") + 
    scale_color_gradient(low = "#BC3C28",high = "#0072B5", limits = c(0, 0.05)) + 
    labs(x = NULL, y = NULL) + 
    guides(size = guide_legend(order = 1)) + 
    scale_size_continuous(range = c(1, 8), 
                          breaks = c(0.03, 0.06, 0.09, 0.12, 0.15, 0.2), 
                          labels = c("0.03", "0.06", "0.09", "0.12", "0.15", "0.2"), 
                          limits = c(0, 0.2)) + 
    scale_y_discrete(position = "right")
  
  
  go.res.dataframe.mf.plot <- go.res.dataframe.mf %>% 
    group_by(Group) %>% 
    slice_head(n = 3) %>% 
    ungroup()
  
  go.res.dataframe.mf.plot <- merge(go.res.dataframe.mf, data.frame(ID = unique(go.res.dataframe.mf.plot$ID)), by = "ID")
  go.res.dataframe.mf.plot$Group <- factor(go.res.dataframe.mf.plot$Group, levels = c("0 h", "3 h", "6 h", "9 h", "12 h", "24 h"))
  go.res.dataframe.mf.plot <- go.res.dataframe.mf.plot[go.res.dataframe.mf.plot$pvalue < 0.05,]
  select.path <- unique(go.res.dataframe.mf.plot$ID)[1: 6]
  go.res.dataframe.mf.plot <- go.res.dataframe.mf.plot[go.res.dataframe.mf.plot$ID %in% select.path, ]
  group.go.plot.mf <- ggplot(go.res.dataframe.mf.plot, aes(Group, Description)) + 
    geom_point(aes(color = pvalue, size = GeneRatio)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "right") + 
    scale_color_gradient(low =  "#BC3C28",high ="#0072B5") + 
    labs(x = NULL, y = NULL) + 
    guides(size = guide_legend(order = 1))
  
  
  go.res.dataframe.cc.plot <- go.res.dataframe.cc %>% 
    group_by(Group) %>% 
    slice_head(n = 3) %>% 
    ungroup()
  
  go.res.dataframe.cc.plot <- merge(go.res.dataframe.cc, data.frame(ID = unique(go.res.dataframe.cc.plot$ID)), by = "ID")
  go.res.dataframe.cc.plot$Group <- factor(go.res.dataframe.cc.plot$Group, levels = c("0 h", "3 h", "6 h", "9 h", "12 h", "24 h"))
  go.res.dataframe.cc.plot <- go.res.dataframe.cc.plot[go.res.dataframe.cc.plot$pvalue < 0.05,]
  select.path <- unique(go.res.dataframe.cc.plot$ID)[1: 6]
  go.res.dataframe.cc.plot <- go.res.dataframe.cc.plot[go.res.dataframe.cc.plot$ID %in% select.path, ]
  group.go.plot.cc <- ggplot(go.res.dataframe.cc.plot, aes(Group, Description)) + 
    geom_point(aes(color = pvalue, size = GeneRatio)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "right") + 
    scale_color_gradient(low =  "#BC3C28",high ="#0072B5") + 
    labs(x = NULL, y = NULL) + 
    guides(size = guide_legend(order = 1))
  
  ggsave(paste0(output.dir, "group_go_bp_cluster", c, ".pdf"), group.go.plot.bp, width = 7, height = 1.5)
  ggsave(paste0(output.dir, "group_go_mf_cluster", c, ".pdf"), group.go.plot.mf, width = 7, height = 1.5)
  ggsave(paste0(output.dir, "group_go_cc_cluster", c, ".pdf"), group.go.plot.cc, width = 7, height = 1.5)
}

