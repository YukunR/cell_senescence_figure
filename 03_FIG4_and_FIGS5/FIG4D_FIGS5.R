# ========================
# ==== FIG4D & FIGS5 =====
# ========================

library(dplyr)
library(reshape2)

go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
gene.data <- read.delim("./res/lo2/gene_data_p.txt")
gene.noms <- read.csv("./res/lo2/Norm/data_noms.csv")
sample.info <- read.delim("./data/sample_info.txt")

gene.data.noms <- merge(gene.noms[-4], gene.data[1: 19], by = "Accession")


getGenefromeGOBackgroundbyPath <- function(path, go.background) {
  res <- go.background[go.background$GO == path, ]
  return(list("Accession" = res$UNIPROT, "PathName" = unique(res$GONAME)))
}

zscoreTransform <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


gene.data.z <- t(apply(gene.data.noms[-c(1: 3)], 1, zscoreTransform))
gene.data.z <- cbind(gene.data.noms[1: 3], gene.data.z)

gene.data.list <- list()

selected.go <- paste0("GO:", 
                      c("0005743", "0005741", "0031966", "0005759", 
                        "0016604", "0016363", "0005730", "0005788",  
                        "0005793", "0033106", "0005797", "0005802", 
                        "0005764"))
for (i in selected.go) {
  select.gene.info.tmp <- getGenefromeGOBackgroundbyPath(i, go.background)
  select.gene.tmp <- gene.data.z[gene.data.z$Accession %in% select.gene.info.tmp$Accession, ]
  select.gene.tmp$GOTerm <- select.gene.info.tmp$PathName
  gene.data.list[[select.gene.info.tmp$PathName]] <- select.gene.tmp
}

# >>> 改变数据格式并绘图 >>>
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(ComplexHeatmap)
library(circlize)
for (i in names(gene.data.list)) {
  gene.data.tmp <- gene.data.list[[i]]
  gene.data.tmp.melt <- gene.data.tmp %>% 
    melt(id.vars = c("Accession", "GeneName", "Description"), 
         variable.name = "Sample", 
         value.name = "Intensity") %>% 
    merge(sample.info, by = "Sample") %>% 
    group_by(Accession, Group) %>% 
    mutate(mean.z = mean(Intensity)) %>% 
    ungroup() %>% 
    dplyr::select(-Intensity, -Sample) %>% 
    unique() %>% 
    mutate(Group = factor(Group, levels = paste(c(0, 3, 6, 9, 12, 24), "h", sep = " ")))
  # gene.data.tmp.melt.fine <- gene.data.tmp.melt %>% 
  #   arrange(mean.z) %>% 
  #   mutate(rank = rank(mean.z, ties.method = "first")) %>% 
  #   group_by(mean.z) %>%
  #   mutate(mean.z = if(n() > 1) mean.z + 1e-3 * (rank - 1) else mean.z) %>%
  #   ungroup() %>%
  #   select(-rank)  # 微调使y不同
  
  average.z.group <- gene.data.tmp.melt %>% 
    dplyr::select(Accession, Group, mean.z) %>% 
    group_by(Group) %>% 
    mutate(mean.z.group = mean(mean.z)) %>% 
    ungroup() %>% 
    unique()
  
  # >>>  绘图 >>>
  fill <- c("0 h" = "#fff99e", "3 h" = "#bfff7f", "6 h" = "#ffbf7f", 
            "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
  color <- c("0 h" = "#ffe00e", "3 h" = "#bfda00", "6 h" = "#ffbf0f", 
             "9 h" = "#ff7f7f", "12 h" = "#717fff", "24 h" = "#9400d3")
  p <- ggplot(gene.data.tmp.melt, aes(Group, mean.z, fill = Group)) +
    geom_boxplot(linewidth = 0.7, alpha = 0.8, outlier.shape = NA) +
    geom_line(aes(group = Accession), 
              color = "grey70", size = 0.5, linetype = "dashed", position = position_dodge2(0.3), 
              alpha = 0.35) +
    geom_point(aes(group = Group, color = Group), position = position_dodge2(0.3), size = 1, alpha = 0.9) + 
    geom_line(data = average.z.group,
              mapping = aes(x = Group, y = mean.z.group, group = 1),
              color = "black", linetype = "twodash", linewidth = 1) + 
    
    scale_fill_manual(values = fill) + 
    scale_color_manual(values = color) + 
    theme_bw() + 
    ylim(c(-2.2, max(gene.data.tmp.melt$mean.z) + 2.3)) + 
    labs(x = "Time", y = "Z-score") + 
    ggtitle(i) +
    theme(plot.title = element_text(size = 14, face = "bold", color = "black"),
          panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))
  
  p + stat_compare_means(comparisons = list(c("0 h", "3 h"), c("3 h", "6 h"), c("6 h", "9 h"), 
                                            c("9 h", "12 h"), c("12 h", "24 h"), c("0 h", "6 h"), 
                                            c("6 h", "24 h"), c("0 h", "24 h")),
                         method = "t.test", paired = T, 
                         step.increase = c(0, 0.05, 0, 0.015, 0, 0.05, 0.05, 0.07))
  ggsave(paste0("./res/lo2/BoxPlot/", i, ".pdf"), width = 6, height = 2.75)
  
  
  average.z.group.long <- dcast(average.z.group, Accession ~ Group, value.var = "mean.z")
  pdf(paste0("./res/lo2/BoxPlot/", i, "_heatmap.pdf"), width = 6, height = 2.25)
  heatmap <- densityHeatmap(average.z.group.long[-1], ylab = "Z-score", 
                            column_names_rot = 0, title = i)
  draw(heatmap)
  dev.off()
  
  
}
