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
# >>> get histone protein >>>
gene.data.histone <- gene.data.z[grepl("histone", tolower(gene.data.z$Description)), ]
gene.data.histone.mod <- gene.data.histone[!(grepl("histone h", tolower(gene.data.histone$Description)) | gene.data.histone$Accession == "O75367"), ]
gene.data.histone <- gene.data.histone[grepl("histone h", tolower(gene.data.histone$Description)) | gene.data.histone$Accession == "O75367", ]
gene.data.histone$GOTerm <- "histone"
gene.data.histone.mod$GOTerm <- "histone modification"
gene.data.list[["histone"]] <- gene.data.histone
gene.data.list[["histone modification"]] <- gene.data.histone.mod
# >>> get other protein by go >>>
selected.go <- paste0("GO:", 
                      c("0003700", "0006397", "0006401", "0003735", 
                        "0042273", "0042274", "0006412", "0036211", 
                        "0030163"))
for (i in selected.go) {
  select.gene.info.tmp <- getGenefromeGOBackgroundbyPath(i, go.background)
  select.gene.tmp <- gene.data.z[gene.data.z$Accession %in% select.gene.info.tmp$Accession, ]
  select.gene.tmp <- select.gene.tmp %>% mutate(GOTerm = select.gene.info.tmp$PathName)
  gene.data.list[[select.gene.info.tmp$PathName]] <- select.gene.tmp
}

# >>> 改变数据格式并绘图 >>>
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
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
    dplyr::select(Group, mean.z) %>% 
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
                         method = "t.test", 
                         paired = T,
                         step.increase = c(0, 0.05, 0, 0.015, 0, 0.05, 0.05, 0.07))
  ggsave(paste0("./res/lo2/BoxPlot/FIG2C/", i, ".pdf"), width = 6, height = 2.75)
}



# >>> 获取GO图 >>> 
load("./go.res.all.RData")
go.res.dataframe <- data.frame()
for (i in group) {
  go.res.tmp <- as.data.frame(go.res[[i]])
  go.res.tmp$Group <- i
  go.res.dataframe <- rbind(go.res.dataframe, go.res.tmp)
}

go.res.dataframe$GeneRatio <- as.numeric(gsub(
  "([0-9]*)/[0-9]*", 
  "\\1", 
  go.res.dataframe$GeneRatio)) / 
  as.numeric(gsub(
    "[0-9]*/([0-9]*)", 
    "\\1", 
    go.res.dataframe$GeneRatio))

go.res.dataframe.plot <- go.res.dataframe %>% 
  dplyr::filter(ID %in% selected.go)

go.res.dataframe.plot <- merge(go.res.dataframe, data.frame(ID = unique(go.res.dataframe.plot$ID)), by = "ID")
go.res.dataframe.plot$Group <- factor(go.res.dataframe.plot$Group, levels = unique(group))
go.res.dataframe.plot <- go.res.dataframe.plot[order(match(go.res.dataframe.plot$ID, selected.go)), ]
go.res.dataframe.plot$Description <- factor(go.res.dataframe.plot$Description, levels = rev(unique(go.res.dataframe.plot$Description)))
group.go.plot <- ggplot(go.res.dataframe.plot, aes(Group, Description)) + 
  geom_point(aes(color = pvalue, size = GeneRatio)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend(order = 1)) + 
  scale_color_gradient(low = "#BC3C28",high = "#0072B5") + 
  scale_size_continuous(range = c(3, 13), 
                        breaks = c(0.02, 0.04, 0.06), 
                        labels = c("0.02", "0.04", "0.06"), 
                        limits = c(0, 0.08))

ggsave("./res/lo2/BoxPlot/go.pdf", group.go.plot, width = 6, height = 10)