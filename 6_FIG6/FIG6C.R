library(dplyr)
library(reshape2)

go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
gene.dat.sign.syn <- read.csv("./res/FIG3/gene_data_sign_syn.csv")
gene.dat.sign.deg <- read.csv("./res/FIG3/gene_data_sign_deg.csv")

gene.data <- rbind(data.frame(gene.dat.sign.syn, group = "Syn"), 
                   data.frame(gene.dat.sign.deg, group = "Deg"))

getGenefromeGOBackgroundbyPath <- function(path, go.background) {
  res <- go.background[go.background$GO == path, ]
  return(list("Accession" = res$UNIPROT, "PathName" = unique(res$GONAME)))
}

gene.data.list <- list()
# >>> get histone protein >>>
gene.data.histone <- gene.data[grepl("histone", tolower(gene.data$Description)) & !grepl("non-histone", tolower(gene.data$Description)), ]
gene.data.histone.mod <- gene.data.histone[!(grepl("histone h", tolower(gene.data.histone$Description)) | gene.data.histone$Accession == "O75367"), ]
gene.data.histone <- gene.data.histone[grepl("histone h", tolower(gene.data.histone$Description)) | gene.data.histone$Accession == "O75367", ]
gene.data.list[["histone"]] <- gene.data.histone
gene.data.list[["histone modification"]] <- gene.data.histone.mod
# >>> get other protein by go >>>
selected.go <- paste0("GO:", 
                      c("0003700", "0006397", "0006401", "0003735", 
                        "0042273", "0042274", "0006412", "0036211", 
                        "0030163"))
for (i in selected.go) {
  select.gene.info.tmp <- getGenefromeGOBackgroundbyPath(i, go.background)
  select.gene.tmp <- gene.data[gene.data$Accession %in% select.gene.info.tmp$Accession, ]
  gene.data.list[[select.gene.info.tmp$PathName]] <- select.gene.tmp
}

# >>> 改变数据格式并绘图 >>>
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
for (i in names(gene.data.list)) {
  gene.data.tmp <- gene.data.list[[i]]
  gene.data.filtered <- gene.data.tmp
  outlier.up <- -9e10
  outlier.down <- 9e10
  for (var in c("S", "Y")) {
    Q1 <- quantile(gene.data.tmp[[var]], 0.25)
    Q3 <- quantile(gene.data.tmp[[var]], 0.75)
    outlier.up <- max(outlier.up, Q3 + (Q3 - Q1) * 1.5)
    outlier.down <- min(outlier.down, Q1 - (Q3 - Q1) * 1.5)
    cat(outlier.up, outlier.down)
  }
  
  gene.data.tmp.melt.filter <- gene.data.filtered %>% 
    select(c(Accession, Description, GeneName, S, Y, group)) %>% 
    melt(id.vars = c("Accession", "GeneName", "Description", "group"), 
         variable.name = "Group", 
         value.name = "tm") %>% 
    mutate(Group = factor(Group, levels = c("Y", "S")), 
           outlier = ifelse(tm > outlier.up | tm < outlier.down, T, F), 
           tm = ifelse(tm > outlier.up, outlier.up * 1.2, ifelse(tm < outlier.down, outlier.down / 1.2, tm)))
  gene.data.tmp.melt <- gene.data.tmp %>% 
    select(c(Accession, Description, GeneName, S, Y, group)) %>% 
    melt(id.vars = c("Accession", "GeneName", "Description", "group"), 
         variable.name = "Group", 
         value.name = "tm") %>% 
    mutate(Group = factor(Group, levels = c("Y", "S")))
  
  average.tm.group <- gene.data.tmp.melt %>% 
    dplyr::select(Group, tm) %>% 
    group_by(Group) %>% 
    mutate(mean.tm.group = mean(tm)) %>% 
    ungroup() %>% 
    select(-tm) %>% 
    unique()
  
  stat.test <- gene.data.tmp.melt %>% 
    group_by(group) %>% 
    t_test(tm ~ Group) %>% 
    add_xy_position(x = "group", dodge = 0.8)
  
  # >>>  绘图 >>>
  fill <- c("Y" = "#fff99e", "S" = "#9400d3")
  color <- c("Y" = "#ffe00e", "S" = "#9400d3")
  p <- ggplot(gene.data.tmp.melt, aes(group, tm)) +
    geom_boxplot(aes(fill = Group), linewidth = 0.7, alpha = 0.8, outlier.shape = NA, outliers = F) +
    geom_point(data = gene.data.tmp.melt.filter, 
               position = position_dodge2(0.6),# x轴方向上的闪避量
               aes(color = Group, shape = outlier), 
               size = 1.5) + 
    scale_fill_manual(values = fill) + 
    scale_color_manual(values = color) + 
    theme_bw() + 
    scale_y_continuous(limits = c(min(gene.data.tmp.melt.filter$tm) - 0.5, max(gene.data.tmp.melt.filter$tm) + 0.5), expand = c(0, 0)) + 
    labs(x = "Group", y = "T half (h)") + 
    ggtitle(i) +
    theme(plot.title = element_text(size = 14, face = "bold", color = "black"),
          panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), 
          legend.position = "bottom")
  
  p + stat_pvalue_manual(stat.test, y.position = outlier.up * 1.1, tip.length = 0) + coord_flip()
  ggsave(paste0("./res/FIG4/", i, ".pdf"), width = 4.5, height = 2.75)
}


source("./GO.R", local = T)
go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
gene.dat.sign.syn <- read.csv("./res/FIG3/gene_data_sign_syn.csv")
gene.dat.sign.deg <- read.csv("./res/FIG3/gene_data_sign_deg.csv")
gene.data.list <- list("syn" = gene.dat.sign.syn, "deg" = gene.dat.sign.deg)
go.res <- list()
for (group in c("syn", "deg")) {
  gene.data.tmp <- gene.data.list[[group]]
  for (sig in c("up", "down")) {
    gene.sig.tmp <- gene.data.tmp[gene.data.tmp$Group == sig, ] %>% 
      dplyr::rename(Gene = Accession)
    go.res[[paste(group, sig, sep = "_")]] <- runGOAnalysis(gene.sig.tmp, 
                                                            go.background)
  }
}

go.res.dataframe <- data.frame()
group <- names(go.res)
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
        axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend(order = 1)) + 
  scale_color_gradient(low = "#BC3C28",high = "#0072B5") + 
  scale_size_continuous(range = c(1, 10), 
                        breaks = c(0.04, 0.08, 0.12, 0.16), 
                        labels = c("0.04", "0.08", "0.12", "0.16"), 
                        limits = c(0, 0.15))

ggsave("./res/FIG4/go.pdf", group.go.plot, width = 6, height = 10)
