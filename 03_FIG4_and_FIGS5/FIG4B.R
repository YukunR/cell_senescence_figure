library(dplyr)
library(stringr)

gene.dat <- read.csv("./res/lo2/gene_data_p.csv")
gene.head.list <- list(
  c("Accession", "S0.1", "S0.2", "S0.3", "S3.1", "S3.2", "S3.3", "FC_3.h_0.h", "p_3.h.0.h"), 
  c("Accession", "S0.1", "S0.2", "S0.3", "S6.1", "S6.2", "S6.3", "FC_6.h_0.h", "p_6.h.0.h"), 
  c("Accession", "S0.1", "S0.2", "S0.3", "S9.1", "S9.2", "S9.3", "FC_9.h_0.h", "p_9.h.0.h"), 
  c("Accession", "S0.1", "S0.2", "S0.3", "S12.1", "S12.2", "S12.3", "FC_12.h_0.h", "p_12.h.0.h"), 
  c("Accession", "S0.1", "S0.2", "S0.3", "S24.1", "S24.2", "S24.3", "FC_24.h_0.h", "p_24.h.0.h")
)


source("./coveragePlot.R", local = T)
source("./volcano_2group.R", local = T)

gene.noms <- read.csv("./res/lo2/Norm/data_noms.csv")
output.dir <- "./res/lo2/GOofAllGene/"

gene.dat.up <- list()
gene.dat.down <- list()
gene.dat.all <- list()
for (i in 1: 5) {
  gene.head.tmp <- gene.head.list[[i]]
  gene.dat.select <- gene.dat %>% select(gene.head.tmp)
  colnames(gene.dat.select)[8: 9] <- c("FC", "p")
  gene.bin <- calculateCoverage(gene.dat.select)
  threshold <- thresholdInference(gene.bin)
  gene.dat.sign <- regulationGrouping(gene.dat.select, 
                                      foldchange.threshold = threshold, 
                                      p = "p", 
                                      p.threshold = 0.05)
  
  gene.dat.sign <- merge(gene.noms[-4], gene.dat.sign, 
                         by = "Accession")
  colnames(gene.dat.sign)[1] <- "Gene"
  write.csv(gene.dat.sign, 
            paste0(output.dir, "gene_data_sign_", str_remove(gene.head.tmp[8], "FC_"), ".csv"), 
            row.names = F)
  gene.dat.up[[i]] <- gene.dat.sign[gene.dat.sign$Group == "up", ]
  gene.dat.down[[i]] <- gene.dat.sign[gene.dat.sign$Group == "down", ]
  gene.dat.all[[i]] <- rbind(gene.dat.sign[gene.dat.sign$Group == "up", ], gene.dat.sign[gene.dat.sign$Group == "down", ])
}

source("./GO.R", local = T)

go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
group <- c("3 h", "6 h", "9 h", "12 h", "24 h")
go.res <- list()
for (i in 1: 5) {
  go.res[[group[i]]] <- runGOAnalysis(gene.dat.all[[i]], go.background)
}
sig <- "all"

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
  group_by(Group) %>% 
  slice_head(n = 3) %>% 
  ungroup()

go.res.dataframe.bp <- go.res.dataframe[go.res.dataframe$Category == "BP", ]

go.res.dataframe.bp.plot <- go.res.dataframe.bp %>% 
  group_by(Group) %>% 
  slice_head(n = 3) %>% 
  ungroup()

go.res.dataframe.bp.plot <- merge(go.res.dataframe.bp, data.frame(ID = unique(go.res.dataframe.bp.plot$ID)), by = "ID")
go.res.dataframe.bp.plot$Group <- factor(go.res.dataframe.bp.plot$Group, levels = unique(group))
go.res.dataframe.bp.plot <- go.res.dataframe.bp.plot[go.res.dataframe.bp.plot$pvalue < 0.05,]
group.go.plot.bp <- ggplot(go.res.dataframe.bp.plot, aes(Group, Description)) + 
  geom_point(aes(color = pvalue, size = GeneRatio)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend(order = 1)) + 
  scale_color_gradient(low = "#BC3C28",high = "#0072B5", limits = c(0, 0.03)) + 
  scale_size_continuous(range = c(1, 8), 
                        breaks = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.2), 
                        labels = c("0.01", "0.03", "0.05", "0.07", "0.09", "0.2"), 
                        limits = c(0, 0.25))


write.csv(go.res.dataframe, paste0(output.dir, "GO_result_", sig, ".csv"), row.names = F)
# ggsave(paste0(output.dir, "group_go_bp_", sig, ".pdf"), group.go.plot.bp, width = 7, height = 2)