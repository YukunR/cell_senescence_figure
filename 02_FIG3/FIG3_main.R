library(dplyr)

# ===================
# == preprocessing ==
# ===================
gene.data <- read.delim("./res/lo2/gene_data_p.txt")
gene.data.select <- gene.data %>% select(c(1: 4, 17: 19))
sample.info <- data.frame(Sample = c(colnames(gene.data.select)[-1]), 
                          Group = c(rep("0 h", 3), rep("24 h", 3)))

output.dir <- "./res/lo2/24_0h/"

source("./p.R")
gene.dat.select <- tTest(gene.data.select, sample.info, control.group.name = "0 h")
write.csv(gene.dat.select, 
          paste0(output.dir, "gene_data_p.csv"), 
          row.names = F)

# ===================
# ====== FIG3 A =====
# ===================
source("./coveragePlot.R", local = T)
gene.bin <- calculateCoverage(gene.dat.select)
threshold <- thresholdInference(gene.bin)
coverage <- coveragePlot(gene.bin)
coverage
ggsave(paste0(output.dir, "Coverage.pdf"), 
       coverage, 
       width = 5, 
       height = 3)


# ===================
# ====== FIG3 B =====
# ===================
# >>> Volcano plot >>>
source("./volcano_2group.R", local = T)
gene.noms <- read.csv("./res/lo2/Norm/data_noms.csv")
gene.dat.sign <- regulationGrouping(gene.dat.select, 
                                    foldchange.threshold = threshold, 
                                    p = "p", 
                                    p.threshold = 0.05)
volcano <- volcanoPlot(gene.dat.sign, 
                       gene.noms = gene.noms, 
                       foldchange.threshold = threshold, 
                       p = "p", 
                       p.threshold = 0.05, 
                       show.annotation = c(20, 20), 
                       sign.add = "P06748")
volcano

gene.dat.sign <- merge(gene.noms[-4], gene.dat.sign, 
                       by.x = "Accession", by.y = "Gene")
write.csv(gene.dat.sign, 
          paste0(output.dir, "gene_data_sign.csv"), 
          row.names = F)
ggsave(paste0(output.dir, "volcano.pdf"), 
       volcano, 
       width = 8, 
       height = 8)


# ===================
# === FIG3 C - G ====
# ===================
# >>> GO KEGG >>>
source("./GO.R", local = T)
source("./KEGG.R", local = T)
output.dir <- "./res/lo2/24_0h/GOKEGG/"
go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
kegg.background <- read.delim("./data/pathfromKegg_hsa.txt")
kegg.background.sig <- read.delim("./data/pathfromKegg_signal_hsa.txt")
colnames(gene.dat.sign)[1] <- "Gene"
gene.up <- gene.dat.sign[gene.dat.sign$Group == "up", ]
gene.down <- gene.dat.sign[gene.dat.sign$Group == "down", ]
gene.all <- rbind(gene.up, gene.down)
gene.list <- list(gene.up, gene.down, gene.all)
gene.sig <- c("up", "down", "all")
for (i in 1: 3) {
  sig <- gene.sig[i]
  gene.dat.tmp <- gene.list[[i]]
  go.res <- runGOAnalysis(gene.dat.tmp, 
                          go.background)
  write.csv(go.res, paste0(output.dir, "go_result_", sig, ".csv"), row.names = F)
  go.barplot <- goBarPlot(go.res)
  go.barplot
  ggsave(paste0(output.dir, "GObarplot_", sig, ".pdf"), 
         go.barplot, 
         width = 8, 
         height = 6)
  ggsave(paste0(output.dir, "GObarplot_", sig, ".tiff"), 
         go.barplot, 
         width = 8, 
         height = 6)
  
  
  
  kegg.res <- runKEGGAnalysis(gene.dat.tmp, kegg.background)
  write.csv(kegg.res, paste0(output.dir, "kegg_result_", sig, ".csv"), row.names = F)
  kegg.plot <- drawKEGGEnrichment(kegg.res)
  ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".pdf"), 
         kegg.plot$dot, 
         width = 6, 
         height = 6)
  ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".tiff"), 
         kegg.plot$dot, 
         width = 6, 
         height = 6)
  
  kegg.res <- runKEGGAnalysis(gene.dat.tmp, kegg.background.sig)
  write.csv(kegg.res, paste0(output.dir, "kegg_result_", sig, ".sig.csv"), row.names = F)
  kegg.plot <- drawKEGGEnrichment(kegg.res)
  ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".sig.pdf"), 
         kegg.plot$dot, 
         width = 6, 
         height = 6)
  ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".sig.tiff"), 
         kegg.plot$dot, 
         width = 6, 
         height = 6)
}