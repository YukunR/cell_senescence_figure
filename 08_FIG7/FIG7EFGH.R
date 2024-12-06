source("./GO.R", local = T)
source("./KEGG.R", local = T)
output.dir <- "./res/new/GOKEGG/"
go.background <- read.csv("./data/all_uniprot_go_background.csv", header = T)
kegg.background <- read.delim("./data/pathfromKegg_hsa.txt")
kegg.background.sig <- read.delim("./data/pathfromKegg_signal_hsa.txt")



gene.dat.sign <- read.csv("./res/new/gene_data_sign.csv")


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
  tryCatch({
    go.barplot <- goBarPlot(go.res)
    go.barplot
    ggsave(paste0(output.dir, "GObarplot_", sig, ".pdf"), 
           go.barplot, 
           width = 8, 
           height = 6)
  }, error = function(e) {e})
  
  
  kegg.res <- runKEGGAnalysis(gene.dat.tmp, kegg.background)
  write.csv(kegg.res, paste0(output.dir, "kegg_result_", sig, ".csv"), row.names = F)
  tryCatch({
    kegg.plot <- drawKEGGEnrichment(kegg.res)
    ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".pdf"), 
           kegg.plot$dot, 
           width = 6, 
           height = 6)
    ggsave(paste0(output.dir, "KEGGbarplot_", sig, ".pdf"), 
           kegg.plot$bar, 
           width = 6, 
           height = 6)
  }, error = function(e) {e})
  
  
  kegg.res <- runKEGGAnalysis(gene.dat.tmp, kegg.background.sig)
  write.csv(kegg.res, paste0(output.dir, "kegg_result_", sig, ".sig.csv"), row.names = F)
  tryCatch({
    kegg.plot <- drawKEGGEnrichment(kegg.res)
    ggsave(paste0(output.dir, "KEGGdotplot_", sig, ".sig.pdf"), 
           kegg.plot$dot, 
           width = 6, 
           height = 6)
    ggsave(paste0(output.dir, "KEGGbarplot_", sig, ".sig.pdf"), 
           kegg.plot$bar, 
           width = 6, 
           height = 6)
  }, error = function(e) {e})
  
}