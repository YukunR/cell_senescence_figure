library(dplyr)
library(tidyr)
load("NPARCFit.RData")


model.dat <- modelMetrics[!is.na(modelMetrics$group), ]
model.dat <- model.dat %>%
  select(id, tm, group) %>%
  spread(key = group, value = tm)
model.dat$FC <- model.dat$SEN / model.dat$YOU
model.dat <- merge(model.dat, select(fStats, c(id, pVal)), by = "id")
model.dat <- na.omit(model.dat)
colnames(model.dat) <- c("Accession", "SEN", "YOU", "FC", "p")

source("./coveragePlot.R")
gene.bin <- calculateCoverage(model.dat)
threshold <- thresholdInference(gene.bin)
coverage <- coveragePlot(gene.bin)
coverage
output.dir <- "./res/new/"
ggsave(paste0(output.dir, "Coverage.pdf"), 
       coverage, 
       width = 5, 
       height = 3)

source("volcano.R", local = T)
colnames(model.dat)[1] <- "Accession"
model.sign <- regulationGrouping(model.dat, 
                                 foldchange.threshold = threshold, 
                                 p = "p", 
                                 p.threshold = 0.05)
gene.noms <- read.csv("./res/NPARC/gene_noms.csv")
volcano <- volcanoPlot(model.sign, 
                       foldchange.threshold = threshold, 
                       p = "p", 
                       sort.by = "p", 
                       p.threshold = 0.05, 
                       gene.noms = gene.noms, 
                       show.annotation = c(20, 20))
volcano <- volcano + coord_cartesian(xlim = c(-0.5, 0.5))


gene.dat.sign <- merge(gene.noms[-3], model.sign, 
                       by = "Accession")
write.csv(gene.dat.sign, 
          paste0("./res/gene_data_sign.csv"), 
          row.names = F)
ggsave("./res/volcano_p.pdf", 
       volcano, 
       width = 10, 
       height = 10)