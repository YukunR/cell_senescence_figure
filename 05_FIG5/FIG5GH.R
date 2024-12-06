library(NPARC)
library(reshape2)
library(dplyr)
library(tidyr)
# load("./degradation.RData")
load("./synthesis.RData")
modelMetrics <- fits.synthesis$metrics  # fits.degredation$metrics  fits.synthesis$metrics 
fStats <- NPARCtest(modelMetrics, dfType = "empirical")
model.dat <- modelMetrics[!is.na(modelMetrics$group), ]
model.dat$tm <- log(2) / model.dat$k

# >>> 根据误差传播原理计算FC和FC的标准偏差 >>>
model.dat.1 <- model.dat %>%
  select(id, tm, group) %>% 
  spread(key = group, value = tm) %>% 
  na.omit() %>% 
  filter(S >= 0 & S <= 2000 & Y >= 0 & Y <= 2000)
model.dat.2 <- model.dat %>%
  select(id, t.half_sd, group) %>% 
  spread(key = group, value = t.half_sd) %>% 
  na.omit() %>% 
  merge(model.dat.1[1])

tm.data <- merge(model.dat.1, model.dat.2, 
                 by = "id", suffixes = c("", "_sd")) %>% 
  rename(c("Accession" = "id")) %>% 
  mutate(FC = S / Y, 
         FC_sd = abs(FC) * sqrt((S_sd / S) ^ 2 + (Y_sd / Y) ^ 2), 
         log2FC = log2(S / Y), 
         log2FC_sd = sqrt((S_sd / (S * log(2))) ^ 2) + (Y_sd / (Y * log(2))) ^ 2)


source("./coveragePlot.R", local = T)
output.dir <- "./res/"
gene.bin <- calculateCoverage(tm.data)
threshold <- thresholdInference(gene.bin)
coverage <- coveragePlot(gene.bin)
coverage
ggsave(paste0(output.dir, "Coverage.pdf"), 
       coverage, 
       width = 5, 
       height = 3)

source("./volcano.R", local = T)
gene.noms <- read.csv("./res/NPARC/gene_noms.csv")
gene.noms$GeneName <- gsub(".*GN=([^ ]+).*", "\\1", gene.noms$Description)
gene.noms$ProteinName <- gsub("(.*) OS=.*", "\\1", gene.noms$Description)
gene.dat.sign <- regulationGrouping(tm.data, 
                                    foldchange.threshold = threshold, 
                                    sd = "log2FC_sd", 
                                    sd.threshold = median(tm.data$log2FC_sd))  # median(tm.data$log2FC_sd)
volcano <- volcanoPlot(gene.dat.sign, 
                       gene.noms = gene.noms, 
                       foldchange.threshold = threshold, 
                       sd = "log2FC_sd", 
                       sd.threshold = median(tm.data$log2FC_sd), 
                       show.annotation = c(20, 20), sign.add = "P06748")
volcano
volcano <- volcano + 
  coord_cartesian(ylim = c(0, 6.5), xlim = c(-7, 7))
volcano
gene.dat.sign <- merge(gene.noms, gene.dat.sign, 
                       by = "Accession")
write.csv(gene.dat.sign, 
          paste0(output.dir, "gene_data_sign_syn.csv"), 
          row.names = F)
ggsave(paste0(output.dir, "volcano_syn.pdf"), 
       volcano, 
       width = 9, 
       height = 9)