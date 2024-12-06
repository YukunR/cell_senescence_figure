library(clusterProfiler)
library(dplyr)
library(enrichplot)

aggregate.background <- read.delim("./data/aggregate_background.txt")
aggregate.background <- aggregate.background %>% 
  mutate(Term = Aggregate)

term2gene <- aggregate.background[, c('Term', 'Accession')]
term2name <- aggregate.background[, c('Term', 'Aggregate')]

gene.data <- read.csv("./res/new/gene_data_sign.csv")
gene.data <- gene.data %>% 
  filter(!Group == "no")

Aggregate.enrich <- enricher(
  gene = gene.data$Accession, 
  TERM2GENE = term2gene, 
  TERM2NAME = term2name, 
  pAdjustMethod = 'BH', 
  pvalueCutoff = 1, 
  qvalueCutoff = 1
)

dot.plot <- enrichplot::dotplot(Aggregate.enrich,
                                x = "Count", 
                                color = "pvalue", 
                                showCategory = 15)

library(ggplot2)
ggsave("./res/new/aggregate.pdf", dot.plot, height = 6, width = 6)

a <- as.data.frame(Aggregate.enrich)
write.csv(a, "./res/new/aggregate.csv", row.names = F)
