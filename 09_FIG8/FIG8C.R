library(dplyr)

gene.proteome <- read.csv("./data/gene_data_sign_24_0_h.csv")
gene.silac <- read.csv("./data/gene_data_sign_SILAC.csv")
gene.tpp <- read.csv("./data/gene_data_sign_TPP.csv")


gene.inter <- intersect(gene.proteome$Accession, gene.silac$Accession)
gene.inter <- intersect(gene.inter, gene.tpp$Accession)


gene.proteome <- filter(gene.proteome, Accession %in% gene.inter)
gene.silac <- filter(gene.silac, Accession %in% gene.inter)
gene.tpp <- filter(gene.tpp, Accession %in% gene.inter)

gene.list <- list("gene.proteome" = gene.proteome, "gene.silac" = gene.silac, "gene.tpp" = gene.tpp)

getGenebyGroup <- function(gene.data, group) {
  try({
    return(gene.data[gene.data$Group %in% group, ]$Accession)
  })
}

group <- c("up", "down")
venn.list <- list()
for (i in names(gene.list)) {
  venn.list[[i]] <- getGenebyGroup(gene.list[[i]], group)
}

library(ggplot2)
library(ggpolypath)


library(ggVennDiagram)
venn <- ggVennDiagram(venn.list, 
                      set_color = c("#fea05a", "#6f4f97", "#8fd01f"), 
                      category.names = c("Quantity", "T half", "Tm"), 
                      label = "count", 
)
venn

ggsave("./res/FIG6/venn.pdf", venn, width = 5, height = 4)

quantity.thalf <- intersect(venn.list$gene.proteome, venn.list$gene.silac)
quantity.tm <- intersect(venn.list$gene.proteome, venn.list$gene.tpp)
thalf.tm <- intersect(venn.list$gene.silac, venn.list$gene.tpp)