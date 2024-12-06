library(clusterProfiler)
library(ggplot2)
library(enrichplot)


#' @param gene.list include "Gene" column
#' @param kegg.background with column -- gene, ontology, ontology name
#' @param go.p.threshold default 1
runKEGGAnalysis <- function(gene.list, kegg.background, go.p.threshold = 1) {
  term2gene <- kegg.background[, c('PATH', 'UNIPROT')]
  term2name <- kegg.background[, c('PATH', 'PATHNAME')]
  KEGGenrich <- enricher(
    gene = gene.list$Gene, 
    TERM2GENE = term2gene, 
    TERM2NAME = term2name, 
    pAdjustMethod = 'BH', 
    pvalueCutoff = 1, 
    qvalueCutoff = 1
  )
  return(KEGGenrich)
}


#' @param enrichment kegg enrichment result
#' @param display.number the number of show.category. default 10
drawKEGGEnrichment <- function(enrichment, display.number = 10) {
  dot.plot <- enrichplot::dotplot(enrichment,
                                  x = "Count", 
                                  color = "pvalue", 
                                  showCategory = display.number)
  bar.plot <- barplot(enrichment,
                      x = "qvalue", 
                      color = "qvalue", 
                      showCategory = display.number)
  return(list(dot = dot.plot, bar = bar.plot))
}
