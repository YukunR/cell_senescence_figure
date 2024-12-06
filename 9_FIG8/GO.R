library(clusterProfiler)
library(scales)
library(ggplot2)


#' @param gene.list include "Gene" column
#' @param go.background with column -- gene, ontology, ontology name
#' @param go.p.threshold default 1
runGOAnalysis <- function(gene.list, go.background, go.p.threshold = 1) {
  term2gene.bp <- go.background[go.background[, "ONTOLOGY"] == "BP", c("GO", "UNIPROT")] # 将属于BP的内容提出，第一列为go.id 第二列为gene.id 
  term2name.bp <- go.background[go.background[, "ONTOLOGY"] == "BP", c("GO", "GONAME")] # 与上一样，第一列为go.id 第二列为go.term
  term2gene.cc <- go.background[go.background[, "ONTOLOGY"] == "CC", c("GO", "UNIPROT")]
  term2name.cc <- go.background[go.background[, "ONTOLOGY"] == "CC", c("GO", "GONAME")]
  term2gene.mf <- go.background[go.background[, "ONTOLOGY"] == "MF", c("GO", "UNIPROT")]
  term2name.mf <- go.background[go.background[, "ONTOLOGY"] == "MF", c("GO", "GONAME")]
  
  # >>> 进行富集 >>>
  BPenrich <- enricher(
    gene = gene.list$Gene, 
    TERM2GENE = term2gene.bp, 
    TERM2NAME = term2name.bp, 
    pAdjustMethod = "BH", 
    pvalueCutoff = go.p.threshold, 
    qvalueCutoff = 1
  )
  
  MFenrich <- enricher(
    gene = gene.list$Gene, 
    TERM2GENE = term2gene.mf, 
    TERM2NAME = term2name.mf, 
    pAdjustMethod = "BH", 
    pvalueCutoff = go.p.threshold, 
    qvalueCutoff = 1
  )
  
  CCenrich <- enricher(
    gene = gene.list$Gene, 
    TERM2GENE = term2gene.cc, 
    TERM2NAME = term2name.cc, 
    pAdjustMethod = "BH", 
    pvalueCutoff = go.p.threshold, 
    qvalueCutoff = 1
  )
  
  # >>> 改变参数 >>>
  BPenrich <- as.data.frame(BPenrich)
  names(BPenrich)[6] <- "padj"
  BPenrich$Category <- rep("BP", time = nrow(BPenrich))
  
  CCenrich <- as.data.frame(CCenrich)
  names(CCenrich)[6] <- "padj"
  CCenrich$Category <- rep("CC", time = nrow(CCenrich))
  
  MFenrich <- as.data.frame(MFenrich)
  names(MFenrich)[6] <- "padj"
  MFenrich$Category <- rep("MF", time = nrow(MFenrich))
  
  GOenrich <- rbind(BPenrich, CCenrich, MFenrich)
  return(GOenrich)
}

#' @param enrichment go enrichment result
#' @param display.number a numeric vector contains number of BP, CC, MF. default 10, 5, 5
goBarPlot <- function(enrichment, display.number = c(10, 5, 5)) {
  go.result.BP <- head(enrichment[enrichment$Category == "BP", ], display.number[1])
  go.result.CC <- head(enrichment[enrichment$Category == "CC", ], display.number[2])
  go.result.MF <- head(enrichment[enrichment$Category == "MF", ], display.number[3])
  go.result <- rbind(go.result.BP, go.result.CC, go.result.MF)
  go.result$GeneRatio <- as.numeric(gsub(
    "([0-9]*)/[0-9]*", 
    "\\1", 
    go.result$GeneRatio)) / 
    as.numeric(gsub(
      "[0-9]*/([0-9]*)", 
      "\\1", 
      go.result$GeneRatio))
  
  # >>> x scaling function >>>
  min.x <- min(log10(go.result$padj))
  max.x <- max(go.result$GeneRatio)
  
  scale.x.0 <- function(x) {
    if (max.x > abs(min.x)) {
      ifelse(x <= 0, x * (max.x / abs(min.x)), x)
    } else {
      ifelse(x > 0, x * (abs(min.x) / max.x), x)}
  }
  scale.x.0.inverse <- function(x) {
    if (max.x > abs(min.x)) {
      ifelse(x <= 0, x / (max.x / abs(min.x)), x)
    } else {
      ifelse(x > 0,  x / (abs(min.x) / max.x), x)}
  }
  trans.x <- trans_new(name = "scale.x.0", 
                       transform = scale.x.0, 
                       inverse = scale.x.0.inverse)
  # >>> label setting function >>>
  setLabel <- function(x) {
    if (x < 0.05) {return(seq(0, 0.05, 0.01))}
    else if (x < 0.06) {return(seq(0, 0.06, 0.02))} 
    else if (x < 0.075) {return(seq(0, 0.075, 0.025))} 
    else if (x < 0.08) {return(seq(0, 0.08, 0.02))} 
    else if (x < 0.1) {return(seq(0, 0.1, 0.02))} 
    else if (x < 0.12) {return(seq(0, 0.12, 0.04))} 
    else if (x < 0.15) {return(seq(0, 0.15, 0.03))} 
    else if (x < 0.18) {return(seq(0, 0.18, 0.03))} 
    else if (x < 0.2) {return(seq(0, 0.2, 0.04))} 
    else if (x < 0.25) {return(seq(0, 0.25, 0.05))} 
    else if (x < 0.3) {return(seq(0, 0.3, 0.06))} 
    else if (x < 0.35) {return(seq(0, 0.35, 0.07))} 
    else if (x < 0.4) {return(seq(0, 0.4, 0.08))} 
    else if (x < 0.45) {return(seq(0, 0.45, 0.09))} 
    else if (x < 0.5) {return(seq(0, 0.5, 0.1))} 
    else if (x < 0.6) {return(seq(0, 0.6, 0.12))} 
    else if (x < 0.7) {return(seq(0, 0.7, 0.14))} 
    else if (x < 0.8) {return(seq(0, 0.8, 0.16))} 
    else if (x < 0.9) {return(seq(0, 0.9, 0.18))} 
    else if (x <= 0.1) {return(seq(0, 0.1, 0.2))}
  }
  
  # >>> draw bar plot >>>
  go.result$Description <- factor(
    go.result$Description,
    levels = rev(go.result$Description),
    ordered = F
  )
  
  enrichPlot <- ggplot() + 
    geom_col(data = go.result, 
             mapping = aes(x = GeneRatio, y = Description, fill = "#fbe5b8")) + 
    geom_col(data = go.result, 
             mapping = aes(x = log10(padj), y = Description, fill = Category, alpha = -log10(padj))) + 
    scale_fill_manual(values = c("#fbe5b8", "#85a8d3", "#c2dcbf", "#F19e9c"), 
                      labels = c("GeneRatio", "BP", "CC", "MF")) + 
    scale_alpha_continuous(range = c(0.2, 1), guide = "none") + 
    scale_x_continuous(trans = trans.x, 
                       limits = c(floor(min.x), max(setLabel(max.x))), 
                       breaks = c(seq(floor(min.x), 0, abs(floor(min.x))/ 5), setLabel(max.x)), 
                       labels = c(-seq(floor(min.x), 0, abs(floor(min.x))/ 5), setLabel(max.x))) + 
    scale_y_discrete(expand = expansion(add = c(1.3, 0.7))) + 
    annotate("text", label = "-log10(padj)", 0.55*floor(min.x), 0.3, size = 3) +
    annotate("text", label = "Gene Ratio", 0.55*max(setLabel(max.x)), 0.3, size = 3) + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          axis.title = element_blank(), 
          text = element_text(size = 10))
  return(enrichPlot)
}

