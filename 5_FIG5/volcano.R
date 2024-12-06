library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#' @description
#' group proteins by their fold change and p value
#' @param gene.dat gene data with FC and p value
#' @param foldchange.threshold 
#' @param p should be one of p, p.adjust or q, default p
#' @param p.threshold
#' FC >= foldchange.threshold and p < p.threshold-- up
#' FC <= 2^(-log2(foldchange.threshold)) and p < p.threshold -- down
#' Others -- no
#' @return gene data with group
regulationGrouping <- function(gene.dat, foldchange.threshold, sd = "log2FC_sd", sd.threshold = 0.05) {
  gene.dat$Group <- ifelse(gene.dat$log2FC >= log2(foldchange.threshold) & gene.dat[[sd]] <= sd.threshold, "up", 
                           ifelse(gene.dat$log2FC <= -log2(foldchange.threshold) & gene.dat[[sd]] <= sd.threshold, "down", 
                                  "no"))
  return(gene.dat)
}

#' @description
#' group proteins by their fold change and p value
#' @param gene.dat gene data with FC and p value
#' @param foldchange.threshold 
#' @param sd should be one of p, p.adjust or q, default p
#' @param sd.threshold
#' @param gene.noms to convert UniprotID to GeneName
#' @param show.annotation a 2 components numeric vector determine the number of 
#'                        up- and down-reuglated gene annotation showed in the volcano plot.
#'                        Only used when gene.noms are provided.
#' @param sort.by should be one of "FC" or "Distance". 
#'                Only used when gene.noms are provided.
#' @param sign.add add genename apart from up and down, in format of Accession
#' @return volcano plot, ggplot2 object
volcanoPlot <- function(gene.dat, 
                        foldchange.threshold, 
                        sd = "log2FC_sd", 
                        sd.threshold = 0.05, 
                        gene.noms = NULL, 
                        show.annotation = c(10, 10), 
                        sort.by = "FC", 
                        sign.add = NULL) {
  volcano.dat <- data.frame(Gene = gene.dat$Accession, 
                            logFC = gene.dat$log2FC, 
                            pval = 1 / gene.dat[[sd]], 
                            Group = gene.dat$Group)
  up.count <- nrow(volcano.dat[which(volcano.dat$Group == "up"), ])
  down.count <- nrow(volcano.dat[which(volcano.dat$Group == "down"), ])
  no.count <- nrow(volcano.dat[which(volcano.dat$Group == "no"), ])
  
  if (!is.null(gene.noms)) {
    if (sort.by == "FC") {
      sort.by <- "logFC"
    } else if (sort.by == "Distance") {
      volcano.dat$Distance <- volcano.dat$logFC^2 + volcano.dat$pval^2
    } else {
      sort.by <- "logFC"
      warning("invalid sort.by argument, use FC instead.")
    }
    
    sign.up <- volcano.dat[which(volcano.dat$Group == "up"), ]
    sign.up <- head(sign.up[order(sign.up[[sort.by]], decreasing = T),], show.annotation[1])
    sign.down <- volcano.dat[which(volcano.dat$Group == "down"), ]
    if (sort.by == "logFC") {
      sign.down <- head(sign.down[order(sign.down[[sort.by]]),], show.annotation[2]) 
    } else {
      sign.down <- head(sign.down[order(sign.down[[sort.by]], decreasing = T),], show.annotation[2]) 
    }
    sign <- rbind(sign.up, sign.down)
    if (!is.null(sign.add)) {
      sign.add <- volcano.dat[volcano.dat$Gene %in% sign.add, ]
      sign <- rbind(sign, sign.add)
    }
    sign.all <- merge(sign, gene.noms, 
                      by.x = "Gene", by.y = "Accession")
  }
  
  theme_set(theme_bw())
  p <- ggplot(volcano.dat, aes(logFC, pval, color = Group)) + 
    geom_point(size = 1, aes(color = Group)) + 
    xlim(-max(abs(volcano.dat$logFC)), max(abs(volcano.dat$logFC))) + 
    labs(x = "log2FC", y = paste0("1 / ", sd), 
         title = paste0("UP ", up.count, "  DOWN ", down.count))
  p <- p + scale_color_manual(values = c("down" = "#0072B5", "no" = 'grey', "up" = '#BC3C28'))+
    geom_hline(yintercept = c(1 / sd.threshold), linetype = 4) +
    geom_vline(xintercept = c(-log2(foldchange.threshold), log2(foldchange.threshold)), linetype = 4)
  p <- p + theme(panel.grid = element_blank()) +
    theme(axis.line = element_line(size = 0))
  if (!is.null(gene.noms)) {
    p <- p + geom_text_repel(  #  标注基因
      data = sign.all, 
      aes(label = GeneName),
      size = 3,
      color = "black",
      segment.color = "black", show.legend = FALSE, 
      min.segment.length = 0)
    p <- p + guides()
  }
  volcano <- p + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 10)) + 
    scale_y_continuous(expand = c(0, 0))
  
  return(volcano)
}