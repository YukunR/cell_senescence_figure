library(dplyr)
library(tidyr)
library(edgeR)
library(reshape2)
library(ggplot2)
library(ggforce) # Accelerating 'ggplot2'
library(gghalves) # Compose Half-Half Plots Using Your Favourite Geoms
library(ggdist) # Visualizations of Distributions and Uncertainty
library(cowplot)


Normalization <- function(output.dir, 
                          gene.dat.na, 
                          sample.info, 
                          mycol, 
                          filter.threshold) {
  # >>> data preparation >>>
  
  gene.dat.na <- separate_rows(gene.dat.na, 
                               c(Accession, GeneName, Description), 
                               sep = ";")
  gene.dat.na <- as.data.frame(gene.dat.na)
  row.names(gene.dat.na) <- gene.dat.na$Accession
  
  gene.dat.noms <- subset(gene.dat.na, 
                          select = c(Accession, GeneName, Description, PG.Qvalue))
  gene.dat.na <- subset(gene.dat.na, 
                        select = -c(Accession, GeneName, Description, PG.Qvalue))
  
  write.csv(gene.dat.noms, 
            file = paste0(output.dir, "data_noms.csv"), 
            row.names = F)
  
  # >>> draw NA percentage >>>
  na.percentage <- colSums(is.na(gene.dat.na)) * 100 / nrow(gene.dat.na)
  names(na.percentage) <- colnames(gene.dat.na)
  pdf(paste0(output.dir, "na_percentage.pdf"), height = 5, width = 7)
  barplot(na.percentage, 
          ylim = c(0, min(max(na.percentage) + 20, 100)), 
          xlab = "Group", 
          ylab = "NA percentage (%)", 
          axisnames = F, 
          las = 2)  # rotate y axis label to horizontal
  text(x = barplot(na.percentage, plot = F), 
       y = na.percentage, 
       labels = round(na.percentage, 2), 
       pos = 3)
  text(barplot(na.percentage, plot = F), 
       par("usr")[3] - 1, 
       labels = names(na.percentage), 
       srt = 45, 
       adj = c(1.1, 1.1), 
       xpd = TRUE, 
       cex = 0.8)
  dev.off()
  
  # >>> normalization by col sums >>>
  norm.factors <- mean(colSums(gene.dat.na, na.rm = T)) / colSums(gene.dat.na, na.rm = T)
  gene.dat.na.norm <- sweep(gene.dat.na, 2, norm.factors, FUN = "*")
  
  # >>> log2 transition and imputation >>>
  gene.dat.na.norm.log <- log2(
    replace(gene.dat.na.norm, is.na(gene.dat.na.norm), 0)
  )
  
  # >>> draw distribution >>>
  pdf(paste0(output.dir, "before log.pdf"), height = 10, width = 10)
  par(mfrow = c(ceiling(ncol(gene.dat.na.norm) / 3), 3))
  for (colname in colnames(gene.dat.na.norm)) {
    hist(gene.dat.na.norm[[colname]], 
         main = colname, 
         xlab = "Intensity")
  }
  dev.off()
  
  pdf(paste0(output.dir, "after log.pdf"), height = 10, width = 10)
  par(mfrow = c(ceiling(ncol(gene.dat.na.norm.log) / 3), 3))
  for (colname in colnames(gene.dat.na.norm.log)) {
    hist(gene.dat.na.norm.log[[colname]], 
         main = colname, 
         xlab = "log2(Intensity)")
  }
  dev.off()
  
  # >>> draw density >>>
  gene.dat.na.log <- log2(
    replace(gene.dat.na, is.na(gene.dat.na), 0)
  )
  gene.dat.na.log.melt <- melt(gene.dat.na.log, 
                               variable.name = "Group", 
                               value.name = "Intensity")
  
  gene.dat.na.norm.log$Accession <- row.names(gene.dat.na.norm.log)
  gene.dat.na.norm.log.melt <- melt(gene.dat.na.norm.log, 
                                    id.vars = "Accession",
                                    variable.name = "Group", 
                                    value.name = "Intensity")
  
  
  # >>> imputation >>>
  gene.dat.na.norm.log.melt.draw <- gene.dat.na.norm.log.melt
  gene.dat.na.norm.log.melt[gene.dat.na.norm.log.melt == -Inf] <- NA
  sample.info <- sample.info %>% 
    add_count(Group, name = "total.count")
  gene.dat.impute <- merge(gene.dat.na.norm.log.melt, 
                           sample.info, 
                           by.x = "Group", 
                           by.y = "Sample", 
                           all.x = T)
  gene.dat.impute <- gene.dat.impute %>% 
    group_by(Accession, Group.y) %>% 
    mutate(non.na.count = sum(!is.na(Intensity))) %>% 
    mutate(non.na.percentage = non.na.count / total.count) %>% 
    ungroup() %>% 
    group_by(Accession) %>%
    mutate(non.na.percentage = min(non.na.percentage))
  gene.dat.impute <- gene.dat.impute[gene.dat.impute$non.na.percentage >= filter.threshold, ]
  
  gene.dat.impute <- dcast(gene.dat.impute, 
                           formula = Accession ~ Group, 
                           value.var = "Intensity")
  row.names(gene.dat.impute) <- gene.dat.impute$Accession
  gene.dat.impute <- subset(gene.dat.impute, select = -c(Accession))
  #' Impute missing values by random numbers drawn from a normal distribution
  #'
  #' Impute missing values by random numbers drawn from a normal distribution that has a down-shifted mean
  #' and shrunken standard deviation from the sample distribution. This is meant to be similar to imputation
  #' in the Perseus software.
  #'
  #' @param object Data frame or matrix containing filtered and log-transformed data.
  #' @param width Scale factor for the standard deviation of imputed distribution relative to the sample standard deviation.
  #' @param downshift Down-shifted the mean of imputed distribution from the sample mean, in units of sample standard deviation.
  #' @param seed Random seed
  #' @return A matrix with imputed values
  #' @export
  imputeNormal <- function(object, width=0.3, downshift=1.8, seed=100) {  # Perseus
    
    if (!is.matrix(object)) object <- as.matrix(object)
    mx <- max(object, na.rm = TRUE)
    mn <- min(object, na.rm = TRUE)
    if (mx - mn > 20) warning("Please make sure the values are log-transformed")
    
    set.seed(seed)
    object <- apply(object, 2, function(temp) {
      temp[!is.finite(temp)] <- NA
      temp.sd <- stats::sd(temp, na.rm=TRUE)
      temp.mean <- mean(temp, na.rm=TRUE)
      shrinked.sd <- width * temp.sd   # shrink sd width
      downshifted.mean <- temp.mean - downshift * temp.sd   # shift mean of imputed values
      n.missing <- sum(is.na(temp))
      temp[is.na(temp)] <- stats::rnorm(n.missing, mean=downshifted.mean, sd=shrinked.sd)
      temp
    })
    return(object)
  }
  
  sub.gene.dat.impute <- lapply(unique(sample.info$Group), function(group) {
    sample.names <- sample.info$Sample[sample.info$Group == group]
    col.names <- paste(sample.names, collapse = "|")
    subset <- gene.dat.impute[, grepl(col.names, colnames(gene.dat.impute))]
    return(subset)
  })
  names(sub.gene.dat.impute) <- unique(sample.info$Group)
  
  sub.gene.dat.imputed <- lapply(sub.gene.dat.impute, imputeNormal)
  gene.dat.imputed <- do.call(cbind, sub.gene.dat.imputed)
  
  gene.dat.save <- as.data.frame(gene.dat.imputed)
  gene.dat.save <- 2^(gene.dat.save)
  gene.dat.save$Accession <- row.names(gene.dat.save)
  gene.dat.save <- gene.dat.save[, c("Accession", setdiff(colnames(gene.dat.save), "Accession")), drop = FALSE]
  write.csv(gene.dat.save, 
            paste0(output.dir, "gene_data_imputation.csv"), 
            row.names = F)
  write.table(gene.dat.save, 
              file = paste0(output.dir, "gene_data_imputation.txt"), 
              row.names = F, 
              sep = "\t", 
              quote = F)
  
  # >>> plot imputated data >>>
  gene.dat.imputed.melt <- melt(as.data.frame(gene.dat.imputed), 
                                variable.name = "Group", 
                                value.name = "Intensity")
  
  drawDistribution <- function(melt.data, sign) {
    violin <- ggplot(melt.data, aes(x=Group, y=Intensity, fill=Group)) +
      geom_half_violin(position = position_nudge(x=0.25), side = "r", width=0.8, color=NA) +
      geom_jitter(aes(fill=Group, colour=Group), shape=21, size=0.5, width=0.2) +
      geom_boxplot(width=0.4, size=1.2, outlier.color=NA, linewidth = 0.4) +
      scale_fill_manual(values = mycol) +
      scale_color_manual(values = mycol) + 
      ylab("log2(Intensity)") + 
      xlab("Group") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = 'black',size=10),
            axis.title = element_text(color = 'black',size=12),
            legend.position = "none") + 
      coord_flip() + 
      scale_y_continuous(expand = c(0, 0))
    ggsave(violin, 
           filename = paste0(output.dir, sign, "_violin.pdf"), 
           width = 7, 
           height = 5)
    
    density <- ggplot(melt.data, aes(x=Intensity, color=Group)) +
      geom_density(show.legend = T, key_glyph = "timeseries", linewidth = 0.7) + 
      scale_color_manual(values = mycol) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = 'black',size=10),
            axis.title = element_text(color = 'black',size=12)) + 
      xlab("log2(Intensity)")
    ggsave(density, 
           filename = paste0(output.dir, sign, "_density.pdf"), 
           width = 7, 
           height = 5)
    figs[[length(figs) + 1]] <<- violin
    figs[[length(figs) + 1]] <<- density
  }
  figs <- list()
  drawDistribution(gene.dat.na.log.melt, sign = "raw")
  drawDistribution(gene.dat.na.norm.log.melt.draw, sign = "norm")
  drawDistribution(gene.dat.imputed.melt, sign = "impute")
  plot_grid(figs[[1]], figs[[2]], figs[[3]], figs[[4]], figs[[5]], figs[[6]], 
            ncol = 2, labels = LETTERS[1:6])
  ggsave(paste0(output.dir, "norm_imputed.pdf"), width = 10, height = 7)
  ggsave(paste0(output.dir, "norm_imputed.tiff"), width = 10, height = 7, dpi = 150)
}


#' @param data intensity matrix, first col = Accession, second col - last col = data
rowNormalization <- function(data) {
  noms.dat <- subset(data, select = Accession)
  gene.dat <- subset(data, select = -Accession)
  sum.target <- 100 * ncol(gene.dat)
  sum.dat <- rowSums(gene.dat)
  calibration.factor <- sum.dat / sum.target
  gene.dat <- sweep(gene.dat, 1, calibration.factor, FUN = "/")
  gene.dat <- cbind(noms.dat, gene.dat)
  return(gene.dat)
}