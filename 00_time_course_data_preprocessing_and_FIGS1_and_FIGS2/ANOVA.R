library(scales)
library(car)
library(MASS)  # 多重比较
library(rstatix)  # 方差不齐anova
library(stringr)
library(reshape2)


anovaTest <- function(gene.dat, sample.info) {
  anovaTestforLine <- function(gene.dat.line) {
    gene.dat.line.melt <- melt(gene.dat.line, 
                               id.vars = "Accession", 
                               variable.name = "Group", 
                               value.name = "Abundance")
    gene.dat.line.melt <- merge(gene.dat.line.melt, 
                                sample.info, 
                                by.x = "Group", 
                                by.y = "Sample")
    gene.dat.line.melt <- subset(gene.dat.line.melt, select = -Group)
    colnames(gene.dat.line.melt)[3] <- "Group"
    
    # >>> 正态性检验 >>>
    levene.res <- leveneTest(Abundance ~ Group, data = gene.dat.line.melt, center = mean)
    levene.pval <- levene.res$`Pr(>F)`[1]
    
    # >>> T 检验 >>>
    if (levene.pval >= 0.05) {
      test.method <- "ANOVA,TukeyHSD"
      aov.model <- aov(Abundance ~ Group, gene.dat.line.melt)
      aov.summary <- summary(aov.model)[[1]]
      aov.pval <- aov.summary$`Pr(>F)`[1]
      mc <- TukeyHSD(aov.model)$Group
      mc.res <- t(as.data.frame(mc))["p adj",]
    } else {
      test.method <- "WelchANOVA,GamesHowell"
      aov.model <- welch_anova_test(gene.dat.line.melt, Abundance ~ Group)
      aov.pval <- aov.model$p
      mc <- games_howell_test(gene.dat.line.melt, Abundance ~ Group)
      mc.res <- t(data.frame(row.names = paste(mc$group2, mc$group1, sep = "-"), 
                             padj = mc$p.adj))["padj",]
    }
    
    # >>> FC 计算 >>>
    fold.change <- data.frame("a")
    sample.group <- unique(sample.info$Group)
    for (i in 1: (length(sample.group)-1)) {
      for (j in (i+1): length(sample.group)) {
        fc.colname <- paste0("FC_", sample.group[j], "_", sample.group[i])
        fold.change[fc.colname] <- mean(gene.dat.line.melt$Abundance[gene.dat.line.melt$Group == sample.group[j]]) / 
          mean(gene.dat.line.melt$Abundance[gene.dat.line.melt$Group == sample.group[i]])
      }
    }
    fold.change <- fold.change[-1]
    
    return(as.data.frame(c(gene.dat.line, Levene_p = levene.pval, 
                           ANOVA_p = aov.pval, mc.res, fold.change, 
                           TestMethod = test.method)))
  }
  
  gene.dat.res <- as.data.frame(c())
  for (i in 1: nrow(gene.dat)) {
    gene.dat.line <- gene.dat[i, ]
    test.res <- anovaTestforLine(gene.dat.line)
    gene.dat.res <- bind_rows(gene.dat.res, test.res)
  }
  
  colnames(gene.dat.res) <- gsub("X(.*)", "p_\\1", colnames(gene.dat.res))
  return(gene.dat.res)
  }