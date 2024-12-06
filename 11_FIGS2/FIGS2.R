output.dir <- "./res/293t/Norm/"
gene.dat.path <- "./data/origin_data_293t.txt"
sample.info.path <- "./data/sample_info.txt"

mycol <- c(rep("#ba3e45", 3), rep("#f2ab6a", 3), rep("#fbfa8b", 3), 
           rep("#91ad4a", 3), rep("#5ca3e8", 3), rep("#f090d0", 3))  # for draw density
filter.threshold <- 0.6  # threshold that determine whether to imputate

# >> Normalize and Imputation >>
source("./Normalization.R", local = T)
gene.dat.na <- read.delim(gene.dat.path)
sample.info <- read.delim(sample.info.path)
Normalization(output.dir, 
              gene.dat.na, 
              sample.info, 
              mycol, 
              filter.threshold)


# >>> PCA config >>>
output.dir <- "./res/293t/PCA_batch/"
gene.imputed.path <- "./res/293t/Norm/gene_data_imputation.txt"

# >> PCA >>
gene.dat.imputed <- read.delim(gene.imputed.path)
source("./PCA.R", local = T)
library(limma)
gene.dat.imputed[-1] <- log10(gene.dat.imputed[-1])
covariates <- c(rep("0 h", 3), rep("3 h", 3), rep("6 h", 3), 
                rep("9 h", 3), rep("12 h", 3), rep("24 h", 3))
gene.dat.imputed[-1] <- removeBatchEffect(gene.dat.imputed[-1], 
                                          batch = c("batch1", "batch1", "batch2", rep("batch1", 15)), 
                                          design = model.matrix(~factor(covariates))) 
gene.dat.imputed[-1] <- 10^gene.dat.imputed[-1]
PCA(gene.dat.imputed, sample.info, output.dir)
write.csv(gene.dat.imputed, 
          "./res/293t/gene_data_remove_batch_effect.csv", 
          row.names = F)

# >>> Calculate FC and p value >>>
source("./ANOVA.R", local = T)
gene.dat.imputed <- read.csv("./res/293t/gene_data_remove_batch_effect.csv")
output.dir <- "./res/293t/"
gene.dat.p <- anovaTest(gene.dat.imputed, sample.info)
write.csv(gene.dat.p, 
          paste0(output.dir, "gene_data_p.csv"), 
          row.names = F)
write.table(gene.dat.p, 
            paste0(output.dir, "gene_data_p.txt"), 
            sep = '\t', 
            quote = F, 
            row.names = F)