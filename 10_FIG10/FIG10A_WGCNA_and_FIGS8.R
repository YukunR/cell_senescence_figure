library(WGCNA)
library(dplyr)
enableWGCNAThreads(nThreads=6)
colSumNormalize <- function(df) {
  col_sums <- colSums(df)
  # 计算所有列总和的平均值
  average_col_sum <- mean(col_sums)
  # 按列归一化，使每列的和为列总和的平均值
  normalized_df <- sweep(df, 2, col_sums, FUN="/") * average_col_sum
}


gene.data.293t <- read.csv("./res/293t/Norm/gene_data_imputation.csv")
gene.data.lo2 <- read.csv("./res/lo2/gene_data_remove_batch_effect.csv")
gene.data.all <- merge(gene.data.293t, gene.data.lo2, by = "Accession", suffixes = c(".293t", ".lo2"))
row.names(gene.data.all) <- gene.data.all$Accession
gene.data.all <- subset(gene.data.all, select = -Accession)

gene.data.norm <- colSumNormalize(gene.data.all)


sample.info <- read.delim("./data/sample_info_all.txt")
gene.dat <- na.omit(gene.data.norm)


calculate_cv <- function(x) {
  sd(x) / mean(x) * 100
}
# 初始化一个逻辑向量，用于标记保留的行
keep_rows <- rep(TRUE, nrow(gene.dat))
# 遍历每个组，计算CV并更新保留的行
for (group in unique(sample.info$Group)) {
  # 构造该组内所有列的名称
  columns <- sample.info$Sample[sample.info$Group == group]
  
  # 计算该组的CV
  group_cv <- apply(gene.dat[, columns], 1, calculate_cv)
  
  # 更新保留的行
  keep_rows <- keep_rows & (group_cv < 40)
}
# 筛选数据
gene.dat <- gene.dat[keep_rows, ]


# remove batch effect
library(limma)
batch.info <- c(rep("batch1", 18), rep("batch2", 18))
gene.dat <- log10(gene.dat)
gene.dat <- removeBatchEffect(gene.dat, batch = batch.info)
# gene.dat <- as.data.frame(10^gene.dat)
datExpr <- t(gene.dat)

# ----------------- soft threshold -----------------------
#soft thresholding power值即软阈值的筛选原则是使构建的网络更符合无标度网络特征 
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# ---------------------- non-scale network -------------------
net = blockwiseModules(datExpr, power = sft$powerEstimate, # 3
                       corType = "pearson",#or bicor
                       TOMType = "unsigned", # minModuleSize = 30, #30
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)
mergedColors = labels2colors(net$colors)
pdf("./res/WGCNA/dendrogram.pdf", width = 7, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE
                    , guideHang = 0.05)
dev.off()
# -------------------- Module trait ----------------------
nSamples <- nrow(datExpr)
datTraits <- sample.info %>% 
  mutate(Time = gsub(".*_(\\d+) h", "\\1", Group)) %>% 
  select(-Group)
row.names(datTraits) <- datTraits$Sample
datTraits <- datTraits[-1]

MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf(file = "./res/WGCNA/Module-trait relationships.pdf", width = 4, height = 7)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# ------------------- gene ---------------------
############################################################
################## select the trait that your interests 
############################################################

Time = as.data.frame(datTraits$Time);
names(Time) = "Time"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Time, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Time), sep="");
names(GSPvalue) = paste("p.GS.", names(Time), sep="");

#####module选择对应显著的module
module = "pink"
column = match(module, modNames);
moduleGenes = mergedColors==module;
pdf(file = "./res/WGCNA/pinkScatterplot.pdf", width = 12, height = 9);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Time",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# Select module probes 选择模块对应的基因
probes = colnames(datExpr)
inModule = (mergedColors==module);
modProbes = probes[inModule];
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= ceiling(0.05 * sum(inModule)))  # 取前5%
module.out.pink.top <- module.out.pink[top, ]
