# >>> 时间序列分析 >>>
library(Mfuzz)
library(dplyr)
gene.dat.p <- read.csv("./res/293t/gene_data_p.csv")
gene.dat.time <- gene.dat.p %>% filter(ANOVA_p < 0.05) %>% select(1: 19)
group.names <- unique(sample.info$Group)
for (group.name in group.names) {
  columns <- names(gene.dat.time) %in% c(sample.info$Sample[sample.info$Group == group.name])
  gene.dat.time[[group.name]] <- rowMeans(gene.dat.time[columns])
}
gene.dat.time <- gene.dat.time %>% 
  select(all_of(c("Accession", group.names)))
row.names(gene.dat.time) <- gene.dat.time$Accession
gene.dat.time <- as.matrix(select(gene.dat.time, -Accession))

mfuzz.class <- new('ExpressionSet', exprs = gene.dat.time)
mfuzz.class <- standardise(mfuzz.class)

cluster.num <- 4
set.seed(2)
mfuzz.cluster <- mfuzz(mfuzz.class, c = cluster.num, m = mestimate(mfuzz.class))
library(RColorBrewer)
color <- colorRampPalette(rev(c("#df6561", "#f2d086", "#95be7e", "#57b1a9")))(1000)
pdf("./res/293t/mFuzz/cluster.pdf", width = 4, height = 10)
mfuzz.plot2(mfuzz.class, cl = mfuzz.cluster, mfrow = c(4, 1), time.labels = colnames(gene.dat.time)[-7], 
            x11 = F, colo = "fancy", centre = T, Xwidth = 4, Xheight = 10)
dev.off()

png("./res/293t/mFuzz/cluster.png", width = 4000, height = 2500, res = 300)
# mfuzz.plot(mfuzz.class, cl = mfuzz.cluster, mfrow = c(3, 3), time.labels = colnames(gene.dat.time), 
#            new.window = F, colo = color)
mfuzz.plot2(mfuzz.class, cl = mfuzz.cluster, mfrow = c(3, 3), time.labels = colnames(gene.dat.time)[-7], 
            x11 = F, colo = "fancy", centre = T)

dev.off()

gene.cluster <- data.frame(cluster = mfuzz.cluster$cluster)
gene.dat.time <- as.data.frame(gene.dat.time)
gene.dat.time$Accession <- row.names(gene.dat.time)
gene.cluster$Accession <- row.names(gene.cluster)
gene.dat.cluster <- merge(gene.dat.p, gene.cluster, by = "Accession")
write.csv(gene.dat.cluster, "./res/293t/mFuzz/gene_dat_cluster.csv", row.names = F)