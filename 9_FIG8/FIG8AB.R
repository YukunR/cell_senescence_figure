library(dplyr)

gene.data.quantity <- read.csv("../cellSenescence/res/lo2/24_0h/gene_data_sign.csv")
gene.data.thalf <- read.csv("../CellSenescenceSilacNew/res/FIG3/tm/gene_data_sign_deg.csv")
gene.data.tm <- read.csv("../CellSenescenceTPP/res/new/gene_data_sign.csv")

combined_data <- inner_join(gene.data.quantity, gene.data.thalf, by = "Accession")


# 计算概率分布
combined_data$mean_S0 <- log10(rowMeans(combined_data[, c("S0.1", "S0.2", "S0.3")]))
combined_data$mean_S24 <- log10(rowMeans(combined_data[, c("S24.1", "S24.2", "S24.3")]))
combined_data$mean_Y <- log2(combined_data$Y)
combined_data$mean_S <- log2(combined_data$S)

normalize_zscore <- function(x) {
  (x - mean(x)) / sd(x)
}

combined_data$mean_S0_norm <- normalize_zscore(combined_data$mean_S0)
combined_data$mean_S24_norm <- normalize_zscore(combined_data$mean_S24)
combined_data$mean_Y_norm <- normalize_zscore(combined_data$mean_Y)
combined_data$mean_S_norm <- normalize_zscore(combined_data$mean_S)

# 离散化数据
combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)


combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)

# 使用entropy包计算熵和互信息
library(entropy)

# 计算概率分布和熵
p_S0_Y <- table(combined_data$S0_binned, combined_data$Y_binned) / nrow(combined_data)
p_S24_S <- table(combined_data$S24_binned, combined_data$S_binned) / nrow(combined_data)

H_S0_Y <- entropy(as.vector(p_S0_Y), unit="log2")
H_S24_S <- entropy(as.vector(p_S24_S), unit="log2")

# 计算互信息
I_S0_Y <- entropy(table(combined_data$S0_binned), unit="log2") + entropy(table(combined_data$Y_binned), unit="log2") - H_S0_Y
I_S24_S <- entropy(table(combined_data$S24_binned), unit="log2") + entropy(table(combined_data$S_binned), unit="log2") - H_S24_S

# 输出结果
cat("互信息 (S0 vs Y):", I_S0_Y, "\n")
cat("互信息 (S24 vs S):", I_S24_S, "\n")


library(ggplot2)
library(ggdensity)
library(reshape2)
library(ggpmisc)

plot.data <- combined_data %>% 
  select(Accession, mean_S0_norm, mean_S24_norm, mean_Y_norm, mean_S_norm)
plot.data1 <- plot.data %>%
  select(c(1, 2, 4)) %>% 
  mutate(Group = "Young")
colnames(plot.data1) <- c("Accession", "Quantity", "T half", "Group")
plot.data2 <- plot.data %>%
  select(c(1, 3, 5)) %>% 
  mutate(Group = "Senescence")
colnames(plot.data2) <- c("Accession", "Quantity", "T half", "Group")
plot.data <- rbind(plot.data1, plot.data2) %>% 
  mutate(Group = factor(Group, levels = c("Young", "Senescence")))

fill <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
color <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
p <- ggplot(plot.data, aes(x=Quantity, y=`T half`, color = Group, fill = Group)) + 
  ggtitle(paste0(I_S0_Y, "   ", I_S24_S)) + 
  # stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + 
  geom_hdr_rug() +
  geom_hdr_lines(linewidth = 0.6) +
  geom_point(alpha = 0.2, shape = ".", size = 0.5) + 
  # labs(x = "quantity", y = "tm") + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "right") + 
  facet_grid(~Group) + 
  stat_smooth(color = "#1868b2", formula = y ~ x,fill = "#1868b2", method = "lm") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 3, # 公式字体大小
    label.x = 0.1,
    label.y = 0.95)
p


ggsave("./res/FIG5/thalf_quantity.pdf", p, width = 8, height = 4)


# ================================================================
# ================================================================
# ================================================================
# ================================================================
# ================================================================

combined_data <- inner_join(gene.data.quantity, gene.data.tm, by = "Accession")


# 计算概率分布
combined_data$mean_S0 <- log10(rowMeans(combined_data[, c("S0.1", "S0.2", "S0.3")]))
combined_data$mean_S24 <- log10(rowMeans(combined_data[, c("S24.1", "S24.2", "S24.3")]))
combined_data$mean_Y <- log2(combined_data$YOU)
combined_data$mean_S <- log2(combined_data$SEN)

normalize_zscore <- function(x) {
  (x - mean(x)) / sd(x)
}

combined_data$mean_S0_norm <- normalize_zscore(combined_data$mean_S0)
combined_data$mean_S24_norm <- normalize_zscore(combined_data$mean_S24)
combined_data$mean_Y_norm <- normalize_zscore(combined_data$mean_Y)
combined_data$mean_S_norm <- normalize_zscore(combined_data$mean_S)

# 离散化数据
combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)


combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)

# 使用entropy包计算熵和互信息
library(entropy)

# 计算概率分布和熵
p_S0_Y <- table(combined_data$S0_binned, combined_data$Y_binned) / nrow(combined_data)
p_S24_S <- table(combined_data$S24_binned, combined_data$S_binned) / nrow(combined_data)

H_S0_Y <- entropy(as.vector(p_S0_Y), unit="log2")
H_S24_S <- entropy(as.vector(p_S24_S), unit="log2")

# 计算互信息
I_S0_Y <- entropy(table(combined_data$S0_binned), unit="log2") + entropy(table(combined_data$Y_binned), unit="log2") - H_S0_Y
I_S24_S <- entropy(table(combined_data$S24_binned), unit="log2") + entropy(table(combined_data$S_binned), unit="log2") - H_S24_S

# 输出结果
cat("互信息 (S0 vs Y):", I_S0_Y, "\n")
cat("互信息 (S24 vs S):", I_S24_S, "\n")


library(ggplot2)
library(ggdensity)
library(reshape2)
library(ggExtra)

plot.data <- combined_data %>% 
  select(Accession, mean_S0_norm, mean_S24_norm, mean_Y_norm, mean_S_norm)
plot.data1 <- plot.data %>%
  select(c(1, 2, 4)) %>% 
  mutate(Group = "Young")
colnames(plot.data1) <- c("Accession", "Quantity", "Tm", "Group")
plot.data2 <- plot.data %>%
  select(c(1, 3, 5)) %>% 
  mutate(Group = "Senescence")
colnames(plot.data2) <- c("Accession", "Quantity", "Tm", "Group")
plot.data <- rbind(plot.data1, plot.data2) %>% 
  mutate(Group = factor(Group, levels = c("Young", "Senescence")))

fill <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
color <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
p <- ggplot(plot.data, aes(x=Quantity, y=`Tm`, color = Group, fill = Group)) + 
  ggtitle(paste0(I_S0_Y, "   ", I_S24_S)) + 
  # stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + 
  geom_hdr_rug() +
  geom_hdr_lines(linewidth = 0.6) +
  geom_point(alpha = 0.2, shape = ".", size = 0.5) + 
  # labs(x = "quantity", y = "tm") + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "right") + 
  facet_grid(~Group) + 
  stat_smooth(color = "#1868b2", formula = y ~ x,fill = "#1868b2", method = "lm") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 3, # 公式字体大小
    label.x = 0.1,
    label.y = 0.95) + 
  coord_cartesian(ylim = c(-5, 7.5))
p
ggsave("./res/FIG5/tm_quantity.pdf", p, width = 8, height = 4)



# TODO: 假如互信息，生成第三幅图
# ================================================================
# ================================================================
# ================================================================
# ================================================================
# ================================================================

combined_data <- inner_join(gene.data.thalf, gene.data.tm, by = "Accession")


# 计算概率分布
combined_data$mean_Y <- log2(combined_data$Y)
combined_data$mean_S <- log2(combined_data$S)
combined_data$mean_S0 <- log2(combined_data$YOU)
combined_data$mean_S24 <- log2(combined_data$SEN)

normalize_zscore <- function(x) {
  (x - mean(x)) / sd(x)
}

combined_data$mean_S0_norm <- normalize_zscore(combined_data$mean_S0)
combined_data$mean_S24_norm <- normalize_zscore(combined_data$mean_S24)
combined_data$mean_Y_norm <- normalize_zscore(combined_data$mean_Y)
combined_data$mean_S_norm <- normalize_zscore(combined_data$mean_S)

# 离散化数据
combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)

# 使用entropy包计算熵和互信息
library(entropy)

# 计算概率分布和熵
p_S0_Y <- table(combined_data$S0_binned, combined_data$Y_binned) / nrow(combined_data)
p_S24_S <- table(combined_data$S24_binned, combined_data$S_binned) / nrow(combined_data)

H_S0_Y <- entropy(as.vector(p_S0_Y), unit="log2")
H_S24_S <- entropy(as.vector(p_S24_S), unit="log2")

# 计算互信息
I_S0_Y <- entropy(table(combined_data$S0_binned), unit="log2") + entropy(table(combined_data$Y_binned), unit="log2") - H_S0_Y
I_S24_S <- entropy(table(combined_data$S24_binned), unit="log2") + entropy(table(combined_data$S_binned), unit="log2") - H_S24_S

# 输出结果
cat("互信息 (S0 vs Y):", I_S0_Y, "\n")
cat("互信息 (S24 vs S):", I_S24_S, "\n")


library(ggplot2)
library(ggdensity)
library(reshape2)
library(ggExtra)

plot.data <- combined_data %>% 
  select(Accession, mean_S0_norm, mean_S24_norm, mean_Y_norm, mean_S_norm)
plot.data1 <- plot.data %>%
  select(c(1, 2, 4)) %>% 
  mutate(Group = "Young")
colnames(plot.data1) <- c("Accession", "Tm", "T half", "Group")
plot.data2 <- plot.data %>%
  select(c(1, 3, 5)) %>% 
  mutate(Group = "Senescence")
colnames(plot.data2) <- c("Accession", "Tm", "T half", "Group")
plot.data <- rbind(plot.data1, plot.data2) %>% 
  mutate(Group = factor(Group, levels = c("Young", "Senescence")))

fill <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
color <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
p <- ggplot(plot.data, aes(x=`T half`, y=`Tm`, color = Group, fill = Group)) + 
  ggtitle(paste0(I_S0_Y, "   ", I_S24_S)) + 
  # stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + 
  geom_hdr_rug() +
  geom_hdr_lines(linewidth = 0.6) +
  geom_point(alpha = 0.2, shape = ".", size = 0.5) + 
  # labs(x = "quantity", y = "tm") + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "right") + 
  facet_grid(~Group) + 
  stat_smooth(color = "#1868b2", formula = y ~ x,fill = "#1868b2", method = "lm") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 3, # 公式字体大小
    label.x = 0.1,
    label.y = 0.95) + 
  coord_cartesian(ylim = c(-4, 5))
ggsave("./res/FIG5/tm_thalf.pdf", p, width = 8, height = 4)



# ================================================================
# ================================================================
# ================================================================
# ================================================================
# ================================================================

gene.data.quantity.lo2 <- read.csv("../cellSenescence/res/lo2/24_0h/gene_data_sign.csv")
gene.data.quantity.293t <- read.csv("../cellSenescence/res/293t/gene_data_p.csv")
combined_data <- merge(gene.data.quantity.lo2, gene.data.quantity.293t, by = "Accession", suffixes = c(".lo2", ".293t"))

library(limma)
combined_data <- combined_data %>% 
  select(Accession, S0.1.lo2, S0.2.lo2, S0.3.lo2, S24.1.lo2, S24.2.lo2, S24.3.lo2, 
         S0.1.293t, S0.2.293t, S0.3.293t, S24.1.293t, S24.2.293t, S24.3.293t)
combined_data[-1] <- removeBatchEffect(log10(combined_data[-1]), batch = c(rep("batch1", 6), rep("batch2", 6)))
combined_data[-1] <- 10^combined_data[-1]
# 计算概率分布
combined_data$mean_S0 <- log10(rowMeans(combined_data[, c("S0.1.lo2", "S0.2.lo2", "S0.3.lo2")]))
combined_data$mean_S24 <- log10(rowMeans(combined_data[, c("S24.1.lo2", "S24.2.lo2", "S24.3.lo2")]))
combined_data$mean_Y <- log10(rowMeans(combined_data[, c("S0.1.293t", "S0.2.293t", "S0.3.293t")]))
combined_data$mean_S <- log10(rowMeans(combined_data[, c("S24.1.293t", "S24.2.293t", "S24.3.293t")]))

normalize_zscore <- function(x) {
  (x - mean(x)) / sd(x)
}

combined_data$mean_S0_norm <- normalize_zscore(combined_data$mean_S0)
combined_data$mean_S24_norm <- normalize_zscore(combined_data$mean_S24)
combined_data$mean_Y_norm <- normalize_zscore(combined_data$mean_Y)
combined_data$mean_S_norm <- normalize_zscore(combined_data$mean_S)

# 离散化数据
combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)


combined_data$S0_binned <- cut(combined_data$mean_S0_norm, breaks=4, labels=FALSE)
combined_data$S24_binned <- cut(combined_data$mean_S24_norm, breaks=4, labels=FALSE)
combined_data$Y_binned <- cut(combined_data$mean_Y_norm, breaks=4, labels=FALSE)
combined_data$S_binned <- cut(combined_data$mean_S_norm, breaks=4, labels=FALSE)

# 使用entropy包计算熵和互信息
library(entropy)

# 计算概率分布和熵
p_S0_Y <- table(combined_data$S0_binned, combined_data$Y_binned) / nrow(combined_data)
p_S24_S <- table(combined_data$S24_binned, combined_data$S_binned) / nrow(combined_data)

H_S0_Y <- entropy(as.vector(p_S0_Y), unit="log2")
H_S24_S <- entropy(as.vector(p_S24_S), unit="log2")

# 计算互信息
I_S0_Y <- entropy(table(combined_data$S0_binned), unit="log2") + entropy(table(combined_data$Y_binned), unit="log2") - H_S0_Y
I_S24_S <- entropy(table(combined_data$S24_binned), unit="log2") + entropy(table(combined_data$S_binned), unit="log2") - H_S24_S

# 输出结果
cat("互信息 (S0 vs Y):", I_S0_Y, "\n")
cat("互信息 (S24 vs S):", I_S24_S, "\n")


library(ggplot2)
library(ggdensity)
library(reshape2)
library(ggExtra)

plot.data <- combined_data %>% 
  select(Accession, mean_S0_norm, mean_S24_norm, mean_Y_norm, mean_S_norm)
plot.data1 <- plot.data %>%
  select(c(1, 2, 4)) %>% 
  mutate(Group = "Young")
colnames(plot.data1) <- c("Accession", "Quantity.LO2", "Quantity.293T", "Group")
plot.data2 <- plot.data %>%
  select(c(1, 3, 5)) %>% 
  mutate(Group = "Senescence")
colnames(plot.data2) <- c("Accession", "Quantity.LO2", "Quantity.293T", "Group")
plot.data <- rbind(plot.data1, plot.data2) %>% 
  mutate(Group = factor(Group, levels = c("Young", "Senescence")))

fill <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
color <- c("Young" = "#ffc00e", "Senescence" = "#9400d3")
p <- ggplot(plot.data, aes(x=`Quantity.LO2`, y=`Quantity.293T`, color = Group, fill = Group)) + 
  ggtitle(paste0(I_S0_Y, "   ", I_S24_S)) + 
  # stat_density_2d(aes(fill = stat(level)), geom = 'polygon') + 
  geom_hdr_rug() +
  geom_hdr_lines(linewidth = 0.6) +
  geom_point(alpha = 0.2, shape = ".", size = 0.5) + 
  # labs(x = "quantity", y = "tm") + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "right") + 
  facet_grid(~Group) + 
  stat_smooth(color = "#1868b2", formula = y ~ x,fill = "#1868b2", method = "lm") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 3, # 公式字体大小
    label.x = 0.1,
    label.y = 0.95)
ggsave("./res/FIG5/lo2_293t_remove_batch.pdf", p, width = 8, height = 4)