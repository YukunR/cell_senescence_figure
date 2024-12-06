library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(NPARC)

load("./degredation.RData")
modelMetrics <- fits.degredation$metrics
fStats <- NPARCtest(modelMetrics, dfType = "empirical")
model.dat <- modelMetrics[!is.na(modelMetrics$group), ]
model.dat$tm <- log(2) / model.dat$k

model.dat <- model.dat %>%
  select(id, v, group) %>% 
  spread(key = group, value = v) %>% 
  na.omit() %>% 
  filter(S < 0  & Y < 0 ) %>% 
  melt(id.vars = "id", variable.name = "group", value.name = "v") %>% 
  mutate(v = log10(abs(v)))

mix_transform <- function(x) {
  ifelse(x <= 100, x, 100 + log10(x / 100) * 100)
}

# 混合尺度的逆转换函数
inv_mix_transform <- function(x) {
  ifelse(x <= 100, x, 10 ^ ((x - 100) / 100) * 100)
}
mixed_trans <- trans_new(name = "mixed",
                         transform = mix_transform,
                         inverse = inv_mix_transform,
                         domain = c(1, Inf))  # 设定适用域，确保没有非正数


pic <- gghistogram(model.dat, x = "v", y = "count", 
                   add = "median", rug = TRUE,
                   fill = "group", color = "group", palette = c("Y" = "#ffe00e", "S" = "#9400d3"),
                   bins = 150, 
                   
) + xlab("v (Intensity/h)") + ylab("Count")
# pic
ggsave("./res/histogram_DEG.pdf", pic, width = 8, height = 6)



fill <- c("Y" = "#fff99e", "S" = "#9400d3")
color <- c("Y" = "#ffe00e", "S" = "#9400d3")
p <- ggplot(model.dat, aes(group, v, fill = group)) +
  geom_boxplot(linewidth = 0.7, alpha = 0.8, outlier.shape = NA, ) +
  geom_line(aes(group = id), 
            color = "grey70", size = 0.5, linetype = "dashed", position = position_dodge2(0.5), 
            alpha = 0.2) +
  geom_point(aes(group = group, color = group), position = position_dodge2(0.5), size = 1, alpha = 0.7) + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  labs(x = "Group", y = "v (Intensity/h)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"),
        panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), 
        panel.background = element_blank())

p <- p + stat_compare_means(comparisons = list(c("Y", "S")),
                            method = "t.test", paired = T)
ggsave("./res/boxplot_DEG.pdf", width = 2.5, height = 4)





load("./synthesis.RData")
modelMetrics <- fits.synthesis$metrics
fStats <- NPARCtest(modelMetrics, dfType = "empirical")
model.dat <- modelMetrics[!is.na(modelMetrics$group), ]
model.dat$tm <- log(2) / model.dat$k

model.dat <- model.dat %>%
  select(id, v, group) %>% 
  spread(key = group, value = v) %>% 
  na.omit() %>% 
  filter(S >= 0 & Y >= 0) %>% 
  melt(id.vars = "id", variable.name = "group", value.name = "v") %>% 
  mutate(v = log10(abs(v)))

mix_transform <- function(x) {
  ifelse(x <= 100, x, 100 + log10(x / 100) * 100)
}

# 混合尺度的逆转换函数
inv_mix_transform <- function(x) {
  ifelse(x <= 100, x, 10 ^ ((x - 100) / 100) * 100)
}
mixed_trans <- trans_new(name = "mixed",
                         transform = mix_transform,
                         inverse = inv_mix_transform,
                         domain = c(1, Inf))  # 设定适用域，确保没有非正数

pic <- gghistogram(model.dat, x = "v", y = "count", 
                   add = "median", rug = TRUE,
                   fill = "group", color = "group", palette = c("Y" = "#ffe00e", "S" = "#9400d3"),
                   bins = 150
) + xlab("v (Intensity/h)") + ylab("Count")


# pic
ggsave("./res/histogram_SYN.pdf", pic, width = 8, height = 6)



fill <- c("Y" = "#fff99e", "S" = "#9400d3")
color <- c("Y" = "#ffe00e", "S" = "#9400d3")
p <- ggplot(model.dat, aes(group, v, fill = group)) +
  geom_boxplot(linewidth = 0.7, alpha = 0.8, outlier.shape = NA, ) +
  geom_line(aes(group = id), 
            color = "grey70", size = 0.5, linetype = "dashed", position = position_dodge2(0.5), 
            alpha = 0.2) +
  geom_point(aes(group = group, color = group), position = position_dodge2(0.5), size = 1, alpha = 0.7) + 
  scale_fill_manual(values = fill) + 
  scale_color_manual(values = color) + 
  labs(x = "Group", y = "v (Intensity/h)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"),
        panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), 
        panel.background = element_blank())

p <- p + stat_compare_means(comparisons = list(c("Y", "S")),
                            method = "t.test", paired = T)
ggsave("./res/boxplot_SYN.pdf", width = 2.5, height = 4)
