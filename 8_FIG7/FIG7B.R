library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

load("NPARCFit.RData")
model.dat <- modelMetrics[!is.na(modelMetrics$group), ]
model.dat <- model.dat %>%
  select(id, tm, group) %>% 
  spread(key = group, value = tm) %>% 
  na.omit() %>% 
  melt(id.vars = "id", variable.name = "group", value.name = "tm")
model.dat <- model.dat[model.dat$tm > 0, ]
pic <- gghistogram(model.dat, x = "tm", y = "count", 
                   add = "mean", rug = TRUE,
                   fill = "group", color = "group", palette = c("YOU" = "#ffe00e", "SEN" = "#9400d3"),
                   bins = 80, 
                   add_density = TRUE # 添加密度曲线
) + scale_x_continuous(limits = c(40, 90), breaks = seq(40, 90, 5))
pic
ggsave("./res/histogram_average.pdf", pic, width = 6, height = 4)

# box plot
library(ggpubr)
library(rstatix)
# 单边pair t检验
stat.test <- model.dat %>%
  t_test(tm ~ group, paired = T) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(x = "group")
bxp <- ggpaired(model.dat, x = "group", y = "tm", fill = "group", line.color = "gray", line.size = 0.4, 
                point.size = 0.8) + 
  scale_fill_manual(values = c("YOU" = "#ffe00e", "SEN" = "#9400d3"))
bxp <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
bxp
ggsave("./res/boxplot.pdf", bxp, width = 3, height = 5)