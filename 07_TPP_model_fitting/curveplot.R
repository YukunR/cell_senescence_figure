# =====================================
# ===== include fig9B, 9F, 10C ========
# =====================================

library(dplyr)
library(ggplot2)
library(NPARC)
library(tidyr)

load("./NPARCFit.RData")
load("./data_df.RData")

# >> 读取数据
gene.dat <- read.csv("./res/new/gene_data_sign.csv") %>% filter(!Group == "no")
gene.dat$GeneName <- gsub(".*GN=([^ ]+).*", "\\1", gene.dat$Description)

width <- options()$width #获取显示界面行宽度
N <- nrow(gene.dat)

for (i in 1: nrow(gene.dat)) {
  
  gene <- gene.dat$Accession[i]
  param <- modelMetrics[modelMetrics$id == gene,]
  data.df <- df[df$Gene == gene, ]
  
  # param.synthesis$a = param.synthesis$b = c(6e8, 6e8, 6e8)
  
  data.df <- data.df %>%
    left_join(param, by = c("Gene" = "id", "Conc" = "group"))
  
  curve_data <- param %>%
    crossing(x = seq(37, 67)) %>%
    mutate(y = (1 - pl)  / (1 + exp((b - a/x))) + pl) %>% 
    mutate(group = ifelse(is.na(group), "Null fit", group))
  
  data.df <- data.df %>% mutate(Group = Conc)
  p <- ggplot() + 
    geom_point(data = data.df, aes(x = Temp, y = RelAbundance, color = Group, shape = Group), alpha = 0.6) +
    geom_line(data = curve_data, aes(x = x, y = y, color = group, linetype=group), size = 1) +
    labs(title = paste0("Melting curve of ", gene.dat$GeneName[i], "\nTm: Young (", 
                        round(gene.dat$YOU[i], 2), " °C), Senescence (", round(gene.dat$SEN[i], 2), " °C)"), 
         x = "Temperature (°C)", y = "Non-denatured fraction") +
    scale_color_manual(values = c("YOU" = "#ffe00e", "SEN" = "#9400d3")) + # 确保加和组有独特颜色
    scale_linetype_manual(values = c("YOU" = "solid", "SEN" = "solid", "Null fit" = "dashed")) + # 使用不同线型区分模型
    theme_minimal() + 
    theme_bw()
  
  ggsave(paste0("./res/melting_curve/", gene.dat$Accession[i], "_", gene.dat$GeneName[i], ".pdf"), width = 6, height = 4)
  cat('[', paste0(rep('#', i/N*width), collapse=''),
      paste0(rep('-', width - i/N*width), collapse=''),
      ']',
      round(i/N*100),'%')
  if(i==N)cat('\nDONE!\n')
  else cat('\r')
} 
