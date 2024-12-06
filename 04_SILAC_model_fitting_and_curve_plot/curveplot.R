# =====================================
# ===== include fig9A, 9E, 10B ========
# =====================================


library(NPARC)
library(dplyr)
library(tidyr)
library(ggplot2)
load("./degredation.RData")
load("./synthesis.RData")
load("./data_df.RData")


# >> read diff protien
gene.dat <- read.csv("./res/FIG3/tm/gene_data_sign_deg.csv")
gene.dat <- gene.dat %>% filter(!Group == "no")

model.matrix.degredation <- fits.degredation$metrics
fStats.degredation <- NPARCtest(model.matrix.degredation, dfType = "empirical")

model.matrix.synthesis <- fits.synthesis$metrics
fStats.synthesis <- NPARCtest(model.matrix.synthesis, dfType = "empirical")

width <- options()$width #获取显示界面行宽度
N <- nrow(gene.dat)

for (i in 1: nrow(gene.dat)) {
  gene <- gene.dat$Accession[i]
  param.synthesis <- model.matrix.synthesis[model.matrix.synthesis$id == gene,]
  data.synthesis <- df.synthesis[df.synthesis$Accession == gene, ]
  
  param.degredation <- model.matrix.degredation[model.matrix.degredation$id == gene,]
  data.degredation <- df.degredation[df.degredation$Accession == gene,]
  # param.synthesis$a = param.synthesis$b = c(6e8, 6e8, 6e8)
  data.degredation <- data.degredation %>%
    left_join(param.degredation, by = c("Accession" = "id", "group" = "group"))
  data.synthesis <- data.synthesis %>%
    left_join(param.synthesis, by = c("Accession" = "id", "group" = "group"))
  
  curve_data1 <- param.degredation %>%
    crossing(x = seq(0, 48)) %>%
    mutate(y = a * exp(-k * x), model = "Degradation")
  
  curve_data2 <- param.synthesis %>%
    crossing(x = seq(0, 48)) %>%
    mutate(y = b - a * exp(-k * x), model = "Synthesis")
  
  curve_data <- bind_rows(curve_data1, curve_data2)
  # 绘制散点图和拟合曲线
  curve_data_m1_y <- curve_data %>% filter(group == "Y", model == "Degradation")
  curve_data_m1_s <- curve_data %>% filter(group == "S", model == "Degradation")
  # Model 2
  curve_data_m2_y <- curve_data %>% filter(group == "Y", model == "Synthesis")
  curve_data_m2_s <- curve_data %>% filter(group == "S", model == "Synthesis")
  
  curve_data_y <- curve_data_m1_y %>%
    left_join(curve_data_m2_y, by = "x") %>%
    mutate(y = y.x + y.y, group = "Y", model = "Sum")
  
  
  curve_data_s <- curve_data_m1_s %>%
    left_join(curve_data_m2_s, by = "x") %>%
    mutate(y = y.x + y.y, group = "S", model = "Sum")
  
  # 合并所有曲线数据
  curve_data_final <- bind_rows(curve_data, curve_data_y, curve_data_s) %>% 
    mutate(Group = model, 
           group = ifelse(is.na(group), "Null fit", ifelse(group=="Y", "Young", "Senescence")))
  
  # 合并原始数据（如果在同一图中绘制散点图）
  data <- bind_rows(data.degredation, data.synthesis) %>% 
    mutate(Model = ifelse(isotope == "L", "Degradation", "Synthesis"), 
           Color = ifelse(group=="Y", "Young", "Senescence"))
  
  p <- ggplot() + 
    geom_point(data = data, aes(x = time, y = Abundance, color = Color, shape = Model), alpha = 0.6) +
    geom_line(data = curve_data_final, aes(x = x, y = y, color = group, linetype = Group), size = 1) +
    labs(title = paste0("Synthesis and degradation of ", gene.dat$GeneName[i], "\nHalf life: Young (", 
                        round(gene.dat$Y[i], 2), " h), Senescence (", round(gene.dat$S[i], 2), " h)"), 
         x = "Time (h)", y = "Abundance") +
    scale_color_manual(values = c("Young" = "#ffe00e", "Senescence" = "#9400d3", "Null fit"="grey70")) + # 确保加和组有独特颜色
    scale_linetype_manual(values = c("Degradation" = "dashed", "Synthesis" = "dotdash", "Sum" = "solid")) + # 使用不同线型区分模型
    theme_minimal() + 
    theme_bw()
  p
  
  ggsave(paste0("./res/NPARC/curve/", gene.dat$Accession[i], "_", gene.dat$GeneName[i], ".pdf"), width = 6, height = 4)
  cat('[', paste0(rep('#', i/N*width), collapse=''),
      paste0(rep('-', width - i/N*width), collapse=''),
      ']',
      round(i/N*100),'%')
  if(i==N)cat('\nDONE!\n')
  else cat('\r')
}

