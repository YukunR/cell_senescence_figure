
library(dplyr)
library(reshape2)
library(stringr)


gene.data <- read.csv("./res/norm_data.csv")
gene.data.colsum <- colSums(gene.data[-c(1: 2)], na.rm = T)
names(gene.data.colsum) <- colnames(gene.data[-c(1: 2)])
gene.data.colsum <- as.data.frame(gene.data.colsum) %>% 
  mutate(Group = row.names(.)) %>% 
  rename(colsum = gene.data.colsum) %>% 
  mutate(rep = gsub(".*_(\\d+)", "\\1", Group), 
         isotope = gsub("[YS]\\d+(\\w+)_.*", "\\1", Group), 
         time = paste(gsub("[YS](\\d+).*", "\\1", Group), "h", sep = " "), 
         group = gsub("([YS]).*", "\\1", Group), 
         group = factor(group, levels = c("Y", "S")), 
         Group = str_split_i(Group, "_", 1), 
         time = factor(time, levels = paste(c(1, 3, 6, 9, 24), "h", sep = " "))) %>% 
  group_by(time, group, rep) %>% 
  mutate(sum_ = sum(colsum)) %>% 
  ungroup() %>% 
  mutate(relcolsum = colsum / sum_)


gene.data.colsum.mean <- gene.data.colsum %>% 
  group_by(Group) %>% 
  mutate(mean.colsum = mean(relcolsum)) %>% 
  ungroup() %>%
  select(-c(Group, rep, colsum, relcolsum, sum_)) %>% 
  unique()

library(ggalluvial)
library(ggplot2)
ggplot(gene.data.colsum.mean, aes(x = time, y = mean.colsum, fill = isotope)) +
  geom_stratum(aes(stratum = isotope), color = NA, width = 0.65) +
  geom_flow(aes(alluvium = isotope), knot.pos = 0.25, width = 0.65, alpha = 0.85) +
  geom_alluvium(aes(alluvium = isotope),
                knot.pos = 0.25, color = "white", width = 0.65, linewidth = 0.5, fill = NA, alpha = 1) + 
  facet_wrap(~ group, nrow = 1, ncol = 2) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_fill_manual(values = c("H" = "#4758a2", "L" = "#fcab58"), 
                    labels = c("H" = "Synthesis", "L" = "Degredation")) + 
  labs(x = "Time", y = "Relative intensity", fill = "Metabolism process") + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        legend.position = "top")

ggsave("./res/total_express.pdf", width = 4.5, height = 3)
