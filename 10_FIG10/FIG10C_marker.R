library(dplyr)

tm.data <- read.csv("./res/new/gene_data_sign.csv")
tm.marker.up <- tm.data %>% 
  filter(Group == "up") %>% 
  mutate(log2FC = log2(FC), 
         rank = (length(Group) + 1) - rank(abs(log2FC))) %>% 
  filter(rank <= 0.05 * ceiling(length(rank)))

tm.marker.down <- tm.data %>% 
  filter(Group == "down") %>% 
  mutate(log2FC = log2(FC), 
         rank = (length(Group) + 1) - rank(abs(log2FC))) %>% 
  filter(rank <= 0.05 * ceiling(length(rank)))


tm.marker <- rbind(tm.marker.down[1: 3], tm.marker.up[1: 3])
write.csv(tm.marker, "./res/tpp_marker.csv", row.names = F)
