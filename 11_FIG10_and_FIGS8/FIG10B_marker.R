library(dplyr)

tm.data <- read.csv("./res/FIG3/tm/gene_data_sign_deg.csv")
tm.marker.up <- tm.data %>% filter(Group == "up") %>% 
  mutate(rank = ((length(log2FC)+1) - rank(abs(log2FC)))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
tm.marker.down <- tm.data %>% filter(Group == "down") %>% 
  mutate(rank = ((length(log2FC)+1) - rank(abs(log2FC)))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
tm.marker <- unique(rbind(tm.marker.up[1: 3], tm.marker.down[1: 3]))

v.deg.data <- read.csv("./res/FIG3/v/gene_data_sign_deg.csv")%>% 
  filter(!Group == "no")
v.syn.data <- read.csv("./res/FIG3/v/gene_data_sign_syn.csv")%>% 
  filter(!Group == "no")


v.deg.marker.up <- v.deg.data %>% 
  filter(Group == "up") %>% 
  mutate(rank = (length(log2FC)+1) - rank(abs(log2FC))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
v.deg.marker.down <- v.deg.data %>% 
  filter(Group == "down") %>% 
  mutate(rank = (length(log2FC)+1) - rank(abs(log2FC))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
v.deg.marker <- rbind(v.deg.marker.down[1: 3], v.deg.marker.up[1: 3])

v.syn.marker.up <- v.syn.data %>%  
  filter(Group == "up") %>% 
  mutate(rank = (length(log2FC)+1) - rank(abs(log2FC))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
v.syn.marker.down <- v.syn.data %>%  
  filter(Group == "down") %>% 
  mutate(rank = (length(log2FC)+1) - rank(abs(log2FC))) %>% 
  filter(rank <= ceiling(0.05 * length(rank)))
v.syn.marker <- rbind(v.syn.marker.down[1: 3], v.syn.marker.up[1: 3])

v.marker <- unique(rbind(v.deg.marker[1: 3], v.syn.marker[1: 3]))


write.csv(marker, "./res/marker.csv", row.names = F) 
