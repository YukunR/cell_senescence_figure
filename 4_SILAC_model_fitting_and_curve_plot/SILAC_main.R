
# -------------------------- 拟合 -----------------------------------
library(dplyr)
library(NPARC)
library(stringr)
library(reshape2)
library(broom)
library(ggplot2)
library(tidyr)
source("./evalModels.R", local = T)
source("./modelFitting.R", local = T)
# >>> config >>>

data.file <- "./res/norm_data.csv"
output.dir <- "./res/NPARC/"

# ------------ do not change code blow -------------------------------
# >>> set output dir >>>
output.dir.hierarchy <- str_split(output.dir, "/")[[1]]
for (i in 2: length(output.dir.hierarchy)) {
  dir.tmp <- paste(output.dir.hierarchy[1: i], collapse = "/")
  if (!dir.exists(dir.tmp)) {
    dir.create(dir.tmp)
  }
}

# >>> read file >>>
gene.data <- read.csv(data.file)
gene.noms <- gene.data[1: 2]
gene.data <- gene.data[-2]
write.csv(gene.noms, 
          paste0(output.dir, "gene_noms.csv"), 
          row.names = F)

# >>> melt data >>>
melt.data <- melt(gene.data, 
                  id.vars = "Accession", 
                  variable.name = "Group", 
                  value.name = "Abundance") %>% 
  mutate(rep = gsub(".*_(\\d+)", "\\1", Group), 
         isotope = gsub("[YS]\\d+(\\w+)_.*", "\\1", Group), 
         time = gsub("[YS](\\d+).*", "\\1", Group), 
         group = gsub("([YS]).*", "\\1", Group))

df <- melt.data %>% 
  select(Accession, rep, time, isotope, group, Abundance) %>% 
  na.omit()
df$time <- as.numeric(df$time)
df <- unique(df)

df <- df %>%
  group_by(Accession, isotope) %>%
  filter(all(c("Y", "S") %in% group) && 
           n_distinct(time[group == "Y"]) >= 2 && n_distinct(time[group == "S"]) >= 2) %>%
  ungroup()

df.synthesis <- df %>% filter(isotope == "H")
df.degredation <- df %>% filter(isotope == "L")
df <- df.degredation

# >>> fit all proteins >>>
BPPARAM <- BiocParallel::SerialParam(progressbar = T)
fits <- MDfit(x = df$time, 
              y = df$Abundance, 
              id = df$Accession, 
              groupsNull = NULL, 
              groupsAlt = df$group, 
              BPPARAM = BPPARAM,
              returnModels = FALSE, 
              synthesis = F)
# >>> estimate degrees of freedom >>>
modelMetrics <- fits$metrics
fStats <- NPARCtest(modelMetrics, dfType = "empirical")

# >>> draw F stats distribution >>>
ggplot(filter(fStats, !is.na(pAdj))) +
  geom_density(aes(x = fStat), fill = "steelblue", alpha = 0.5) +
  geom_line(aes(x = fStat, y = df(fStat, df1 = df1, df2 = df2)), color = "darkred", size = 1.5) +
  theme_bw() +
  # Zoom in to small values to increase resolution for the proteins under H0:
  xlim(c(0, 10))
ggsave(paste0(output.dir, "F_stats_distribution.pdf"))

# >>> draw p distribution >>>
ggplot(filter(fStats, !is.na(pAdj))) +
  geom_histogram(aes(x = pVal, y = ..density..), fill = "steelblue", alpha = 0.5, boundary = 0, bins = 30) +
  geom_line(aes(x = pVal, y = dunif(pVal)), color = "darkred", size = 1.5) +
  theme_bw()
ggsave(paste0(output.dir, "p_distribution.pdf"))

save(fits, modelMetrics, fStats, file = "./synthesis.RData")
# save(fits, modelMetrics, fStats, file = "./degradation.RData")
save(df, file = "./data_df.RData")
