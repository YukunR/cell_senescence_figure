library(dplyr)
library(NPARC)
library(stringr)
library(reshape2)
library(broom)
library(ggplot2)

# >>> config >>>
data.dir <- "./data/"
data.file <- c("rep1", "rep2")
output.dir <- "./res/NPARC//"
Temp0 <- 37.0
TempCount <- 10
drugConcentration <- "SEN"
drugName <- "Group"

PairedTtest.padj <- 0.5
NPARC.padj <- 0.05

# underdevelopment
filterDown <- F  # whether filtering lower stability after drug treatment

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
file.list <- list()
for (i in 1: length(data.file)) {
  cur.data <- read.delim(paste0(data.dir, data.file[i], ".txt"))
  file.list[[data.file[i]]] <- cur.data
}


# >>> quality control >>>
source("./qualityControl.R", local = T)
noms.list <- list()
dat.list <- list()
for (i in 1: length(data.file)) {
  gene.dat.tmp <- qualityControlDDA(file.list[[data.file[i]]])
  gene.dat.tmp <- separateNominalColDDA(gene.dat.tmp)
  gene.noms.tmp <- gene.dat.tmp$gene.noms
  gene.dat.tmp <- gene.dat.tmp$gene.dat
  noms.list[[data.file[i]]] <- gene.noms.tmp
  dat.list[[data.file[i]]] <- gene.dat.tmp
}
gene.noms <- unique(rbind(noms.list$rep1, noms.list$rep2))  # , noms.list$rep2
write.csv(gene.noms, 
          paste0(output.dir, "gene_noms.csv"), 
          row.names = F)

# >>> melt data >>>
melt.data <- melt(dat.list) %>% 
  mutate(Temp = as.double(gsub("temp(\\d+(\\.\\d+)?)_.*", "\\1", variable))) %>%
  mutate(Conc = str_extract(variable, "(?<=conc)\\w+")) %>% 
  mutate(L1 = str_extract(L1, "(?<=rep)[0-9]+"))
colnames(melt.data)[1: 4] <- c("Gene", "variable", "Abundance", "Replicate")
melt.data <- melt.data[-2]

# >>> calculate relative abundance >>>
melt.data <- melt.data %>% 
  group_by(Gene, Replicate, Conc) %>%
  mutate(RelAbundance = Abundance/Abundance[Temp == Temp0]) %>% 
  ungroup()

# >>> filter NA >>>
df <- melt.data
df <- df %>% 
  group_by(Gene) %>% 
  filter(!is.na(RelAbundance)) %>% 
  filter(length(RelAbundance) == 4*TempCount) %>%  # 4*TempCount
  ungroup()

# >>> filter lower stability after drug treatment >>>
if (filterDown == T) {
  df <- df %>% 
    arrange(Temp) %>%
    group_by(Gene) %>% 
    filter(mean(RelAbundance[which(Conc == drugConcentration)][7: 14]) >= 
             mean(RelAbundance[which(Conc == 0)][7: 14])) %>% 
    ungroup()
}


# >>> fit all proteins >>>
BPPARAM <- BiocParallel::SerialParam(progressbar = FALSE)
fits <- NPARCfit(x = df$Temp, 
                 y = df$RelAbundance, 
                 id = df$Gene, 
                 groupsNull = NULL, 
                 groupsAlt = df$Conc, 
                 BPPARAM = BPPARAM,
                 returnModels = FALSE)

# >>> estimate degrees of freedom >>>
modelMetrics <- fits$metrics 
fStats <- NPARCtest(modelMetrics, dfType = "empirical")
save(modelMetrics, fStats, file = "NPARCFit.RData")
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


# >>> select top proteins affected by drug >>>
topHits <- fStats %>% 
  filter(pVal <= NPARC.padj) %>%
  dplyr::select(id, fStat, pVal, pAdj) %>%
  merge(gene.noms, by.x = "id", by.y = "Accession") %>%
  arrange(-fStat) %>% 
  subset(select = c(id, GeneName, ProteinName, Description, pVal, pAdj, fStat))
colnames(topHits)[1] <- "Accession"
write.csv(topHits, 
          file = paste0(output.dir, "topHits.csv"), 
          row.names = F)
save(modelMetrics, fStats, file = "./NPARCFit.RData")
save(df, file = "./data_df.RData")