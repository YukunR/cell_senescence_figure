
# ---------------------- Normalization ------------------------------
#  >>> 读取文件 >>>
young.data <- read.delim("./data/you.txt")
sene.data <- read.delim("./data/sene.txt")

normalization <- function(data) {
  # >>> 预处理 >>> 
  data.noms <- data[1: 4]
  data <- data[-(1: 4)]
  # >>> 合并轻链和重链 >>> 
  data.colsum <- colSums(data, na.rm = T)
  data.colsum.target <- c()
  for (i in seq(1, length(data.colsum), by = 2)) {
    # 计算两列的和并赋值给新dataframe的对应列
    sum.tmp <- data.colsum[i] + data.colsum[i + 1]
    data.colsum.target[i] <- sum.tmp
    data.colsum.target[i + 1] <- sum.tmp
  }
  names(data.colsum.target) <- names(data.colsum)
  
  # >>> 列总和归一化 >>>
  norm.factors <- mean(data.colsum.target, na.rm = T) / data.colsum.target
  data.norm <- sweep(data, 2, norm.factors, FUN = "*")
  data.norm <- cbind(data.noms, data.norm)
  return(data.norm)
}

young.data.norm <- normalization(young.data)
sene.data.norm <- normalization(sene.data)

# ----------------------- 合并两个数据 ------------------------------
library(dplyr)
# >>> 质控 >>>
young.data.norm <- young.data.norm %>% 
  filter(Sum.PEP.Score >= 5 & X..Unique.Peptides > 1) %>% 
  select(-Sum.PEP.Score, -X..Unique.Peptides)
sene.data.norm <- sene.data.norm %>% 
  filter(Sum.PEP.Score >= 5 & X..Unique.Peptides > 1) %>% 
  select(-Sum.PEP.Score, -X..Unique.Peptides)
data.all <- merge(young.data.norm, sene.data.norm, 
                  by = c("Accession", "Description"))
write.csv(data.all, "./res/norm_data.csv", 
          row.names = F)