library(dplyr)
library(stringr)


# >>> quality control >>>
#' @param gene.dat gene data
#' @details
#' gene dat should and should only have columns: Accession, Description, Sum.PEP.Score, Unique.Peptides, ...
#' 
#' @param SPE.threshold Sum PEP Score threshold, default 5, p.s. > 5
#' @param UP.threshold Unique Peptides threshold, default 1, p.s. >= 1
qualityControlDDA <- function(gene.dat, SPE.threshold = 5, UP.threshold = 1) {
  gene.dat <- gene.dat %>% 
    # na.omit() %>%
    filter(Sum.PEP.Score > SPE.threshold & X..Unique.Peptides >= UP.threshold)
  return(gene.dat)
}

#' @param gene.dat gene data
#' @details
#' gene dat should and should only have columns: Accession, Description, Sum.PEP.Score, Unique.Peptides, ...
#' ... means Abundance data
#' 
#' @return a list of 2 dataframe, gene.nomes and gene.dat
separateNominalColDDA <- function(gene.dat) {
  gene.noms <- subset(gene.dat, select = c(Accession, Description))
  gene.dat <- subset(gene.dat, select = -c(Description, Sum.PEP.Score, X..Unique.Peptides))
  
  # >>> extract gene name from gene description >>>
  gene.noms$GeneName <- gsub(".*GN=([^ ]+).*", "\\1", gene.noms$Description)
  gene.noms$ProteinName <- gsub("(.*) OS=.*", "\\1", gene.noms$Description)
  
  return(list(gene.noms = gene.noms, gene.dat = gene.dat))
}