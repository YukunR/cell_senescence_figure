library(dplyr)

#' Normalize pulldown data by NPM1 bait protein abundance
#'
#' @description
#' Adjusts all senescent samples (both SC and SN) so that the mean NPM1 abundance
#' in the senescent experimental group (SN) matches that in the young experimental
#' group (YN). This corrects for differences in pulldown efficiency between age groups,
#' allowing downstream ΔFC comparison to reflect true changes in protein interactions
#' rather than differences in how much NPM1 was captured.
#'
#' The normalization is performed in log2 space:
#'   delta = mean(log2_NPM1_YN_samples) - mean(log2_NPM1_SN_samples)
#'   All senescent samples (SC + SN) receive: value + delta
#'
#' @param log2_data Log2-transformed expression matrix (rownames = Accession or GeneName)
#' @param protein_annotation Data frame with Accession and GeneName columns
#' @param npm1_gene_id Gene name to locate NPM1 row (default "NPM1")
#' @param npm1_accession UniProt accession fallback for NPM1 (default "P06748" for human NPM1)
#' @param young_expr_samples Sample names for young experimental (NPM1 pulldown) group
#'   (default c("YN1","YN2","YN3"))
#' @param senescent_expr_samples Sample names for senescent experimental group
#'   (default c("SN1","SN2","SN3"))
#' @param senescent_all_samples All senescent samples to be adjusted (expr + ctrl)
#'   (default c("SC1","SC2","SC3","SN1","SN2","SN3"))
#'
#' @return List with elements:
#'   \item{normalized_data}{Adjusted log2 expression matrix}
#'   \item{delta}{The log2-scale adjustment applied to senescent samples}
#'   \item{npm1_values_young}{NPM1 log2 values in young experimental samples}
#'   \item{npm1_values_senescent}{NPM1 log2 values in senescent experimental samples}
#'   \item{npm1_row}{Row name identifying NPM1 in the matrix}
#'
#' @export
normalize_by_npm1 <- function(log2_data,
                              protein_annotation,
                              npm1_gene_id       = "NPM1",
                              npm1_accession     = "P06748",
                              young_expr_samples = c("YN1","YN2","YN3"),
                              senescent_expr_samples = c("SN1","SN2","SN3"),
                              senescent_all_samples  = c("SC1","SC2","SC3",
                                                         "SN1","SN2","SN3")) {

  # --- Locate NPM1 row ---
  npm1_row <- find_npm1_row(log2_data, protein_annotation, npm1_gene_id, npm1_accession)

  cat("NPM1 identified as row:", npm1_row, "\n")

  # --- Validate sample columns exist ---
  missing_y <- setdiff(young_expr_samples,     colnames(log2_data))
  missing_s <- setdiff(senescent_expr_samples, colnames(log2_data))
  missing_a <- setdiff(senescent_all_samples,  colnames(log2_data))

  if (length(missing_y) > 0) {
    stop(paste("Young experimental samples not found in data:", paste(missing_y, collapse = ", ")))
  }
  if (length(missing_s) > 0) {
    stop(paste("Senescent experimental samples not found in data:", paste(missing_s, collapse = ", ")))
  }
  # Only warn for the senescent adjustment set (SC may genuinely lack NPM1)
  if (length(missing_a) > 0) {
    warning(paste("Some senescent samples not found:", paste(missing_a, collapse = ", ")))
    senescent_all_samples <- intersect(senescent_all_samples, colnames(log2_data))
  }

  # --- Extract NPM1 values ---
  npm1_vals_young <- as.numeric(log2_data[npm1_row, young_expr_samples])
  npm1_vals_senescent <- as.numeric(log2_data[npm1_row, senescent_expr_samples])

  # Check for NAs in NPM1 values
  if (any(is.na(npm1_vals_young))) {
    na_samples <- young_expr_samples[is.na(npm1_vals_young)]
    stop(paste(
      "NPM1 is NA (missing) in young experimental sample(s):",
      paste(na_samples, collapse = ", "),
      "\nNPM1 must be detected in all experimental samples for normalization.",
      "\nConsider lowering na_threshold or checking the input data."
    ))
  }
  if (any(is.na(npm1_vals_senescent))) {
    na_samples <- senescent_expr_samples[is.na(npm1_vals_senescent)]
    stop(paste(
      "NPM1 is NA (missing) in senescent experimental sample(s):",
      paste(na_samples, collapse = ", "),
      "\nNPM1 must be detected in all experimental samples for normalization."
    ))
  }

  # --- Compute scaling delta ---
  mean_npm1_young     <- mean(npm1_vals_young,     na.rm = TRUE)
  mean_npm1_senescent <- mean(npm1_vals_senescent, na.rm = TRUE)
  delta <- mean_npm1_young - mean_npm1_senescent

  # --- Report ---
  cat("\n=== NPM1 Pulldown Efficiency Normalization ===\n")
  cat(sprintf("  Young NPM1 (log2): %s  →  mean = %.4f\n",
              paste(round(npm1_vals_young, 3), collapse = ", "), mean_npm1_young))
  cat(sprintf("  Senescent NPM1 (log2): %s  →  mean = %.4f\n",
              paste(round(npm1_vals_senescent, 3), collapse = ", "), mean_npm1_senescent))
  cat(sprintf("  Adjustment delta (log2): %.4f  (= 2^%.4f = %.3f-fold in linear)\n",
              delta, delta, 2^delta))
  cat(sprintf("  Applying +%.4f to %d senescent samples: %s\n",
              delta, length(senescent_all_samples),
              paste(senescent_all_samples, collapse = ", ")))
  cat("==============================================\n\n")

  if (abs(delta) > 3) {
    warning(sprintf(
      "NPM1 normalization delta is large (%.2f, = %.1f-fold). ",
      delta, 2^abs(delta)
    ))
  }

  # --- Apply adjustment ---
  normalized_data <- log2_data
  normalized_data[, senescent_all_samples] <- normalized_data[, senescent_all_samples] + delta

  return(list(
    normalized_data        = normalized_data,
    delta                  = delta,
    npm1_values_young      = setNames(npm1_vals_young, young_expr_samples),
    npm1_values_senescent  = setNames(npm1_vals_senescent, senescent_expr_samples),
    npm1_row               = npm1_row,
    mean_npm1_young        = mean_npm1_young,
    mean_npm1_senescent    = mean_npm1_senescent
  ))
}

#' Locate the NPM1 row in a log2 expression matrix
#'
#' @param log2_data Log2 expression matrix (rownames = Accession)
#' @param protein_annotation Data frame with Accession and GeneName columns
#' @param npm1_gene_id Gene name to search (default "NPM1")
#' @param npm1_accession UniProt accession fallback (default "P06748")
#'
#' @return Row name (character) identifying NPM1
#'
#' @keywords internal
find_npm1_row <- function(log2_data, protein_annotation,
                          npm1_gene_id = "NPM1",
                          npm1_accession = "P06748") {
  row_names <- rownames(log2_data)

  # Strategy 1: match by GeneName in annotation
  if (!is.null(protein_annotation) && "GeneName" %in% colnames(protein_annotation)) {
    gene_match <- protein_annotation$Accession[
      toupper(protein_annotation$GeneName) == toupper(npm1_gene_id)
    ]
    gene_match <- intersect(gene_match, row_names)
    if (length(gene_match) == 1) {
      return(gene_match)
    }
    if (length(gene_match) > 1) {
      warning(paste("Multiple rows match gene name", npm1_gene_id,
                    "- using first:", gene_match[1]))
      return(gene_match[1])
    }
  }

  # Strategy 2: rownames directly match gene name (in case rownames are gene names)
  if (toupper(npm1_gene_id) %in% toupper(row_names)) {
    idx <- which(toupper(row_names) == toupper(npm1_gene_id))
    return(row_names[idx[1]])
  }

  # Strategy 3: match by UniProt accession
  if (npm1_accession %in% row_names) {
    return(npm1_accession)
  }

  # Strategy 4: partial accession match (e.g., "P06748-1" isoform)
  partial_match <- row_names[startsWith(row_names, npm1_accession)]
  if (length(partial_match) >= 1) {
    warning(paste("Using partial accession match for NPM1:", partial_match[1]))
    return(partial_match[1])
  }

  stop(paste(
    "NPM1 not found in expression data.",
    "\nSearched for gene name:", npm1_gene_id,
    "\nSearched for accession:", npm1_accession,
    "\nPlease check that NPM1 passed the QC and NA filters.",
    "\nAvailable row names (first 10):", paste(head(row_names, 10), collapse = ", ")
  ))
}
