library(dplyr)
library(multiUS)

# ======================================================================
# Pulldown-specific imputation
# For proteins with high NA rates (those falling in Perseus zone),
# impute with a random value near 0 in log2 space (= ~1 in linear),
# representing undetected background signal.
# ======================================================================

#' Background-level imputation for pulldown data
#'
#' @description
#' Replaces NA values with random draws from N(mean=0, sd=0.2) in log2 space.
#' This corresponds to a linear intensity of ~1 (background/undetected level),
#' appropriate for co-IP data where absence of a protein reflects genuine
#' non-interaction rather than technical limitation.
#'
#' @param expression_matrix Matrix of log2-transformed values (proteins × samples)
#' @param mean_val Center of imputation distribution in log2 space (default 0 = linear 1)
#' @param sd_val SD of imputation distribution in log2 space (default 0.2)
#' @param seed Random seed for reproducibility (default 100)
#'
#' @return Matrix with NA values replaced
#'
#' @keywords internal
impute_background_level <- function(expression_matrix,
                                    mean_val = 0,
                                    sd_val   = 0.2,
                                    seed     = 100) {
  if (!is.matrix(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }

  n_missing <- sum(is.na(expression_matrix))
  if (n_missing == 0) {
    return(expression_matrix)
  }

  set.seed(seed)
  expression_matrix[is.na(expression_matrix)] <- rnorm(n_missing,
                                                        mean = mean_val,
                                                        sd   = sd_val)

  cat(sprintf("Background imputation: %d NAs filled with N(mean=%.1f, sd=%.2f) in log2 space\n",
              n_missing, mean_val, sd_val))
  return(expression_matrix)
}

#' Filter and impute pulldown data
#'
#' @description
#' Adapts the two-threshold filtering and imputation strategy for pulldown (co-IP) data.
#' Unlike standard proteomics, proteins absent from IgG control samples represent
#' true bait interactors; these are imputed with a background-level value near 0
#' in log2 space rather than the standard downshifted distribution.
#'
#' Two-threshold system (per group, NA ratio of each group evaluated):
#'   - NA ratio < threshold_low  → KNN imputation
#'   - threshold_low <= NA ratio < threshold_high → Background imputation (N(0, 0.2) log2)
#'   - NA ratio >= threshold_high → Protein discarded
#'
#' @param log2_data Log2-transformed matrix (proteins × samples, rownames = Accession)
#' @param sample_info Sample metadata data frame with Sample and Group columns
#' @param filter_threshold Two-value threshold vector c(low, high), default c(0.6, 0.9)
#' @param sample_col Column name for sample IDs in sample_info (default "Sample")
#' @param group_col Column name for group IDs in sample_info (default "Group")
#' @param knn_k K for KNN imputation (default 10)
#' @param background_mean Mean for background imputation in log2 space (default 0)
#' @param background_sd SD for background imputation in log2 space (default 0.2)
#' @param output_dir Output directory for saving intermediate results
#'
#' @return Imputed log2 matrix (proteins × samples)
#'
#' @export
filter_and_impute_pulldown <- function(log2_data,
                                       sample_info,
                                       filter_threshold = c(0.6, 0.9),
                                       sample_col       = "Sample",
                                       group_col        = "Group",
                                       knn_k            = 10,
                                       background_mean  = 0,
                                       background_sd    = 0.2,
                                       output_dir       = "./") {
  # Validate inputs
  if (length(filter_threshold) != 2) {
    stop("filter_threshold must be a 2-element vector c(low_threshold, high_threshold)")
  }
  threshold_low  <- filter_threshold[1]
  threshold_high <- filter_threshold[2]

  if (threshold_low > threshold_high) {
    stop(sprintf("filter_threshold: first value (%.2f) must be <= second value (%.2f)",
                 threshold_low, threshold_high))
  }

  required_cols <- c(sample_col, group_col)
  missing_cols  <- setdiff(required_cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in sample_info:", paste(missing_cols, collapse = ", ")))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Replace -Inf (from log2(0)) with NA
  log2_data[log2_data == -Inf] <- NA

  # --- Calculate per-protein NA rate in each group ---
  groups <- unique(sample_info[[group_col]])

  protein_max_na <- setNames(
    rep(0, nrow(log2_data)),
    rownames(log2_data)
  )

  for (group in groups) {
    samples_in_group <- sample_info[[sample_col]][sample_info[[group_col]] == group]
    valid_samples    <- intersect(samples_in_group, colnames(log2_data))
    if (length(valid_samples) == 0) next

    group_data <- log2_data[, valid_samples, drop = FALSE]
    na_ratio   <- rowSums(is.na(group_data)) / length(valid_samples)
    protein_max_na <- pmax(protein_max_na, na_ratio)
  }

  # --- Apply two-threshold filtering ---
  keep_proteins   <- names(protein_max_na)[protein_max_na < threshold_high]
  discarded_n     <- nrow(log2_data) - length(keep_proteins)

  cat(sprintf(
    "Filtering (two-threshold: %.2f / %.2f):\n  Original: %d proteins\n  Discarded (NA >= %.2f): %d\n  Retained: %d\n",
    threshold_low, threshold_high,
    nrow(log2_data), threshold_high,
    discarded_n, length(keep_proteins)
  ))

  filtered_data <- log2_data[keep_proteins, , drop = FALSE]

  # --- Impute per group ---
  imputed_data <- filtered_data

  for (group in groups) {
    samples_in_group <- sample_info[[sample_col]][sample_info[[group_col]] == group]
    valid_samples    <- intersect(samples_in_group, colnames(filtered_data))
    if (length(valid_samples) == 0) next

    group_data <- filtered_data[, valid_samples, drop = FALSE]
    n_samples  <- length(valid_samples)
    na_ratios  <- protein_max_na[rownames(group_data)]

    knn_proteins        <- rownames(group_data)[na_ratios < threshold_low]
    background_proteins <- rownames(group_data)[na_ratios >= threshold_low &
                                                na_ratios < threshold_high]

    cat(sprintf(
      "  Group '%s': %d proteins KNN (NA < %.2f), %d proteins background imputation (%.2f <= NA < %.2f)\n",
      group, length(knn_proteins), threshold_low,
      length(background_proteins), threshold_low, threshold_high
    ))

    # KNN imputation for low-NA proteins
    if (length(knn_proteins) > 0 && any(is.na(group_data[knn_proteins, , drop = FALSE]))) {
      knn_data    <- group_data[knn_proteins, , drop = FALSE]
      knn_imputed <- multiUS::KNNimp(knn_data, k = min(knn_k, n_samples - 1))
      knn_imputed <- as.data.frame(knn_imputed)
      colnames(knn_imputed) <- valid_samples
      imputed_data[knn_proteins, valid_samples] <- knn_imputed
    }

    # Background imputation for high-NA proteins
    if (length(background_proteins) > 0) {
      bg_data    <- as.matrix(group_data[background_proteins, , drop = FALSE])
      bg_imputed <- impute_background_level(bg_data,
                                            mean_val = background_mean,
                                            sd_val   = background_sd)
      imputed_data[background_proteins, valid_samples] <- bg_imputed
    }
  }

  cat("Pulldown imputation completed.\n")
  return(imputed_data)
}

# ======================================================================
# Pairwise FC and ΔFC calculation
# ======================================================================

#' Calculate pairwise log2 fold changes for pulldown data
#'
#' @description
#' Calculates per-replicate log2 fold change between experimental (pulldown) and
#' control (IgG) samples within each age group.
#'
#' For each protein and replicate i:
#'   log2FC_Young_i     = log2_YN_i - log2_YC_i
#'   log2FC_Senescent_i = log2_SN_i - log2_SC_i
#'
#' @param log2_data Log2 expression matrix (proteins × samples, rownames = Accession)
#' @param pairs_young List of c(expr_sample, ctrl_sample) for young replicates
#'   Default: list(c("YN1","YC1"), c("YN2","YC2"), c("YN3","YC3"))
#' @param pairs_senescent List of c(expr_sample, ctrl_sample) for senescent replicates
#'   Default: list(c("SN1","SC1"), c("SN2","SC2"), c("SN3","SC3"))
#'
#' @return Data frame with columns:
#'   Accession, log2FC_Y1, log2FC_Y2, log2FC_Y3, log2FC_S1, log2FC_S2, log2FC_S3
#'
#' @export
calculate_pairwise_fc <- function(log2_data,
                                  pairs_young = list(
                                    c("YN1","YC1"), c("YN2","YC2"), c("YN3","YC3")
                                  ),
                                  pairs_senescent = list(
                                    c("SN1","SC1"), c("SN2","SC2"), c("SN3","SC3")
                                  )) {
  result <- data.frame(Accession = rownames(log2_data),
                       stringsAsFactors = FALSE)

  # Young pairwise FCs
  for (i in seq_along(pairs_young)) {
    expr_s <- pairs_young[[i]][1]
    ctrl_s <- pairs_young[[i]][2]

    if (!expr_s %in% colnames(log2_data)) stop(paste("Column not found:", expr_s))
    if (!ctrl_s %in% colnames(log2_data)) stop(paste("Column not found:", ctrl_s))

    col_name <- paste0("log2FC_Y", i)
    result[[col_name]] <- as.numeric(log2_data[, expr_s]) - as.numeric(log2_data[, ctrl_s])
  }

  # Senescent pairwise FCs
  for (i in seq_along(pairs_senescent)) {
    expr_s <- pairs_senescent[[i]][1]
    ctrl_s <- pairs_senescent[[i]][2]

    if (!expr_s %in% colnames(log2_data)) stop(paste("Column not found:", expr_s))
    if (!ctrl_s %in% colnames(log2_data)) stop(paste("Column not found:", ctrl_s))

    col_name <- paste0("log2FC_S", i)
    result[[col_name]] <- as.numeric(log2_data[, expr_s]) - as.numeric(log2_data[, ctrl_s])
  }

  cat(sprintf("Pairwise FC calculated for %d proteins\n", nrow(result)))
  cat(sprintf("  Young replicates: %d pairs\n", length(pairs_young)))
  cat(sprintf("  Senescent replicates: %d pairs\n", length(pairs_senescent)))

  return(result)
}

#' Calculate ΔFC and perform statistical testing
#'
#' @description
#' For each protein, computes the delta log2 fold change between senescent and young
#' pulldown experiments (ΔFC = mean_log2FC_Senescent - mean_log2FC_Young), and tests
#' whether the two groups of FC values differ using an independent Welch t-test.
#'
#' Interpretation:
#'   ΔFC > 0: protein interaction with NPM1 is increased in senescent cells
#'   ΔFC < 0: protein interaction with NPM1 is decreased in senescent cells
#'
#' @param pairwise_fc_data Data frame from calculate_pairwise_fc()
#' @param protein_annotation Data frame with Accession, GeneName, Description columns
#' @param n_young Number of young replicates (default 3, expects log2FC_Y1..Yn columns)
#' @param n_senescent Number of senescent replicates (default 3, expects log2FC_S1..Sn columns)
#' @param min_valid_replicates Minimum non-NA FCs required per group for t-test (default 2)
#'
#' @return Data frame with columns:
#'   Accession, GeneName, Description, mean_log2FC_Young, mean_log2FC_Senescent,
#'   delta_log2FC, log2FC_Y1..Yn, log2FC_S1..Sn, p_value, regulation
#'
#' @export
calculate_delta_fc <- function(pairwise_fc_data,
                               protein_annotation,
                               n_young            = 3,
                               n_senescent        = 3,
                               min_valid_replicates = 2) {
  young_cols     <- paste0("log2FC_Y", seq_len(n_young))
  senescent_cols <- paste0("log2FC_S", seq_len(n_senescent))

  missing_y <- setdiff(young_cols,     colnames(pairwise_fc_data))
  missing_s <- setdiff(senescent_cols, colnames(pairwise_fc_data))
  if (length(missing_y) > 0) stop(paste("Missing young FC columns:", paste(missing_y, collapse = ", ")))
  if (length(missing_s) > 0) stop(paste("Missing senescent FC columns:", paste(missing_s, collapse = ", ")))

  n_proteins <- nrow(pairwise_fc_data)
  cat(sprintf("Calculating ΔFC and running t-tests for %d proteins...\n", n_proteins))

  # Extract FC matrices
  fc_young     <- as.matrix(pairwise_fc_data[, young_cols,     drop = FALSE])
  fc_senescent <- as.matrix(pairwise_fc_data[, senescent_cols, drop = FALSE])

  # Per-protein summary statistics
  mean_fc_young     <- rowMeans(fc_young,     na.rm = TRUE)
  mean_fc_senescent <- rowMeans(fc_senescent, na.rm = TRUE)
  delta_log2fc      <- mean_fc_senescent - mean_fc_young

  # Per-protein Welch t-test
  p_values <- vapply(seq_len(n_proteins), function(i) {
    y_vals <- fc_young[i, ]
    s_vals <- fc_senescent[i, ]

    y_valid <- y_vals[!is.na(y_vals)]
    s_valid <- s_vals[!is.na(s_vals)]

    if (length(y_valid) < min_valid_replicates || length(s_valid) < min_valid_replicates) {
      return(NA_real_)
    }

    tryCatch(
      t.test(s_valid, y_valid, var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))

  n_na_tests <- sum(is.na(p_values))
  if (n_na_tests > 0) {
    warning(sprintf("%d proteins had insufficient replicates for t-test (< %d valid values per group) and received NA p-value.",
                    n_na_tests, min_valid_replicates))
  }

  # Build result table
  result <- data.frame(
    Accession          = pairwise_fc_data$Accession,
    mean_log2FC_Young  = round(mean_fc_young, 4),
    mean_log2FC_Senescent = round(mean_fc_senescent, 4),
    delta_log2FC       = round(delta_log2fc, 4),
    p_value            = p_values,
    stringsAsFactors   = FALSE
  )

  # Add individual FC columns
  result <- cbind(result, pairwise_fc_data[, c(young_cols, senescent_cols)])

  # Merge annotation
  if (!is.null(protein_annotation)) {
    result <- merge(result, protein_annotation[, c("Accession","GeneName","Description")],
                    by = "Accession", all.x = TRUE)
    # Reorder columns
    front_cols <- c("Accession","GeneName","Description",
                    "mean_log2FC_Young","mean_log2FC_Senescent","delta_log2FC",
                    "p_value")
    fc_detail_cols <- c(young_cols, senescent_cols)
    result <- result[, c(front_cols, fc_detail_cols), drop = FALSE]
  }

  cat(sprintf("ΔFC calculation complete.\n"))
  cat(sprintf("  Proteins with valid t-test: %d\n", sum(!is.na(p_values))))

  return(result)
}

#' Add significance flags to ΔFC results
#'
#' @param delta_fc_results Data frame from calculate_delta_fc()
#' @param delta_fc_threshold Log2-scale ΔFC threshold (default 1.0, = 2-fold)
#' @param p_threshold P-value threshold (default 0.05)
#' @param require_positive_means Whether both mean log2FC values must be > 0 for significance (default TRUE)
#'
#' @return Data frame with added columns: significant, regulation
#'
#' @export
annotate_delta_fc <- function(delta_fc_results,
                              delta_fc_threshold = 1.0,
                              p_threshold        = 0.05,
                              require_positive_means = TRUE) {
  result <- delta_fc_results

  positive_means <- !require_positive_means |
    (result$mean_log2FC_Young > 0 & result$mean_log2FC_Senescent > 0)

  result$significant <- !is.na(result$p_value) &
    result$p_value < p_threshold &
    abs(result$delta_log2FC) >= delta_fc_threshold &
    positive_means

  result$regulation <- ifelse(
    result$significant & result$delta_log2FC >= delta_fc_threshold,
    "increased_in_senescent",
    ifelse(
      result$significant & result$delta_log2FC <= -delta_fc_threshold,
      "decreased_in_senescent",
      "no_change"
    )
  )

  n_up   <- sum(result$regulation == "increased_in_senescent", na.rm = TRUE)
  n_down <- sum(result$regulation == "decreased_in_senescent", na.rm = TRUE)
  n_ns   <- nrow(result) - n_up - n_down

  cat(sprintf("\n=== ΔFC Significance Summary ===\n"))
  cat(sprintf("  Thresholds: |Δlog2FC| >= %.2f AND p_value < %.3f%s\n",
              delta_fc_threshold, p_threshold,
              if (require_positive_means) " AND mean_log2FC_Young/Senescent > 0" else ""))
  cat(sprintf("  Increased in senescent: %d proteins\n", n_up))
  cat(sprintf("  Decreased in senescent: %d proteins\n", n_down))
  cat(sprintf("  No significant change:  %d proteins\n", n_ns))
  cat(sprintf("================================\n\n"))

  return(result)
}
