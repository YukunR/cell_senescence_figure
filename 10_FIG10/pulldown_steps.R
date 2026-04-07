# ============================================================
# ==== NPM1 Pulldown Proteomics Analysis Pipeline ====
# ============================================================
# Adapted from the DIA analysis pipeline (analysis_steps.R)
# for DDA co-IP (pulldown) data from ProteomeDiscoverer.
#
# Key differences from the standard DIA pipeline:
#   1. Data import includes QC filtering (PEP score, unique peptides)
#   2. NPM1 pulldown efficiency normalization between age groups
#      is performed BEFORE imputation
#   3. Imputation uses background-level values (N(0,0.2) log2)
#      for high-NA proteins instead of Perseus downshift
#   4. Differential analysis uses pairwise FC (per replicate ctrl/expr)
#      and ΔFC = FC_Senescent - FC_Young with independent t-test

source("utils/checkpoint.R")

# ===========================
# ========= Config ==========
# ===========================

get_pulldown_config <- function() {
  config <- list()

  config$base_dir <- if (exists("base_dir", envir = .GlobalEnv)) {
    get("base_dir", envir = .GlobalEnv)
  } else {
    "./res_pulldown/"
  }
  config$protein_expr_file <- if (
    exists("protein_expr_file", envir = .GlobalEnv)
  ) {
    get("protein_expr_file", envir = .GlobalEnv)
  } else {
    "./data/origin_data.txt"
  }
  config$sample_info_file <- if (
    exists("sample_info_file", envir = .GlobalEnv)
  ) {
    get("sample_info_file", envir = .GlobalEnv)
  } else {
    "./data/sample_info.txt"
  }
  config$go_background_file <- if (
    exists("go_background_file", envir = .GlobalEnv)
  ) {
    get("go_background_file", envir = .GlobalEnv)
  } else {
    "./data/all_uniprot_go_background.csv"
  }
  config$kegg_background_file <- if (
    exists("kegg_background_file", envir = .GlobalEnv)
  ) {
    get("kegg_background_file", envir = .GlobalEnv)
  } else {
    "./data/pathfromKegg_hsa.txt"
  }
  config$custom_colors <- if (exists("custom_colors", envir = .GlobalEnv)) {
    get("custom_colors", envir = .GlobalEnv)
  } else {
    NULL
  }

  # QC thresholds (applied during data import)
  config$qc_pep_score_threshold <- if (
    exists("qc_pep_score_threshold", envir = .GlobalEnv)
  ) {
    get("qc_pep_score_threshold", envir = .GlobalEnv)
  } else {
    5
  }
  config$qc_unique_peptides <- if (
    exists("qc_unique_peptides", envir = .GlobalEnv)
  ) {
    get("qc_unique_peptides", envir = .GlobalEnv)
  } else {
    1
  }

  # Normalization / imputation
  config$normalization_method <- if (
    exists("normalization_method", envir = .GlobalEnv)
  ) {
    get("normalization_method", envir = .GlobalEnv)
  } else {
    "global"
  }
  config$use_common_proteins_for_norm <- if (
    exists("use_common_proteins_for_norm", envir = .GlobalEnv)
  ) {
    get("use_common_proteins_for_norm", envir = .GlobalEnv)
  } else {
    FALSE
  }
  config$na_threshold <- if (exists("na_threshold", envir = .GlobalEnv)) {
    get("na_threshold", envir = .GlobalEnv)
  } else {
    c(0.6, 0.9)
  }

  # NPM1 normalization
  config$npm1_gene_id <- if (exists("npm1_gene_id", envir = .GlobalEnv)) {
    get("npm1_gene_id", envir = .GlobalEnv)
  } else {
    "NPM1"
  }
  config$npm1_accession <- if (exists("npm1_accession", envir = .GlobalEnv)) {
    get("npm1_accession", envir = .GlobalEnv)
  } else {
    "P06748"
  }
  config$young_expr_samples <- if (
    exists("young_expr_samples", envir = .GlobalEnv)
  ) {
    get("young_expr_samples", envir = .GlobalEnv)
  } else {
    c("YN1", "YN2", "YN3")
  }
  config$senescent_expr_samples <- if (
    exists("senescent_expr_samples", envir = .GlobalEnv)
  ) {
    get("senescent_expr_samples", envir = .GlobalEnv)
  } else {
    c("SN1", "SN2", "SN3")
  }
  config$senescent_all_samples <- if (
    exists("senescent_all_samples", envir = .GlobalEnv)
  ) {
    get("senescent_all_samples", envir = .GlobalEnv)
  } else {
    c("SC1", "SC2", "SC3", "SN1", "SN2", "SN3")
  }

  # Pairwise FC pairs
  config$pairs_young <- if (exists("pairs_young", envir = .GlobalEnv)) {
    get("pairs_young", envir = .GlobalEnv)
  } else {
    list(c("YN1", "YC1"), c("YN2", "YC2"), c("YN3", "YC3"))
  }
  config$pairs_senescent <- if (exists("pairs_senescent", envir = .GlobalEnv)) {
    get("pairs_senescent", envir = .GlobalEnv)
  } else {
    list(c("SN1", "SC1"), c("SN2", "SC2"), c("SN3", "SC3"))
  }

  # ΔFC thresholds
  config$delta_fc_threshold <- if (
    exists("delta_fc_threshold", envir = .GlobalEnv)
  ) {
    get("delta_fc_threshold", envir = .GlobalEnv)
  } else {
    1.0
  }
  config$p_threshold <- if (exists("p_threshold", envir = .GlobalEnv)) {
    get("p_threshold", envir = .GlobalEnv)
  } else {
    0.05
  }
  config$require_positive_means <- if (
    exists("require_positive_means", envir = .GlobalEnv)
  ) {
    get("require_positive_means", envir = .GlobalEnv)
  } else {
    TRUE
  }

  return(config)
}

# =======================================================
# ========= Step 1: Import, Normalize, Impute ===========
# =======================================================

step_pulldown_normalization <- function(workspace, config = NULL) {
  if (is.null(config)) {
    config <- get_pulldown_config()
  }

  log_message(
    workspace,
    "Starting pulldown data import, normalization, and imputation"
  )

  source("core/data_import.R", local = TRUE)
  source("core/normalization.R", local = TRUE)
  source("core/npm1_normalization.R", local = TRUE)
  source("core/pairwise_fc.R", local = TRUE)

  norm_output_dir <- file.path(workspace$base_dir, "norm_results")
  if (!dir.exists(norm_output_dir)) {
    dir.create(norm_output_dir, recursive = TRUE)
  }

  # --- 0. Import & QC ---
  cat("\n--- [0] Data import and QC filtering ---\n")
  protein_data <- import_pd_data(
    file_path = config$protein_expr_file,
    qc_pep_score_threshold = config$qc_pep_score_threshold,
    qc_unique_peptides = config$qc_unique_peptides
  )

  # Create sample_info.txt if it doesn't exist
  if (!file.exists(config$sample_info_file)) {
    cat("sample_info.txt not found — creating automatically.\n")
    create_pulldown_sample_info(output_file = config$sample_info_file)
  }
  sample_info <- read.delim(config$sample_info_file, stringsAsFactors = FALSE)

  # Separate annotation from expression
  separated_data <- separate_protein_data(
    protein_data,
    handle_duplicates = "first",
    output_dir = norm_output_dir
  )
  protein_annotation <- separated_data$annotation_data
  expression_data <- separated_data$expression_data # linear, rownames = Accession

  calculate_na_percentage(expression_data, output_dir = norm_output_dir)

  # --- 1. Median normalization ---
  cat("\n--- [1] Median normalization ---\n")
  normalized_data <- normalize_by_median(
    expression_data,
    sample_info,
    normalization_method = config$normalization_method,
    use_common_proteins = config$use_common_proteins_for_norm
  )

  # --- 2. Log2 transform ---
  cat("\n--- [2] Log2 transformation ---\n")
  log2_data <- log2_transform(normalized_data)

  # --- 3. NPM1 pulldown efficiency normalization ---
  cat("\n--- [3] NPM1 pulldown efficiency normalization ---\n")
  npm1_result <- normalize_by_npm1(
    log2_data,
    protein_annotation,
    npm1_gene_id = config$npm1_gene_id,
    npm1_accession = config$npm1_accession,
    young_expr_samples = config$young_expr_samples,
    senescent_expr_samples = config$senescent_expr_samples,
    senescent_all_samples = config$senescent_all_samples
  )
  log2_npm1_normalized <- npm1_result$normalized_data

  # Save NPM1 normalization report
  npm1_report <- data.frame(
    metric = c(
      "delta_log2",
      "mean_log2_NPM1_Young",
      "mean_log2_NPM1_Senescent",
      paste0("NPM1_", names(npm1_result$npm1_values_young)),
      paste0("NPM1_", names(npm1_result$npm1_values_senescent))
    ),
    value = c(
      npm1_result$delta,
      npm1_result$mean_npm1_young,
      npm1_result$mean_npm1_senescent,
      npm1_result$npm1_values_young,
      npm1_result$npm1_values_senescent
    )
  )
  write.csv(
    npm1_report,
    file = file.path(norm_output_dir, "npm1_normalization_report.csv"),
    row.names = FALSE
  )

  # --- 4. Filter & Imputation (pulldown-specific) ---
  cat("\n--- [4] Filtering and imputation ---\n")
  imputed_data <- filter_and_impute_pulldown(
    log2_npm1_normalized,
    sample_info,
    filter_threshold = config$na_threshold,
    output_dir = norm_output_dir
  )

  # Visualize normalization workflow
  if (is.null(config$custom_colors)) {
    custom_colors <- generate_sample_colors(sample_info)$group_colors
  } else {
    custom_colors <- config$custom_colors
  }

  visualize_normalization_workflow(
    raw_data = expression_data,
    normalized_data = normalized_data,
    imputed_data = imputed_data,
    sample_info = sample_info,
    custom_colors = custom_colors,
    output_dir = norm_output_dir
  )

  save_processed_data(imputed_data, norm_output_dir)

  saveRDS(
    list(
      imputed_data = imputed_data,
      protein_annotation = protein_annotation,
      npm1_result = npm1_result,
      sample_info = sample_info,
      custom_colors = custom_colors,
      config = config
    ),
    file.path(workspace$base_dir, "normalization_results.rds")
  )

  log_message(workspace, "Normalization step completed")
  return("Normalization completed")
}

# ================================================
# ========= Step 2: PCA (Quality Check) ==========
# ================================================

step_pulldown_pca <- function(workspace) {
  log_message(workspace, "Starting PCA analysis")
  source("core/pca.R", local = TRUE)

  norm_results <- readRDS(file.path(
    workspace$base_dir,
    "normalization_results.rds"
  ))
  pca_output_dir <- file.path(workspace$base_dir, "pca_results")
  if (!dir.exists(pca_output_dir)) {
    dir.create(pca_output_dir, recursive = TRUE)
  }

  pca_input <- norm_results$imputed_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Accession")

  run_comprehensive_pca(
    expression_data = pca_input,
    sample_info = norm_results$sample_info,
    group_colors = norm_results$custom_colors,
    output_dir = pca_output_dir,
    plot_width = 10,
    plot_height = 12
  )

  log_message(workspace, "PCA analysis completed")
  return("PCA completed")
}

# =====================================================
# ========= Step 3: Pairwise FC and ΔFC ==============
# =====================================================

step_pairwise_fc_analysis <- function(workspace, config = NULL) {
  if (is.null(config)) {
    config <- get_pulldown_config()
  }
  log_message(workspace, "Starting pairwise FC and ΔFC analysis")

  source("core/pairwise_fc.R", local = TRUE)

  norm_results <- readRDS(file.path(
    workspace$base_dir,
    "normalization_results.rds"
  ))

  dfc_output_dir <- file.path(workspace$base_dir, "delta_fc_results")
  if (!dir.exists(dfc_output_dir)) {
    dir.create(dfc_output_dir, recursive = TRUE)
  }

  imputed_data <- norm_results$imputed_data
  protein_annotation <- norm_results$protein_annotation

  # --- Pairwise FC ---
  cat("\n--- Calculating pairwise FC ---\n")
  pairwise_fc <- calculate_pairwise_fc(
    imputed_data,
    pairs_young = config$pairs_young,
    pairs_senescent = config$pairs_senescent
  )

  write.csv(
    pairwise_fc,
    file = file.path(dfc_output_dir, "pairwise_fc.csv"),
    row.names = FALSE
  )

  # --- ΔFC + t-test ---
  cat("\n--- Calculating ΔFC and running t-tests ---\n")
  delta_fc_raw <- calculate_delta_fc(
    pairwise_fc,
    protein_annotation,
    n_young = length(config$pairs_young),
    n_senescent = length(config$pairs_senescent)
  )

  delta_fc <- annotate_delta_fc(
    delta_fc_raw,
    delta_fc_threshold = config$delta_fc_threshold,
    p_threshold = config$p_threshold,
    require_positive_means = config$require_positive_means
  )

  # Selection is driven by the finalized deltaFC annotation:
  # raw Welch p-value, |Δlog2FC| threshold, and positive mean FCs in both groups.
  write.csv(
    delta_fc,
    file = file.path(dfc_output_dir, "delta_fc_results.csv"),
    row.names = FALSE
  )

  # --- Volcano plot (ΔFC vs -log10 p_value) ---
  cat("\n--- Generating volcano plot ---\n")
  source("core/volcano_plot.R", local = TRUE)

  # Pass only statistical columns to avoid GeneName column conflict in the
  # left_join inside create_volcano_plot (both data and gene_annotations have GeneName).
  # gene_annotations provides the GeneName labels; annotation cols in main data would clash.
  stat_cols <- c(
    "Accession",
    "delta_log2FC",
    "p_value",
    "mean_log2FC_Young",
    "mean_log2FC_Senescent",
    "regulation"
  )
  volcano_data <- delta_fc[,
    intersect(stat_cols, colnames(delta_fc)),
    drop = FALSE
  ]
  volcano_data$log2_fold_change <- volcano_data$delta_log2FC

  volcano_results <- create_volcano_plot(
    volcano_data,
    fc_column = "log2_fold_change",
    p_column = "p_value",
    fc_threshold = config$delta_fc_threshold,
    p_threshold = config$p_threshold,
    gene_annotations = protein_annotation,
    annotation_counts = c(up = 10, down = 10),
    sort_method = "fc",
    regulation_column = "regulation" # use pre-annotated regulation (respects positive means filter)
  )

  # Relabel axes for clarity
  volcano_plot <- volcano_results$plot +
    ggplot2::labs(
      x = expression(paste(Delta, "log"[2], "FC (Senescent vs Young)")),
      y = expression(-log[10](italic(p))),
      title = "NPM1 Interaction Changes: Senescent vs Young"
    )

  ggsave(
    filename = file.path(dfc_output_dir, "volcano_delta_fc.pdf"),
    plot = volcano_plot,
    width = 7,
    height = 7
  )

  # --- Scatter plot: mean_log2FC_Young vs mean_log2FC_Senescent ---
  cat("\n--- Generating scatter plot ---\n")
  scatter_plot <- create_fc_scatter_plot(delta_fc, config)
  ggsave(
    filename = file.path(dfc_output_dir, "scatter_fc_young_vs_senescent.pdf"),
    plot = scatter_plot,
    width = 7,
    height = 7
  )

  saveRDS(
    list(
      pairwise_fc = pairwise_fc,
      delta_fc = delta_fc,
      protein_annotation = protein_annotation,
      config = config
    ),
    file.path(workspace$base_dir, "delta_fc_results.rds")
  )

  log_message(workspace, "Pairwise FC and ΔFC analysis completed")
  return("ΔFC analysis completed")
}

#' Create scatter plot of mean log2FC Young vs Senescent
#' @keywords internal
create_fc_scatter_plot <- function(delta_fc, config) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package 'ggrepel' is required for the scatter plot. Install with: install.packages('ggrepel')"
    )
  }
  library(ggplot2)

  sig_data <- delta_fc[delta_fc$significant == TRUE, ]
  nonsig_data <- delta_fc[delta_fc$significant != TRUE, ]

  p <- ggplot(delta_fc, aes(x = mean_log2FC_Young, y = mean_log2FC_Senescent)) +
    geom_point(data = nonsig_data, color = "grey70", alpha = 0.5, size = 1) +
    geom_point(
      data = sig_data[sig_data$regulation == "increased_in_senescent", ],
      color = "#E74C3C",
      alpha = 0.8,
      size = 1.5
    ) +
    geom_point(
      data = sig_data[sig_data$regulation == "decreased_in_senescent", ],
      color = "#3498DB",
      alpha = 0.8,
      size = 1.5
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    geom_abline(
      slope = 1,
      intercept = config$delta_fc_threshold,
      linetype = "dotted",
      color = "grey40"
    ) +
    geom_abline(
      slope = 1,
      intercept = -config$delta_fc_threshold,
      linetype = "dotted",
      color = "grey40"
    ) +
    labs(
      x = expression("mean log"[2] * "FC (Young NPM1/IgG)"),
      y = expression("mean log"[2] * "FC (Senescent NPM1/IgG)"),
      title = "Pairwise FC comparison: Young vs Senescent"
    ) +
    theme_classic(base_size = 12)

  if (nrow(sig_data) > 0 && "GeneName" %in% colnames(sig_data)) {
    top_genes <- sig_data[
      order(abs(sig_data$delta_log2FC), decreasing = TRUE),
    ]
    top_genes <- head(top_genes, 20)
    p <- p +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = GeneName),
        size = 2.5,
        max.overlaps = 20
      )
  }

  return(p)
}

# =====================================================
# ========= Step 4: Enrichment Analysis ==============
# =====================================================

step_pulldown_enrichment <- function(workspace, config = NULL) {
  if (is.null(config)) {
    config <- get_pulldown_config()
  }
  log_message(workspace, "Starting enrichment analysis on ΔFC results")

  source("core/enrichment_analysis.R", local = TRUE)

  dfc_results <- readRDS(file.path(workspace$base_dir, "delta_fc_results.rds"))
  delta_fc <- dfc_results$delta_fc

  enrich_base_dir <- file.path(workspace$base_dir, "enrichment_results")

  go_background <- read.csv(config$go_background_file)
  kegg_background <- read.delim(config$kegg_background_file)

  regulation_groups <- list(
    increased = delta_fc[
      delta_fc$regulation == "increased_in_senescent",
      "Accession"
    ],
    decreased = delta_fc[
      delta_fc$regulation == "decreased_in_senescent",
      "Accession"
    ],
    all_sig = delta_fc[delta_fc$significant == TRUE, "Accession"]
  )

  for (reg_type in names(regulation_groups)) {
    proteins <- regulation_groups[[reg_type]]
    if (length(proteins) == 0) {
      cat(sprintf(
        "No proteins in group '%s', skipping enrichment.\n",
        reg_type
      ))
      next
    }

    enrich_dir <- file.path(enrich_base_dir, reg_type)
    if (!dir.exists(enrich_dir)) {
      dir.create(enrich_dir, recursive = TRUE)
    }

    cat(sprintf(
      "\n--- Enrichment for '%s' (%d proteins) ---\n",
      reg_type,
      length(proteins)
    ))
    tryCatch(
      {
        run_combined_analysis(
          gene_list = proteins,
          go_background = go_background,
          kegg_background = kegg_background,
          output_dir = enrich_dir
        )
      },
      error = function(e) {
        warning(paste(
          "Enrichment failed for group",
          reg_type,
          ":",
          conditionMessage(e)
        ))
      }
    )
  }

  log_message(workspace, "Enrichment analysis completed")
  return("Enrichment completed")
}

# =====================================================
# ========= Step 5: STRING Network Analysis ===========
# =====================================================

step_pulldown_network <- function(workspace, config = NULL) {
  if (is.null(config)) {
    config <- get_pulldown_config()
  }
  log_message(workspace, "Starting STRING network analysis")

  source("plots/network_plot.R", local = TRUE)

  dfc_results <- readRDS(file.path(workspace$base_dir, "delta_fc_results.rds"))
  delta_fc <- dfc_results$delta_fc
  protein_annotation <- dfc_results$protein_annotation

  network_output_dir <- file.path(workspace$base_dir, "network_results")
  if (!dir.exists(network_output_dir)) {
    dir.create(network_output_dir, recursive = TRUE)
  }

  # Use enrichment results from all_sig group
  enrichment_dir <- file.path(
    workspace$base_dir,
    "enrichment_results",
    "all_sig"
  )

  cat("\n--- Building STRING interaction network ---\n")
  tryCatch(
    {
      network_list <- build_string_network(
        delta_fc = delta_fc,
        protein_annotation = protein_annotation,
        npm1_accession = config$npm1_accession,
        enrichment_dir = enrichment_dir,
        species = 9606,
        score_threshold = 400
      )

      # Save network data
      saveRDS(network_list, file.path(network_output_dir, "string_network.rds"))

      if (!is.null(network_list$node_data)) {
        write.csv(
          network_list$node_data,
          file.path(network_output_dir, "network_nodes.csv"),
          row.names = FALSE
        )
      }

      if (!is.null(network_list$edge_data)) {
        write.csv(
          network_list$edge_data,
          file.path(network_output_dir, "network_edges.csv"),
          row.names = FALSE
        )
      }

      # Plot network (ggraph version with improved NPM1 separation)
      if (!is.null(network_list$graph)) {
        cat("\n--- Generating network plot ---\n")
        p <- plot_string_network(
          network_list,
          npm1_accession = config$npm1_accession,
          layout = "fr",
          node_size_range = c(5, 25),
          edge_width_scale = 0.01
        )

        # Use ggsave with cairo_pdf to handle fonts better
        output_path <- file.path(network_output_dir, "string_network.pdf")
        ggsave(
          output_path,
          plot = p,
          width = 16,
          height = 12,
          device = cairo_pdf,
          family = "sans"
        )

        cat("Network plot saved to:", output_path, "\n")
      }
    },
    error = function(e) {
      warning(paste("STRING network analysis failed:", conditionMessage(e)))
      cat("Note: STRINGdb package may not be installed. Install with:\n")
      cat("  BiocManager::install('STRINGdb')\n")
    }
  )

  log_message(workspace, "STRING network analysis completed")
  return("Network analysis completed")
}

# ============================================
# ========= Main Pipeline Function ===========
# ============================================

run_pulldown_analysis <- function(
  project_name = "npm1_pulldown",
  force_rerun_list = NULL
) {
  config <- get_pulldown_config()
  workspace <- create_workspace(project_name, config$base_dir)

  log_message(workspace, "Starting NPM1 pulldown analysis pipeline")
  log_message(workspace, paste("Expression file:", config$protein_expr_file))
  log_message(workspace, paste("Output directory:", config$base_dir))

  checkpoint <- load_checkpoint(workspace)

  if (length(checkpoint$completed_steps) > 0) {
    cat(
      "Resuming from checkpoint. Completed steps:",
      paste(checkpoint$completed_steps, collapse = ", "),
      "\n"
    )
  } else {
    cat("Starting fresh analysis.\n")
  }

  if (!is.null(force_rerun_list)) {
    force_rerun_steps(force_rerun_list, project_name, config$base_dir)
  }

  # Step 1: Import + Normalize + NPM1 norm + Impute
  cat("\n=== Step 1: Import, Normalization, NPM1 Correction, Imputation ===\n")
  execute_step(
    workspace = workspace,
    step_name = "normalization",
    step_function = step_pulldown_normalization,
    output_files = c("normalization_results.rds"),
    config_to_track = list(
      normalization_method = config$normalization_method,
      na_threshold = config$na_threshold,
      npm1_gene_id = config$npm1_gene_id,
      qc_pep_score_threshold = config$qc_pep_score_threshold,
      qc_unique_peptides = config$qc_unique_peptides
    ),
    cleanup_patterns = c("norm_results/.*"),
    config = config
  )

  # Step 2: PCA
  cat("\n=== Step 2: PCA (Quality Check) ===\n")
  execute_step(
    workspace = workspace,
    step_name = "pca",
    step_function = step_pulldown_pca,
    output_files = c("pca_results/pca_biplot_PC1_PC2.pdf"),
    dependencies = "normalization",
    cleanup_patterns = c("pca_results/.*")
  )

  # Step 3: Pairwise FC and ΔFC analysis
  cat("\n=== Step 3: Pairwise FC and ΔFC Analysis ===\n")
  execute_step(
    workspace = workspace,
    step_name = "delta_fc_analysis",
    step_function = step_pairwise_fc_analysis,
    output_files = c(
      "delta_fc_results.rds",
      "delta_fc_results/delta_fc_results.csv",
      "delta_fc_results/volcano_delta_fc.pdf"
    ),
    dependencies = "normalization",
    config_to_track = list(
      pairs_young = config$pairs_young,
      pairs_senescent = config$pairs_senescent,
      delta_fc_threshold = config$delta_fc_threshold,
      p_threshold = config$p_threshold,
      require_positive_means = config$require_positive_means
    ),
    cleanup_patterns = c("delta_fc_results/.*"),
    config = config
  )

  # Step 4: Enrichment
  cat("\n=== Step 4: Enrichment Analysis ===\n")
  execute_step(
    workspace = workspace,
    step_name = "enrichment",
    step_function = step_pulldown_enrichment,
    output_files = NULL,
    dependencies = "delta_fc_analysis",
    cleanup_patterns = c("enrichment_results/.*"),
    config = config
  )

  # Step 5: STRING Network
  cat("\n=== Step 5: STRING Network Analysis ===\n")
  execute_step(
    workspace = workspace,
    step_name = "network",
    step_function = step_pulldown_network,
    output_files = NULL,
    dependencies = "enrichment",
    cleanup_patterns = c("network_results/.*"),
    config = config
  )

  log_message(workspace, "NPM1 pulldown analysis pipeline completed!")
  cat("\n=== Pipeline Completed! Results in:", config$base_dir, "===\n")

  return(list(
    workspace = workspace,
    status = "completed",
    results_dir = config$base_dir,
    config = config
  ))
}

# =====================================
# ========= Utility Functions =========
# =====================================

rerun_pulldown_step <- function(step_names, project_name = "npm1_pulldown") {
  run_pulldown_analysis(project_name, force_rerun_list = step_names)
}

check_pulldown_status <- function(project_name = "npm1_pulldown") {
  config <- get_pulldown_config()
  check_project_status(project_name, config$base_dir)
}
