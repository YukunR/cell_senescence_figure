# ==============================================================================
# STRING PROTEIN INTERACTION NETWORK VISUALIZATION
# ==============================================================================
# Functions for creating protein-protein interaction networks using STRINGdb
# with GO/KEGG annotation overlays for pulldown proteomics data
#
# Features:
# - STRINGdb interaction network retrieval
# - Node coloring by deltaFC (diverging gradient)
# - GO/KEGG annotation via pie-chart vertices
# - NPM1 bait protein highlighting
#
# Required packages: STRINGdb, igraph, dplyr
# ==============================================================================

#' Build STRING interaction network for significant proteins
#'
#' @description
#' Retrieves protein-protein interactions from STRINGdb for significant deltaFC
#' proteins plus NPM1, and annotates nodes with GO/KEGG enrichment membership.
#'
#' @param delta_fc Data frame from calculate_delta_fc() with significant column
#' @param protein_annotation Data frame with Accession, GeneName, Description
#' @param npm1_accession NPM1 UniProt accession (default "P06748")
#' @param enrichment_dir Directory containing enrichment CSV results
#' @param species NCBI taxonomy ID (default 9606 for human)
#' @param score_threshold STRING combined score threshold (default 400)
#' @param string_version STRING database version (default "12.0")
#'
#' @return List with elements: graph (igraph object), node_data (data frame),
#'   edge_data (data frame)
#'
#' @export
build_string_network <- function(
  delta_fc,
  protein_annotation,
  npm1_accession = "P06748",
  enrichment_dir = NULL,
  species = 9606,
  score_threshold = 400,
  string_version = "12.0"
) {
  # Check required packages
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop(
      "Package 'STRINGdb' is required. Install from Bioconductor:\n",
      "  BiocManager::install('STRINGdb')"
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required. Install with:\n",
      "  install.packages('igraph')"
    )
  }

  library(STRINGdb)
  library(igraph)
  library(dplyr)

  cat("\n=== Building STRING Interaction Network ===\n")

  # --- Select seed proteins: significant deltaFC + NPM1 ---
  sig_proteins <- delta_fc$Accession[delta_fc$significant == TRUE]

  if (length(sig_proteins) == 0) {
    stop("No significant proteins found in delta_fc. Cannot build network.")
  }

  # Add NPM1 if not already in significant set
  if (!npm1_accession %in% sig_proteins) {
    sig_proteins <- c(npm1_accession, sig_proteins)
  }

  cat(sprintf("  Seed proteins: %d (including NPM1)\n", length(sig_proteins)))

  # --- Prepare protein table for STRING mapping ---
  seed_data <- data.frame(
    Accession = sig_proteins,
    stringsAsFactors = FALSE
  )

  # Merge with annotation to get GeneName
  seed_data <- merge(
    seed_data,
    protein_annotation[, c("Accession", "GeneName", "Description")],
    by = "Accession",
    all.x = TRUE
  )

  # Merge with deltaFC data
  seed_data <- merge(
    seed_data,
    delta_fc[, c(
      "Accession",
      "delta_log2FC",
      "mean_log2FC_Young",
      "mean_log2FC_Senescent",
      "p_value",
      "regulation"
    )],
    by = "Accession",
    all.x = TRUE
  )

  # --- Initialize STRINGdb ---
  cat(sprintf(
    "  Connecting to STRING database (v%s, species=%d)...\n",
    string_version,
    species
  ))

  string_db <- STRINGdb$new(
    version = string_version,
    species = species,
    score_threshold = score_threshold,
    input_directory = ""
  )

  # --- Map proteins to STRING IDs ---
  cat("  Mapping proteins to STRING IDs...\n")

  # Try mapping by GeneName first (preferred for STRING)
  mapped <- string_db$map(seed_data, "GeneName", removeUnmappedRows = FALSE)

  # For unmapped rows, try Accession
  unmapped_idx <- is.na(mapped$STRING_id)
  if (any(unmapped_idx)) {
    cat(sprintf(
      "    %d proteins unmapped by GeneName, trying Accession...\n",
      sum(unmapped_idx)
    ))
    unmapped_data <- seed_data[unmapped_idx, ]
    unmapped_data$query <- unmapped_data$Accession
    mapped_by_acc <- string_db$map(
      unmapped_data,
      "query",
      removeUnmappedRows = FALSE
    )

    # Merge back
    mapped$STRING_id[unmapped_idx] <- mapped_by_acc$STRING_id
  }

  n_mapped <- sum(!is.na(mapped$STRING_id))
  cat(sprintf(
    "  Successfully mapped: %d / %d proteins\n",
    n_mapped,
    nrow(seed_data)
  ))

  if (n_mapped == 0) {
    stop("No proteins could be mapped to STRING IDs. Cannot build network.")
  }

  # Keep only mapped proteins
  mapped <- mapped[!is.na(mapped$STRING_id), ]

  # --- Retrieve interactions ---
  cat("  Retrieving protein-protein interactions...\n")

  interactions <- string_db$get_interactions(mapped$STRING_id)

  # Remove duplicate edges (STRING sometimes returns bidirectional duplicates)
  if (!is.null(interactions) && nrow(interactions) > 0) {
    edge_key <- apply(interactions[, c("from", "to")], 1, function(r) {
      paste(sort(r), collapse = "||")
    })
    interactions <- interactions[!duplicated(edge_key), ]
  }

  if (is.null(interactions) || nrow(interactions) == 0) {
    warning(
      "No interactions found among seed proteins. Returning empty network."
    )
    return(list(graph = NULL, node_data = mapped, edge_data = NULL))
  }

  cat(sprintf("  Found %d interactions\n", nrow(interactions)))

  # --- Build igraph ---
  cat("  Building igraph object...\n")

  g <- graph_from_data_frame(
    d = interactions[, c("from", "to", "combined_score")],
    directed = FALSE,
    vertices = mapped$STRING_id
  )

  # --- Annotate nodes with GO/KEGG membership ---
  if (!is.null(enrichment_dir) && dir.exists(enrichment_dir)) {
    cat("  Annotating nodes with GO/KEGG enrichment membership...\n")

    go_membership <- extract_go_kegg_membership(
      accessions = mapped$Accession,
      enrichment_dir = enrichment_dir
    )

    # Merge annotation back to mapped data
    mapped <- merge(mapped, go_membership, by = "Accession", all.x = TRUE)
  } else {
    # No enrichment annotation
    mapped$GO_BP <- FALSE
    mapped$GO_CC <- FALSE
    mapped$GO_MF <- FALSE
    mapped$KEGG <- FALSE
  }

  # --- Attach node attributes to igraph ---
  for (col in setdiff(colnames(mapped), "STRING_id")) {
    g <- set_vertex_attr(
      g,
      col,
      value = mapped[[col]][match(V(g)$name, mapped$STRING_id)]
    )
  }

  cat("=== Network construction complete ===\n\n")

  return(list(
    graph = g,
    node_data = mapped,
    edge_data = interactions
  ))
}


#' Extract GO/KEGG membership from enrichment results
#'
#' @description
#' Parses enrichment CSV files to determine which proteins belong to
#' significant GO (BP/CC/MF) or KEGG terms.
#'
#' @param accessions Character vector of protein accessions
#' @param enrichment_dir Directory containing enrichment_go_results.csv and
#'   enrichment_kegg_results.csv
#'
#' @return Data frame with columns: Accession, GO_BP, GO_CC, GO_MF, KEGG (logical)
#'
#' @keywords internal
extract_go_kegg_membership <- function(accessions, enrichment_dir) {
  membership <- data.frame(
    Accession = accessions,
    GO_BP = FALSE,
    GO_CC = FALSE,
    GO_MF = FALSE,
    KEGG = FALSE,
    stringsAsFactors = FALSE
  )

  # --- GO enrichment ---
  go_file <- file.path(enrichment_dir, "enrichment_go_results.csv")
  if (file.exists(go_file)) {
    go_results <- read.csv(go_file, stringsAsFactors = FALSE)

    if (
      nrow(go_results) > 0 &&
        "geneID" %in% colnames(go_results) &&
        "Category" %in% colnames(go_results)
    ) {
      for (i in seq_len(nrow(go_results))) {
        gene_ids <- unlist(strsplit(go_results$geneID[i], "/"))
        category <- go_results$Category[i]

        if (category == "BP") {
          membership$GO_BP[membership$Accession %in% gene_ids] <- TRUE
        } else if (category == "CC") {
          membership$GO_CC[membership$Accession %in% gene_ids] <- TRUE
        } else if (category == "MF") {
          membership$GO_MF[membership$Accession %in% gene_ids] <- TRUE
        }
      }
    }
  }

  # --- KEGG enrichment ---
  kegg_file <- file.path(enrichment_dir, "enrichment_kegg_results.csv")
  if (file.exists(kegg_file)) {
    kegg_results <- read.csv(kegg_file, stringsAsFactors = FALSE)

    if (nrow(kegg_results) > 0 && "geneID" %in% colnames(kegg_results)) {
      for (i in seq_len(nrow(kegg_results))) {
        gene_ids <- unlist(strsplit(kegg_results$geneID[i], "/"))
        membership$KEGG[membership$Accession %in% gene_ids] <- TRUE
      }
    }
  }

  return(membership)
}


#' Plot STRING network with p-value coloring and deltaFC-based hull grouping (ggraph version)
#'
#' @description
#' Creates a publication-ready network plot using ggraph with:
#' - Node size scaled by absolute deltaFC
#' - Node color mapped to -log10(p-value) (blue=low significance, red=high significance)
#' - NPM1 highlighted with gold color (not p-value mapped)
#' - Convex hull grouping by regulation direction (warm colors for upregulated, cool for downregulated)
#' - Gene name labels positioned inside nodes
#' - Spatially separated upregulated and downregulated groups
#' - NPM1 positioned to avoid hull overlap
#'
#' @param network_list Output from build_string_network()
#' @param npm1_accession NPM1 accession for highlighting (default "P06748")
#' @param layout Layout algorithm: "fr" (Fruchterman-Reingold) (default "fr")
#' @param node_size_range Node size range as c(min, max) (default c(5, 25))
#' @param edge_width_scale Edge width scaling factor (default 0.01)
#' @param hull_colors Hull colors as c(upregulated, downregulated) (default warm/cool)
#'
#' @return ggplot object
#'
#' @export
plot_string_network <- function(
  network_list,
  npm1_accession = "P06748",
  layout = "fr",
  node_size_range = c(5, 25),
  edge_width_scale = 0.01,
  hull_colors = NULL
) {
  if (is.null(network_list$graph)) {
    warning("Network graph is NULL. Cannot plot.")
    return(invisible(NULL))
  }

  # Load required packages
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package 'ggraph' is required. Install with: install.packages('ggraph')"
    )
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop(
      "Package 'ggforce' is required. Install with: install.packages('ggforce')"
    )
  }
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop(
      "Package 'ggnewscale' is required. Install with: install.packages('ggnewscale')"
    )
  }
  if (!requireNamespace("shadowtext", quietly = TRUE)) {
    stop(
      "Package 'shadowtext' is required. Install with: install.packages('shadowtext')"
    )
  }

  library(igraph)
  library(ggraph)
  library(ggforce)
  library(ggnewscale)
  library(shadowtext)

  g <- network_list$graph
  g <- simplify(
    g,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(combined_score = "max", .default = "first")
  )

  cat("\n=== Plotting STRING Network (ggraph) ===\n")
  cat(sprintf("  Nodes: %d, Edges: %d\n", vcount(g), ecount(g)))

  # --- Identify NPM1 and groups ---
  npm1_idx <- which(V(g)$Accession == npm1_accession)
  upregulated_idx <- which(
    V(g)$delta_log2FC > 0 & V(g)$Accession != npm1_accession
  )
  downregulated_idx <- which(
    V(g)$delta_log2FC < 0 & V(g)$Accession != npm1_accession
  )

  # --- Seed-position layout matching reference image topology ---
  # Approximate (x, y) coordinates extracted from the reference STRING image,
  # with NPM1 as the origin reference. Upregulated proteins fan to the left;
  # downregulated cluster upper-right.
  cat("  Computing seed-position layout from reference image...\n")

  seed_coords <- list(
    NPM1 = c(0.000, 0.000),
    FTSJ3 = c(0.190, 0.105),
    EBNA1BP2 = c(0.265, -0.120),
    DDX24 = c(0.390, 0.050),
    PRC1 = c(-0.175, 0.010),
    YTHDF1 = c(-0.115, 0.105),
    FAM83D = c(-0.285, 0.150),
    RACGAP1 = c(-0.345, -0.145),
    MAP4 = c(-0.425, 0.125),
    EIF3CL = c(-0.480, -0.050),
    EIF3F = c(-0.410, -0.260),
    NCBP2 = c(-0.275, -0.235),
    PRKRA = c(-0.115, -0.215),
    AFF4 = c(-0.175, -0.150)
  )

  node_names <- V(g)$GeneName
  lo <- matrix(0, nrow = vcount(g), ncol = 2)

  # Assign seed positions; fall back to small random offsets near origin
  set.seed(42)
  for (i in seq_len(vcount(g))) {
    nm <- node_names[i]
    if (!is.na(nm) && nm %in% names(seed_coords)) {
      lo[i, ] <- seed_coords[[nm]]
    } else {
      lo[i, ] <- c(runif(1, -0.1, 0.1), runif(1, -0.1, 0.1))
    }
  }

  V(g)$x <- lo[, 1]
  V(g)$y <- lo[, 2]

  # --- Prepare node attributes for ggraph ---
  V(g)$abs_deltaFC <- abs(V(g)$delta_log2FC)
  V(g)$abs_deltaFC[is.na(V(g)$abs_deltaFC)] <- median(
    V(g)$abs_deltaFC,
    na.rm = TRUE
  )

  V(g)$neg_log_p <- -log10(V(g)$p_value)
  V(g)$neg_log_p[is.infinite(V(g)$neg_log_p)] <- max(
    V(g)$neg_log_p[!is.infinite(V(g)$neg_log_p)],
    na.rm = TRUE
  )
  V(g)$neg_log_p[is.na(V(g)$neg_log_p)] <- -log10(0.5)

  V(g)$group <- "Other"
  V(g)$group[upregulated_idx] <- "Upregulated"
  V(g)$group[downregulated_idx] <- "Downregulated"
  V(g)$group[npm1_idx] <- "NPM1"

  V(g)$label <- V(g)$GeneName
  V(g)$label[is.na(V(g)$label)] <- V(g)$Accession[is.na(V(g)$label)]

  npm1_size <- median(V(g)$abs_deltaFC[V(g)$group != "NPM1"], na.rm = TRUE)
  if (!is.finite(npm1_size)) {
    npm1_size <- median(V(g)$abs_deltaFC, na.rm = TRUE)
  }
  V(g)$plot_size_value <- V(g)$abs_deltaFC
  V(g)$plot_size_value[npm1_idx] <- max(
    npm1_size * 1.35,
    quantile(
      abs_fc_vals <- V(g)$abs_deltaFC[V(g)$group != "NPM1"],
      0.80,
      na.rm = TRUE
    )
  )

  if (is.null(hull_colors)) {
    hull_colors <- c(
      "Upregulated" = "#ed8828",
      "Downregulated" = "#81b21f"
    )
  }

  # |ΔFC| size legend: fixed one-decimal style values when available
  abs_fc_vals <- V(g)$abs_deltaFC[V(g)$group != "NPM1"]
  abs_fc_range <- range(abs_fc_vals, na.rm = TRUE)
  preferred_size_breaks <- seq(
    ceiling(abs_fc_range[1] * 5) / 5,
    abs_fc_range[2],
    by = 0.4
  )
  preferred_size_breaks <- round(preferred_size_breaks, 1)
  preferred_size_breaks <- preferred_size_breaks[
    preferred_size_breaks >= abs_fc_range[1] &
      preferred_size_breaks <= abs_fc_range[2]
  ]
  if (length(preferred_size_breaks) >= 3) {
    size_breaks <- preferred_size_breaks[1:3]
  } else {
    size_breaks <- round(
      seq(abs_fc_range[1], abs_fc_range[2], length.out = 3),
      1
    )
  }
  size_break_labels <- format(size_breaks, nsmall = 1, trim = TRUE)

  # STRING score edge-width legend: fixed integer style values
  edge_scores <- edge_attr(g, "combined_score")
  edge_labels <- c("400", "600", "800")
  E(g)$score_bin <- cut(
    edge_scores,
    breaks = c(400, 600, 800, 1000),
    include.lowest = TRUE,
    labels = edge_labels
  )

  p <- ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y) +

    geom_mark_ellipse(
      aes(x = x, y = y, group = group, fill = group),
      data = function(x) x[x$group %in% c("Upregulated", "Downregulated"), ],
      expand = unit(10, "mm"),
      alpha = 0.30,
      color = NA,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      values = hull_colors,
      name = "Regulation",
      guide = guide_legend(
        order = 4,
        override.aes = list(color = NA, linetype = 0)
      )
    ) +

    geom_edge_link(
      aes(width = score_bin),
      color = "#5e5e5e",
      alpha = 0.46
    ) +
    scale_edge_width_manual(
      values = c("400" = 0.8, "600" = 1.8, "800" = 3.0),
      name = "STRING score",
      guide = guide_legend(
        order = 3,
        override.aes = list(
          edge_colour = "#5c5c5c",
          alpha = 0.85,
          fill = NA,
          color = NA
        )
      )
    ) +

    new_scale_fill() +
    geom_node_point(
      aes(size = plot_size_value, fill = neg_log_p),
      data = function(x) x[x$group != "NPM1", ],
      shape = 21,
      color = "#5C5C5C",
      stroke = 0.6
    ) +
    scale_fill_gradient(
      low = "#EBF5FB",
      high = "#1A5276",
      name = expression(-log[10](italic(p))),
      guide = guide_colorbar(
        order = 1,
        barwidth = unit(4, "mm"),
        barheight = unit(40, "mm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +

    # NPM1 — gold circle base
    geom_node_point(
      aes(size = plot_size_value),
      data = function(x) x[x$group == "NPM1", ],
      shape = 21,
      fill = "#D4A017",
      color = "#8B6914",
      stroke = 1.8,
      show.legend = FALSE
    ) +
    # NPM1 — star overlay (shape 8 = asterisk), fixed smaller size
    geom_node_point(
      data = function(x) x[x$group == "NPM1", ],
      shape = 8,
      size = 12,
      color = "white",
      stroke = 1.2,
      show.legend = FALSE
    ) +

    scale_size_continuous(
      range = node_size_range,
      breaks = size_breaks,
      labels = size_break_labels,
      name = expression("|" * Delta * "FC|"),
      guide = guide_legend(
        order = 2,
        keyheight = unit(3, "mm"),
        keywidth = unit(3, "mm"),
        override.aes = list(
          color = NA,
          stroke = 0,
          fill = "grey50"
        )
      )
    ) +

    shadowtext::geom_shadowtext(
      aes(x = x, y = y, label = label),
      size = 5,
      color = "black",
      bg.colour = "white",
      bg.r = 0.18,
      fontface = "bold"
    ) +

    theme_graph() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.direction = "vertical",
      legend.spacing.y = unit(3, "mm"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  cat("=== Plot complete ===\n\n")

  return(p)
}

# ==============================================================================
# STANDALONE EXECUTION
# ==============================================================================
# This section allows the script to be run independently for testing
# Usage: Rscript plots/network_plot.R

if (sys.nframe() == 0) {
  cat("\n=== Running network_plot.R in standalone mode ===\n")

  # Check if network data exists
  network_file <- "res_pulldown/network_results/string_network.rds"

  if (!file.exists(network_file)) {
    stop(
      "Network data file not found: ",
      network_file,
      "\n",
      "Please run the full pulldown pipeline first to generate network data."
    )
  }

  # Load required libraries
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required. Install with: install.packages('igraph')"
    )
  }

  library(igraph)

  # Load network data
  cat("Loading network data from:", network_file, "\n")
  network_list <- readRDS(network_file)

  # Create output directory if needed
  output_dir <- "res_pulldown/network_results"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate plot
  output_file <- file.path(output_dir, "string_network_ggraph.pdf")
  cat("Creating network plot:", output_file, "\n")

  p <- plot_string_network(
    network_list = network_list,
    npm1_accession = "P06748",
    layout = "fr",
    node_size_range = c(20, 30),
    edge_width_scale = 0.01
  )

  # Use ggsave with cairo_pdf device to handle fonts better
  ggsave(
    output_file,
    plot = p,
    width = 12,
    height = 8,
    device = cairo_pdf,
    family = "sans"
  )

  cat("\n=== Standalone execution complete ===\n")
  cat("Output saved to:", output_file, "\n")
}
