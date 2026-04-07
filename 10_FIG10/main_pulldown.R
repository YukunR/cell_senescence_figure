# ===========================================
# ==== NPM1 Pulldown Analysis - Main ====
# ===========================================
# DDA co-IP proteomics pipeline for NPM1 pulldown experiment
# Data from ProteomeDiscoverer (MS260425)
#
# Experiment:
#   - Young cells: YC1-3 (IgG control), YN1-3 (anti-NPM1 pulldown)
#   - Senescent cells: SC1-3 (IgG control), SN1-3 (anti-NPM1 pulldown)
#   - Paired design: YCi <-> YNi, SCi <-> SNi
#
# Goal: Identify proteins with changed NPM1 interaction in senescence (ΔFC)

# ===========================
# ========= Config ===========
# ===========================

# --- File paths ---
base_dir <- "./res_pulldown/"
protein_expr_file <- "./data/origin_data.txt" # PD Excel copied as tab-delimited
sample_info_file <- "./data/sample_info.txt" # auto-created if not present
go_background_file <- "./data/all_uniprot_go_background.csv"
kegg_background_file <- "./data/pathfromKegg_hsa.txt" # human

# --- QC filters (applied during import) ---
qc_pep_score_threshold <- 5 # Minimum Sum PEP Score
qc_unique_peptides <- 1 # Minimum # Unique Peptides

# --- Standard normalization ---
normalization_method <- "within_group" # "global" or "within_group"
use_common_proteins_for_norm <- TRUE
na_threshold <- c(0.6, 1) # Two-threshold: KNN / background imputation / discard

# --- NPM1 normalization ---
# Equalizes pulldown efficiency between young and senescent age groups
npm1_gene_id <- "NPM1"
npm1_accession <- "P06748" # Human NPM1 UniProt ID
young_expr_samples <- c("YN1", "YN2", "YN3")
senescent_expr_samples <- c("SN1", "SN2", "SN3")
senescent_all_samples <- c("SC1", "SC2", "SC3", "SN1", "SN2", "SN3") # All senescent samples to adjust

# --- Pairwise FC pairs (one-to-one within each age group) ---
pairs_young <- list(c("YN1", "YC1"), c("YN2", "YC2"), c("YN3", "YC3"))
pairs_senescent <- list(c("SN1", "SC1"), c("SN2", "SC2"), c("SN3", "SC3"))

# --- ΔFC thresholds ---
# ΔFC = mean_log2FC_Senescent - mean_log2FC_Young (independent Welch t-test)
delta_fc_threshold <- 1.0 # log2 scale (= 2-fold linear ΔFC)
p_threshold <- 0.05 # raw p-value threshold
require_positive_means <- TRUE # both mean_log2FC_Young and mean_log2FC_Senescent must be > 0

# Optional: custom colors for groups (uncomment to use)
# custom_colors <- c(YC = "#4ECDC4", YN = "#45B7D1", SC = "#FFA07A", SN = "#FF6B6B")

# ===========================
# ========= Examples =========
# ===========================

# # Run full pipeline
# result <- run_pulldown_analysis()
#
# # Check status
# check_pulldown_status()
#
# # Resume after interruption (same command)
# result <- run_pulldown_analysis()
#
# # Force rerun specific steps
# rerun_pulldown_step("normalization")
# rerun_pulldown_step(c("delta_fc_analysis", "enrichment"))

# ===========================
# ========= Run ==============
# ===========================

source("pulldown_steps.R")
result <- run_pulldown_analysis()
