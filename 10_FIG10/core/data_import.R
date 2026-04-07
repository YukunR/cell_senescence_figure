library(dplyr)
library(stringr)

# Sample mapping: file ID -> sample name
# Based on MS260425-readme.txt:
#   F1=YC1, F2=YC2, F3=YC3 (Young Control)
#   F4=YN1, F5=YN2, F6=YN3 (Young NPM1 pulldown)
#   F7=SC1, F8=SC2, F9=SC3 (Senescent Control)
#   F10=SN1, F11=SN2, F12=SN3 (Senescent NPM1 pulldown)
SAMPLE_MAP <- c(
  "F1"  = "YC1", "F2"  = "YC2", "F3"  = "YC3",
  "F4"  = "YN1", "F5"  = "YN2", "F6"  = "YN3",
  "F7"  = "SC1", "F8"  = "SC2", "F9"  = "SC3",
  "F10" = "SN1", "F11" = "SN2", "F12" = "SN3"
)

#' Import ProteomeDiscoverer output (copied as tab-delimited txt)
#'
#' @description
#' Reads the tab-delimited txt file copied from ProteomeDiscoverer Excel output,
#' applies QC filters, parses gene names from Description, and renames abundance columns.
#'
#' @param file_path Path to origin_data.txt (tab-delimited PD output)
#' @param qc_pep_score_threshold Minimum Sum PEP Score to keep a protein (default 5)
#' @param qc_unique_peptides Minimum number of unique peptides (default 1)
#' @param sample_map Named character vector mapping F-number codes to sample names.
#'   Default uses the MS260425 mapping (F1=YC1, ..., F12=SN3).
#'
#' @return Data frame with columns: Accession, GeneName, Description, YC1...SN3
#'
#' @export
import_pd_data <- function(file_path,
                           qc_pep_score_threshold = 5,
                           qc_unique_peptides = 1,
                           sample_map = SAMPLE_MAP) {
  cat("Reading ProteomeDiscoverer output from:", file_path, "\n")

  # Must use check.names = FALSE because column names contain #, spaces, and numbers
  data <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

  cat("Total proteins before QC filtering:", nrow(data), "\n")

  # --- QC Filtering ---
  pep_score_col <- "Sum PEP Score"
  unique_pep_col <- "# Unique Peptides"

  missing_qc_cols <- setdiff(c(pep_score_col, unique_pep_col), colnames(data))
  if (length(missing_qc_cols) > 0) {
    stop(paste("QC columns not found in data:", paste(missing_qc_cols, collapse = ", "),
               "\nPlease check that the file was copied correctly from ProteomeDiscoverer."))
  }

  data[[pep_score_col]] <- suppressWarnings(as.numeric(data[[pep_score_col]]))
  data[[unique_pep_col]] <- suppressWarnings(as.numeric(data[[unique_pep_col]]))

  before_filter <- nrow(data)
  data <- data[!is.na(data[[pep_score_col]]) &
               !is.na(data[[unique_pep_col]]) &
               data[[pep_score_col]] >= qc_pep_score_threshold &
               data[[unique_pep_col]] >= qc_unique_peptides, ]
  after_filter <- nrow(data)

  cat(sprintf(
    "QC filter applied (Sum PEP Score >= %g AND # Unique Peptides >= %g):\n",
    qc_pep_score_threshold, qc_unique_peptides
  ))
  cat(sprintf("  %d proteins kept, %d removed\n", after_filter, before_filter - after_filter))

  # --- Extract Gene Names from Description ---
  if (!"Description" %in% colnames(data)) {
    stop("Column 'Description' not found. Please verify the input file structure.")
  }
  if (!"Accession" %in% colnames(data)) {
    stop("Column 'Accession' not found. Please verify the input file structure.")
  }

  # Parse GN= field from Description
  # PD Description format: "PROTEIN_HUMAN Full name OS=Homo sapiens OX=9606 GN=NPM1 PE=1 SV=2"
  # Use [^;\\s]+ to stop at semicolon or whitespace (avoids capturing trailing ";")
  gn_extracted <- str_match(data$Description, "GN=([^;\\s]+)")[, 2]

  # Fall back to Accession for proteins without a gene name
  data$GeneName <- ifelse(is.na(gn_extracted) | nchar(trimws(gn_extracted)) == 0,
                          data$Accession,
                          gn_extracted)

  n_with_gn <- sum(!is.na(gn_extracted) & nchar(trimws(gn_extracted)) > 0)
  cat(sprintf(
    "Gene names parsed: %d with GN= field, %d using Accession as fallback\n",
    n_with_gn, nrow(data) - n_with_gn
  ))

  # --- Identify and Rename Raw Abundance Columns ---
  # Raw abundance columns match: "Abundance: F<number>:" (NOT "Abundances (Scaled)")
  all_cols <- colnames(data)
  # Match columns like "Abundance: F1: Sample" but not "Abundances (Scaled): F1: Sample"
  # Key: starts with "Abundance:" (no 's') and contains F<digit>
  abundance_col_idx <- grep("^Abundance: F\\d+", all_cols)

  if (length(abundance_col_idx) == 0) {
    stop(paste(
      "No raw abundance columns found (pattern 'Abundance: F<N>').",
      "Found columns:\n", paste(all_cols, collapse = "\n ")
    ))
  }

  abundance_col_names <- all_cols[abundance_col_idx]
  cat("Found", length(abundance_col_names), "raw abundance columns\n")

  # Extract F number from each column name
  f_codes <- str_extract(abundance_col_names, "F\\d+")

  # Map to sample names
  unmapped <- setdiff(f_codes, names(sample_map))
  if (length(unmapped) > 0) {
    warning(paste(
      "Some file codes have no sample mapping and will be skipped:",
      paste(unmapped, collapse = ", ")
    ))
  }

  mapped_idx <- f_codes %in% names(sample_map)
  abundance_col_names_mapped <- abundance_col_names[mapped_idx]
  f_codes_mapped <- f_codes[mapped_idx]
  sample_names_mapped <- sample_map[f_codes_mapped]

  # Build clean expression data frame
  expr_df <- as.data.frame(data[, abundance_col_names_mapped, drop = FALSE])
  colnames(expr_df) <- sample_names_mapped

  # Ensure columns are numeric
  expr_df <- as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(x))))

  # Sort columns in expected order (YC1-3, YN1-3, SC1-3, SN1-3)
  expected_order <- unname(sample_map)
  existing_order <- intersect(expected_order, colnames(expr_df))
  expr_df <- expr_df[, existing_order, drop = FALSE]

  # Build final output
  result <- data.frame(
    Accession   = data$Accession,
    GeneName    = data$GeneName,
    Description = data$Description,
    stringsAsFactors = FALSE
  )
  result <- cbind(result, expr_df)

  # Remove rows where Accession is NA or empty
  result <- result[!is.na(result$Accession) & nchar(trimws(result$Accession)) > 0, ]

  cat(sprintf("Import complete: %d proteins x %d samples\n", nrow(result), ncol(expr_df)))
  cat("Samples:", paste(existing_order, collapse = ", "), "\n")

  return(result)
}

#' Create sample_info.txt for the pulldown experiment
#'
#' @description
#' Generates the sample metadata file required by the normalization pipeline.
#' Creates Sample, Group, AgeGroup, Type, and Pair columns.
#'
#' @param output_file Path to write sample_info.txt (default: "./data/sample_info.txt")
#'
#' @return Invisible data frame of sample information
#'
#' @export
create_pulldown_sample_info <- function(output_file = "./data/sample_info.txt") {
  sample_info <- data.frame(
    Sample   = c("YC1","YC2","YC3","YN1","YN2","YN3",
                 "SC1","SC2","SC3","SN1","SN2","SN3"),
    Group    = c("YC","YC","YC","YN","YN","YN",
                 "SC","SC","SC","SN","SN","SN"),
    AgeGroup = c("Young","Young","Young","Young","Young","Young",
                 "Senescent","Senescent","Senescent","Senescent","Senescent","Senescent"),
    Type     = c("Control","Control","Control","Expr","Expr","Expr",
                 "Control","Control","Control","Expr","Expr","Expr"),
    Pair     = c(1,2,3,1,2,3,1,2,3,1,2,3),
    stringsAsFactors = FALSE
  )

  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(sample_info, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Sample info written to:", output_file, "\n")

  return(invisible(sample_info))
}
