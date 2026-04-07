# ==============================================================================
# 3D INTERACTIVE PLOT: Young vs Senescent FC + p-value (Standalone)
# ==============================================================================
# Axes:   X = mean_log2FC_Young, Y = mean_log2FC_Senescent, Z = -log10(p_value)
# Color:  delta_log2FC (blue < 0, white = 0, red > 0)
#         Grey (low opacity) for: p_value >= 0.05 OR either mean < 0
#
# Usage:
#   setwd("d:/ResearchProject/NPM1Pulldown_TianxiangWang")
#   source("plots/3d_delta_fc.R")
#
# Output:
#   res_pulldown/delta_fc_results/delta_fc_3d.html  (always)
#   res_pulldown/delta_fc_results/delta_fc_3d.pdf   (requires: RSelenium, rsvg)
#   res_pulldown/delta_fc_results/delta_fc_3d.png   (requires: RSelenium, rsvg)
# ==============================================================================

if (!requireNamespace("plotly", quietly = TRUE)) {
  stop("Install plotly: install.packages('plotly')")
}
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  stop("Install htmlwidgets: install.packages('htmlwidgets')")
}

library(plotly)
library(htmlwidgets)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
p_threshold <- 0.05 # must match main_pulldown.R
delta_fc_threshold <- 1.0 # must match main_pulldown.R
xy_zoom <- 5 # |mean_log2FC| > xy_zoom shown as diamonds at boundary

# ==============================================================================
# LOAD DATA
# ==============================================================================
data_path <- "res_pulldown/delta_fc_results/delta_fc_results.csv"
output_dir <- "res_pulldown/delta_fc_results"

if (!file.exists(data_path)) {
  stop("Data not found: ", data_path)
}

df <- read.csv(data_path, stringsAsFactors = FALSE)
df$neg_log10_p <- -log10(df$p_value)

cat("Loaded", nrow(df), "proteins\n")

# ==============================================================================
# SPLIT INTO GREY vs COLORED POINTS
# Grey: p_value >= 0.05 OR mean_log2FC_Young < 0 OR mean_log2FC_Senescent < 0
# ==============================================================================
grey_mask <- df$p_value >= p_threshold |
  df$mean_log2FC_Young < 0 |
  df$mean_log2FC_Senescent < 0 |
  abs(df$delta_log2FC) < delta_fc_threshold
df_grey <- df[grey_mask, ]
df_colored <- df[!grey_mask, ]

n_up <- sum(df$regulation == "increased_in_senescent")
n_down <- sum(df$regulation == "decreased_in_senescent")

cat("Grey points (p>=0.05 or mean<0):", nrow(df_grey), "\n")
cat("Colored points:", nrow(df_colored), "\n")

# ==============================================================================
# LABELS: top 5 increased + top 5 decreased by p_value (among significant)
# ==============================================================================
sig_df <- df[df$regulation != "no_change", ]
label_up <- sig_df[sig_df$regulation == "increased_in_senescent", ]
label_up <- label_up[order(label_up$p_value), ][
  seq_len(min(5, nrow(label_up))),
]
label_dn <- sig_df[sig_df$regulation == "decreased_in_senescent", ]
label_dn <- label_dn[order(label_dn$p_value), ][
  seq_len(min(5, nrow(label_dn))),
]
label_df <- rbind(label_up, label_dn)

# ==============================================================================
# ZOOM: clamp x/y to ±xy_zoom, mark out-of-range with diamond symbol
# ==============================================================================
clamp_xy <- function(d, lim, size_inrange) {
  d$x_plot <- pmin(pmax(d$mean_log2FC_Young, -lim * 0.975), lim * 0.975)
  d$y_plot <- pmin(pmax(d$mean_log2FC_Senescent, -lim * 0.975), lim * 0.975)
  d$out_range <- abs(d$mean_log2FC_Young) > lim |
    abs(d$mean_log2FC_Senescent) > lim
  d$symbol <- ifelse(d$out_range, "diamond", "circle")
  d$pt_size <- ifelse(d$out_range, 2, size_inrange) # out-of-range smallest
  d
}
df_grey <- clamp_xy(df_grey, xy_zoom, size_inrange = 4) # grey: medium
df_colored <- clamp_xy(df_colored, xy_zoom, size_inrange = 6) # colored: largest

clamp1 <- function(x, lim) pmin(pmax(x, -lim * 0.975), lim * 0.975)
label_df$x_plot <- clamp1(label_df$mean_log2FC_Young, xy_zoom)
label_df$y_plot <- clamp1(label_df$mean_log2FC_Senescent, xy_zoom)

# ==============================================================================
# HOVER TEXT
# ==============================================================================
make_hover <- function(d) {
  paste0(
    "<b>",
    d$GeneName,
    "</b> (",
    d$Accession,
    ")<br>",
    "Young log2FC: ",
    round(d$mean_log2FC_Young, 3),
    "<br>",
    "Senescent log2FC: ",
    round(d$mean_log2FC_Senescent, 3),
    "<br>",
    "delta_log2FC: ",
    round(d$delta_log2FC, 3),
    "<br>",
    "p-value: ",
    signif(d$p_value, 3),
    "<br>",
    "Regulation: ",
    d$regulation
  )
}

# ==============================================================================
# DIVERGING COLOR SCALE (blue=min, white=0, red=max), centered at 0
# ==============================================================================
max_delta <- max(abs(df_colored$delta_log2FC), na.rm = TRUE)
if (max_delta == 0) {
  max_delta <- 1
} # guard against edge case

# Map 0 to midpoint of the color scale
mid_frac <- (-(-max_delta)) / (2 * max_delta) # = 0.5 always when symmetric

colorscale_diverging <- list(
  list(0, "#3498DB"), # blue  (most negative)
  list(mid_frac, "white"), # white (zero)
  list(1, "#E74C3C") # red   (most positive)
)

# ==============================================================================
# BUILD PLOT
# ==============================================================================
fig <- plot_ly() %>%

  # --- Trace 1: Grey points (clamped x/y, diamond if out-of-range) ---
  add_trace(
    data = df_grey,
    type = "scatter3d",
    mode = "markers",
    x = ~x_plot,
    y = ~y_plot,
    z = ~neg_log10_p,
    marker = list(
      color = "lightgrey",
      symbol = ~symbol,
      size = ~pt_size,
      opacity = 0.5,
      line = list(width = 0)
    ),
    text = make_hover(df_grey),
    hoverinfo = "text",
    name = "p >= 0.05 or mean < 0",
    showlegend = TRUE
  ) %>%

  # --- Trace 2: Colored points (delta_log2FC gradient, clamped x/y) ---
  add_trace(
    data = df_colored,
    type = "scatter3d",
    mode = "markers",
    x = ~x_plot,
    y = ~y_plot,
    z = ~neg_log10_p,
    marker = list(
      color = ~delta_log2FC,
      colorscale = colorscale_diverging,
      cmin = -max_delta,
      cmax = max_delta,
      symbol = ~symbol,
      size = ~pt_size,
      opacity = 0.85,
      line = list(width = 0),
      colorbar = list(
        title = "delta_log2FC",
        titleside = "right",
        tickformat = ".1f",
        tickfont = list(size = 14),
        titlefont = list(size = 14),
        len = 0.6, # colorbar height as fraction of plot height
        x = 0.9 # pull colorbar closer to the 3D scene
      )
    ),
    text = make_hover(df_colored),
    hoverinfo = "text",
    name = "Significant (colored by delta_log2FC)",
    showlegend = TRUE
  ) %>%

  # --- Trace 3: Gene name labels (clamped positions) ---
  add_trace(
    data = label_df,
    type = "scatter3d",
    mode = "text",
    x = ~x_plot,
    y = ~y_plot,
    z = ~neg_log10_p,
    text = ~GeneName,
    textfont = list(size = 10, color = "black"),
    hoverinfo = "skip",
    name = "Labels",
    showlegend = FALSE
  ) %>%

  layout(
    title = list(
      text = paste0(
        "3D: Young vs Senescent FC  |  Z = -log10(p-value)<br>",
        "<sup>Increased: ",
        n_up,
        "   Decreased: ",
        n_down,
        "   (grey: p\u22650.05 or mean log2FC < 0)</sup>"
      ),
      font = list(size = 18),
      x = 0.5,
      xanchor = "center"
    ),
    scene = list(
      xaxis = list(
        title = "mean log2FC (Young)",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        range = c(-xy_zoom, xy_zoom)
      ),
      yaxis = list(
        title = "mean log2FC (Senescent)",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        range = c(-xy_zoom, xy_zoom)
      ),
      zaxis = list(
        title = "-log10(p-value)",
        titlefont = list(size = 14),
        tickfont = list(size = 12)
      ),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 0.8))
    ),
    legend = list(
      x = 0.02,
      y = 0.55, # move legend down to sit beside the scene
      font = list(size = 14),
      bgcolor = "rgba(255,255,255,0.7)",
      bordercolor = "lightgrey",
      borderwidth = 1
    ),
    margin = list(l = 0, r = 120, t = 80, b = 0)
  )

# ==============================================================================
# SAVE
# ==============================================================================
html_path <- file.path(output_dir, "delta_fc_3d.html")
saveWidget(
  fig,
  file = normalizePath(html_path, mustWork = FALSE),
  selfcontained = TRUE
)
cat("Saved: delta_fc_3d.html\n")

# --- PDF + PNG via kaleido (no browser/network required) ---
# Requires: install.packages("reticulate") + conda install -c conda-forge python-kaleido
# If this hangs: open delta_fc_3d.html in Chrome and use the camera toolbar button to save PNG.
# tryCatch({
#   Sys.setenv(RETICULATE_PYTHON = "D:/anaconda3/python.exe")
#   reticulate::use_condaenv("base", required = TRUE)
#   # Use forward slashes: normalizePath() returns backslashes on Windows which
#   # Python misinterprets as unicode escapes (e.g. \N, \r) inside string literals.
#   to_fwd <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))
#   plotly::save_image(fig,
#     file   = to_fwd(file.path(output_dir, "delta_fc_3d.png")),
#     width  = 1800, height = 1500)
#   cat("Saved: delta_fc_3d.png\n")
#   plotly::save_image(fig,
#     file   = to_fwd(file.path(output_dir, "delta_fc_3d.pdf")),
#     width  = 900, height = 750)
#   cat("Saved: delta_fc_3d.pdf\n")
# }, error = function(e) {
#   message("Static export failed: ", conditionMessage(e))
#   message("Tip: open delta_fc_3d.html in Chrome and use the camera icon to save PNG.")
# })
