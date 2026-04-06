# ==============================================================================
# SCATTER FC PLOT: Young vs Senescent mean log2FC (Standalone)
# ==============================================================================
# Usage:
#   setwd("d:/ResearchProject/NPM1Pulldown_TianxiangWang")
#   source("plots/scatter_fc.R")
#
# Output: res_pulldown/delta_fc_results/scatter_fc_standalone.pdf
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(scales)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
delta_fc_threshold <- 1.0 # log2 scale (must match main_pulldown.R)
xy_zoom <- 5 # |mean_log2FC| > xy_zoom shown as triangles at boundary

COLORS <- c(
  "increased_in_senescent" = "#E74C3C",
  "decreased_in_senescent" = "#3498DB",
  "no_change" = "#95A5A6"
)
COLOR_LABELS <- c(
  "increased_in_senescent" = "Increased in Senescent",
  "decreased_in_senescent" = "Decreased in Senescent",
  "no_change" = "No Change"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================
data_path <- "res_pulldown/delta_fc_results/delta_fc_results.csv"
output_dir <- "res_pulldown/delta_fc_results"

if (!file.exists(data_path)) {
  stop("Data not found: ", data_path)
}

df <- read.csv(data_path, stringsAsFactors = FALSE)

cat("Loaded", nrow(df), "proteins\n")
cat("Regulation summary:\n")
print(table(df$regulation))

# ==============================================================================
# ANNOTATIONS: top 20 significant genes by |delta_log2FC|
# ==============================================================================
sig_df <- df[df$regulation != "no_change", ]
top20 <- sig_df[order(abs(sig_df$delta_log2FC), decreasing = TRUE), ][
  seq_len(min(20, nrow(sig_df))),
]

# ==============================================================================
# PLOT
# ==============================================================================
xy_lim <- c(-xy_zoom, xy_zoom)

# A point is "in range" only if BOTH axes are within zoom limits
in_range <- abs(df$mean_log2FC_Young) <= xy_zoom &
  abs(df$mean_log2FC_Senescent) <= xy_zoom
df_in <- df[in_range, ]
df_out <- df[!in_range, ]

# Clamp out-of-range coordinates to just inside boundary
clamp <- function(x, lim) pmin(pmax(x, -lim * 0.975), lim * 0.975)
df_out$x_plot <- clamp(df_out$mean_log2FC_Young, xy_zoom)
df_out$y_plot <- clamp(df_out$mean_log2FC_Senescent, xy_zoom)

scatter_plot <- ggplot(mapping = aes(color = regulation)) +
  # Normal points (circles) — no_change first (background)
  geom_point(
    data = df_in[df_in$regulation == "no_change", ],
    aes(x = mean_log2FC_Young, y = mean_log2FC_Senescent),
    shape = 16,
    alpha = 0.5,
    size = 1.2
  ) +
  geom_point(
    data = df_in[df_in$regulation != "no_change", ],
    aes(x = mean_log2FC_Young, y = mean_log2FC_Senescent),
    shape = 16,
    alpha = 0.8,
    size = 2.0
  ) +
  # Out-of-range points (triangles at boundary)
  {
    if (nrow(df_out) > 0) {
      geom_point(
        data = df_out,
        aes(x = x_plot, y = y_plot),
        shape = 17,
        size = 0.8,
        alpha = 0.85
      )
    }
  } +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = delta_fc_threshold,
    linetype = "dotted",
    color = "gray40",
    linewidth = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = -delta_fc_threshold,
    linetype = "dotted",
    color = "gray40",
    linewidth = 0.5
  ) +
  {
    if (nrow(top20) > 0) {
      # Use clamped coords for labels of out-of-range genes
      top20$x_label <- clamp(top20$mean_log2FC_Young, xy_zoom)
      top20$y_label <- clamp(top20$mean_log2FC_Senescent, xy_zoom)
      geom_text_repel(
        data = top20,
        aes(x = x_label, y = y_label, label = GeneName),
        size = 2.5,
        color = "black",
        bg.color = "white",
        bg.r = 0.08,
        segment.color = "gray30",
        segment.size = 0.3,
        max.overlaps = 20,
        seed = 42
      )
    }
  } +
  scale_color_manual(
    values = COLORS,
    labels = COLOR_LABELS,
    name = "Regulation"
  ) +
  scale_x_continuous(
    name = expression("mean log"[2] * "FC (Young NPM1/IgG)"),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = 0, add = 0)
  ) +
  scale_y_continuous(
    name = expression("mean log"[2] * "FC (Senescent NPM1/IgG)"),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = 0, add = 0)
  ) +
  coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
  labs(
    title = "Pairwise FC Comparison: Young vs Senescent",
    subtitle = paste0(
      "Dotted lines: y = x \u00b1 ",
      delta_fc_threshold,
      " (delta_log2FC threshold)"
    )
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray60"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave(
  filename = file.path(output_dir, "scatter_fc_standalone.pdf"),
  plot = scatter_plot,
  width = 7,
  height = 5
)
cat("Saved: scatter_fc_standalone.pdf\n")
