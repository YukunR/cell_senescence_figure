# ==============================================================================
# VOLCANO PLOT: Delta FC (Standalone)
# ==============================================================================
# Usage:
#   setwd("d:/ResearchProject/NPM1Pulldown_TianxiangWang")
#   source("plots/volcano_delta_fc.R")
#
# Output: res_pulldown/delta_fc_results/volcano_delta_fc_standalone.pdf
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(scales)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
delta_fc_threshold <- 1.0 # log2 scale (must match main_pulldown.R)
p_threshold <- 0.05
x_zoom <- 5 # |delta_log2FC| > x_zoom shown as triangles at boundary

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
df$neg_log10_p <- -log10(df$p_value)

cat("Loaded", nrow(df), "proteins\n")
cat("Regulation summary:\n")
print(table(df$regulation))

# ==============================================================================
# ANNOTATIONS: top 10 up + top 10 down by |delta_log2FC|
# ==============================================================================
sig_df <- df[df$regulation != "no_change", ]
top_up <- sig_df[sig_df$regulation == "increased_in_senescent", ]
top_up <- top_up[order(top_up$delta_log2FC, decreasing = TRUE), ][
  seq_len(min(10, nrow(top_up))),
]
top_down <- sig_df[sig_df$regulation == "decreased_in_senescent", ]
top_down <- top_down[order(top_down$delta_log2FC), ][
  seq_len(min(10, nrow(top_down))),
]
annot_df <- rbind(top_up, top_down)

# ==============================================================================
# PLOT
# ==============================================================================
n_up <- sum(df$regulation == "increased_in_senescent")
n_down <- sum(df$regulation == "decreased_in_senescent")
n_ns <- sum(df$regulation == "no_change")

# Split into in-range and out-of-range points
df_in <- df[abs(df$delta_log2FC) <= x_zoom, ]
df_out <- df[abs(df$delta_log2FC) > x_zoom, ]

# Clamp out-of-range x to just inside the boundary (so triangles are visible)
df_out$x_plot <- sign(df_out$delta_log2FC) * x_zoom * 0.975

x_limits <- c(-x_zoom, x_zoom)

volcano_plot <- ggplot(mapping = aes(color = regulation)) +
  # Normal points (circles): no_change (medium) and significant (largest)
  geom_point(
    data = df_in[df_in$regulation == "no_change", ],
    aes(x = delta_log2FC, y = neg_log10_p),
    shape = 16,
    size = 1.2,
    alpha = 0.7
  ) +
  geom_point(
    data = df_in[df_in$regulation != "no_change", ],
    aes(x = delta_log2FC, y = neg_log10_p),
    shape = 16,
    size = 2.0,
    alpha = 0.7
  ) +
  # Out-of-range points (triangles at boundary)
  {
    if (nrow(df_out) > 0) {
      geom_point(
        data = df_out,
        aes(x = x_plot, y = neg_log10_p),
        shape = 17,
        size = 0.8,
        alpha = 0.85
      )
    }
  } +
  geom_hline(
    yintercept = -log10(p_threshold),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(-delta_fc_threshold, delta_fc_threshold),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.5
  ) +
  {
    if (nrow(annot_df) > 0) {
      # Clamp label x positions for out-of-range annotated genes
      annot_df$x_label <- pmin(
        pmax(annot_df$delta_log2FC, -x_zoom * 0.975),
        x_zoom * 0.975
      )
      geom_text_repel(
        data = annot_df,
        aes(x = x_label, y = neg_log10_p, label = GeneName),
        size = 3,
        color = "black",
        bg.color = "white",
        bg.r = 0.1,
        segment.color = "gray30",
        segment.size = 0.3,
        min.segment.length = 0.1,
        max.overlaps = 20,
        force = 2,
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
    name = bquote(Delta * log[2] ~ "FC  [Senescent - Young]"),
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = 0, add = 0)
  ) +
  scale_y_continuous(
    name = expression(-log[10] * italic(p)),
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = 0, add = 0)
  ) +
  coord_cartesian(xlim = x_limits) +
  ggtitle(
    label = expression(Delta * log[2] * "FC Volcano Plot"),
    subtitle = bquote(
      Delta * "FC threshold: \u00b1" * .(delta_fc_threshold) ~
        " | P-value threshold:" ~ .(p_threshold)
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
  ) +
  annotate(
    "text",
    x = -x_zoom + 0.1 * 2 * x_zoom,
    y = max(df$neg_log10_p, na.rm = TRUE) * 0.95,
    label = paste0(
      "Increased: ",
      n_up,
      " (",
      round(n_up / nrow(df) * 100, 1),
      "%)\n",
      "Decreased: ",
      n_down,
      " (",
      round(n_down / nrow(df) * 100, 1),
      "%)\n",
      "No change: ",
      n_ns,
      " (",
      round(n_ns / nrow(df) * 100, 1),
      "%)"
    ),
    hjust = 0,
    vjust = 1,
    size = 3,
    color = "gray40",
    fontface = "italic"
  )

ggsave(
  filename = file.path(output_dir, "volcano_delta_fc_standalone.pdf"),
  plot = volcano_plot,
  width = 7,
  height = 5
)
cat("Saved: volcano_delta_fc_standalone.pdf\n")
