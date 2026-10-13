library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Load data from CSV
res_combined <- read.csv("output/tables/deg_interaction_20261009.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(padj_interaction))

# 2. Filter out extreme numerical boundary artifacts (clamped LFC caps)
plot_df <- res_combined %>%
  filter(abs(LFC_Capillary) <= 8 & abs(LFC_LEC) <= 8)

# 3. Define key target genes to label
top_label_genes <- c("DCDC2", "JCHAIN", "DCHS1", "NACC1", "ZNF316", "TAF10", "SMC6", "LDLRAD2")

# 4. Generate Publication-Ready Plot
plot_scatter_final <- ggplot(
  plot_df, 
  aes(x = LFC_Capillary, y = LFC_LEC, color = resilience_type)
) +
  geom_point(alpha = 0.45, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(values = c(
    "Preserved/Buffered in LECs (Protected from Loss)" = "#2b83ba",
    "Blunted/Protected in LECs (Protected from Surge)" = "#abdda4",
    "Hyper-induced in LECs"                             = "#d7191c",
    "Hyper-suppressed in LECs"                          = "#fdae61",
    "Moderate / Unclassified Shift"                    = "grey80"
  )) +
  geom_text_repel(
    data = filter(plot_df, gene %in% top_label_genes),
    aes(label = gene),
    color = "black",
    fontface = "bold",
    size = 4,
    box.padding = 0.8,
    point.padding = 0.1,
    force = 4,
    max.overlaps = 30
  ) +
  coord_fixed(ratio = 1, xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  theme_classic() +
  labs(
    title = "Endothelial Aging Divergence: Capillary vs. Lymphatic",
    # subtitle = "Symmetric 1:1 Scale (Excluding Sparse Boundary Artifacts)",
    x = "Capillary Aging LFC \n (Aged vs. Young BEC)",
    y = "Lymphatic Aging LFC \n (Aged vs. Young LEC)",
    color = "Resilience Category"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.5),
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 13),
    # Add explicit outer margins to safeguard against display clipping
    plot.margin = margin(t = 10, r = 25, b = 10, l = 25, unit = "pt")
  ) +
  # Center legend title on top and arrange items in 2 clean columns
  guides(
    color = guide_legend(
      ncol = 2, 
      byrow = TRUE, 
      title.position = "top", 
      title.hjust = 0.5
    )
  )

print(plot_scatter_final)

# Export at optimal width to maintain layout spacing
# ggsave("output/figures/scatter_resilience_final.png", plot_scatter_final, width = 8, height = 8.5, dpi = 300)