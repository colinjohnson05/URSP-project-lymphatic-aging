library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. Load data from CSV
res_combined <- read.csv("output/tables/deg_interaction_20261009.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(padj_interaction))

# 2. Assign 3-class significance categories
volcano_df <- res_combined %>%
  filter(abs(beta_interaction) <= 8) %>%
  mutate(
    volcano_category = case_when(
      padj_interaction < 0.05 & beta_interaction > 0  ~ "Preserved / Higher in LEC Aging",
      padj_interaction < 0.05 & beta_interaction < 0  ~ "Blunted / Lower in LEC Aging",
      TRUE                                             ~ "Not Significant (FDR ≥ 0.05)"
    ),
    volcano_category = factor(
      volcano_category, 
      levels = c("Not Significant (FDR ≥ 0.05)", "Preserved / Higher in LEC Aging", "Blunted / Lower in LEC Aging")
    )
  ) %>%
  arrange(volcano_category)

# 3. Select target genes for labels
key_target_genes <- c("DCDC2", "JCHAIN", "DCHS1", "LDLRAD2", "NACC1", "ZNF316", "TAF10", "SMC6")
label_df <- volcano_df %>% filter(gene %in% key_target_genes & padj_interaction < 0.05)

# 4. Generate Square Aspect Volcano Plot
plot_volcano_square <- ggplot(
  volcano_df, 
  aes(x = beta_interaction, y = -log10(padj_interaction), color = volcano_category)
) +
  geom_point(
    aes(alpha = volcano_category, size = volcano_category), 
    shape = 16
  ) +
  scale_alpha_manual(values = c("Not Significant (FDR ≥ 0.05)" = 0.25, "Preserved / Higher in LEC Aging" = 0.85, "Blunted / Lower in LEC Aging" = 0.85)) +
  scale_size_manual(values = c("Not Significant (FDR ≥ 0.05)" = 1.0, "Preserved / Higher in LEC Aging" = 2.0, "Blunted / Lower in LEC Aging" = 2.0)) +
  scale_color_manual(values = c(
    "Preserved / Higher in LEC Aging" = "#2b83ba",
    "Blunted / Lower in LEC Aging"   = "#d7191c",
    "Not Significant (FDR ≥ 0.05)"   = "#d9d9d9"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey50", linewidth = 0.5) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 3.8,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.1,
    segment.color = "grey30",
    segment.size = 0.4,
    force = 5,
    max.overlaps = 25
  ) +
  coord_cartesian(xlim = c(-6, 6)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  theme_classic() +
  labs(
    title = "Lineage Resilience Interaction Volcano Plot",
    # subtitle = "Divergence in Aging Trajectory (LEC vs. BEC)",
    x = expression("Interaction Coefficient (" * beta["interaction"] * ")"),
    y = expression("-log"[10] * "(Adjusted P-Value)"),
    color = "Significance & Direction"
  ) +
  theme(
    aspect.ratio = 1, # Sets 1:1 square aspect ratio for the plot canvas
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.5),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15, unit = "pt")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 3.5, alpha = 1),
      title.position = "top",
      title.hjust = 0.5
    ),
    alpha = "none",
    size = "none"
  )

print(plot_volcano_square)

# Save with matching square proportions
# ggsave("output/figures/volcano_square.png", plot_volcano_square, width = 6.5, height = 7.2, dpi = 300)