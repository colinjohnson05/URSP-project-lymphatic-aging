library(dplyr)
library(ggplot2)
library(stringr)

# 1. Load GSEA results from CSV
gsea_results <- read.csv("output/tables/gsea_go_bp_interaction.csv", stringsAsFactors = FALSE)

# 2. Filter and apply strict string formatting
gsea_plot_df <- gsea_results %>%
  filter(!is.na(padj) & padj < 0.10) %>%
  mutate(
    # Strip prefix and convert to title case
    pathway_clean = gsub("^GOBP_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    pathway_clean = str_to_title(pathway_clean),
    
    # Restore correct uppercase capitalization for genetic acronyms
    pathway_clean = gsub("\\bMrna\\b", "mRNA", pathway_clean),
    pathway_clean = gsub("\\bRna\\b", "RNA", pathway_clean),
    
    # Wrap text cleanly across multiple lines
    pathway_clean = str_wrap(pathway_clean, width = 38),
    Direction = ifelse(NES > 0, "Preserved / Higher in LEC Aging", "Blunted / Lower in LEC Aging")
  ) %>%
  arrange(desc(NES))

# 3. Render Publication-Ready Plot
plot_gsea_bp <- ggplot(gsea_plot_df, aes(x = reorder(pathway_clean, NES), y = NES, fill = Direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Preserved / Higher in LEC Aging" = "#2b83ba",
    "Blunted / Lower in LEC Aging"   = "#d7191c"
  )) +
  theme_classic() +
  labs(
    title = "Top Enriched GO Biological Processes",
    subtitle = "Lineage Resilience Interaction GSEA (LEC vs. BEC Aging)",
    x = "",
    y = "Normalized Enrichment Score (NES)",
    fill = "Trajectory"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 8.5, color = "black"),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15, unit = "pt")
  )

print(plot_gsea_bp)

# Export finalized figure
ggsave("output/figures/gsea_go_bp_final.svg", plot = plot_gsea_bp, width = 7.5, height = 6.0, units = "in")
