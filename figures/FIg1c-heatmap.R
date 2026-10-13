library(dplyr)
library(pheatmap)
library(tibble)
library(grid)
library(svglite)
library(ggplot2)
library(gridExtra)

# 1. Release any lingering graphics device file handles
graphics.off()

# 2. Load interaction results from CSV
res_combined <- read.csv("output/tables/deg_interaction_20261009.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(padj_interaction))

# 3. Extract top 30 significant genes and select Capillary vs. Lymphatic LFC columns
top30_genes <- res_combined %>%
  slice_min(order_by = padj_interaction, n = 30)

heatmap_matrix <- top30_genes %>%
  column_to_rownames(var = "gene") %>%
  select(LFC_Capillary, LFC_LEC)

colnames(heatmap_matrix) <- c("Capillary Aging", "Lymphatic Aging")

# 4. PRINT INTERACTIVELY TO THE PLOTS PANE
ph <- pheatmap(
  as.matrix(heatmap_matrix),
  scale = "none",
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
  breaks = seq(-4, 4, length.out = 101),  # Center symmetric color scale at 0
  cluster_cols = FALSE,                  # Keep Capillary vs Lymphatic order fixed
  cluster_rows = TRUE,                   # Group genes with similar aging trends
  display_numbers = FALSE,               # Clean boxes without text overlay
  angle_col = "0",                       # Keep column titles horizontal
  fontsize_row = 9.5,
  fontsize_col = 11,
  main = "Top 30 Lineage-Resilient Endothelial Genes",
  border_color = "white"                 # Crisp grid lines
)

# Display in active session
print(ph)
# ph_gtable <- ph[[4]]
# 
# ggsave(
#   "output/figures/heatmap_top30_clean.svg",
#   plot = grid.arrange(ph_gtable, nrow = 1, ncol = 1),
#   width = 5.5,
#   height = 7.5,
#   units = "in"
#   )

# 5. SAVE TO FILE SAFELY (Optional)
# Uncomment the block below to export after reviewing in the Plots pane
# png("output/figures/heatmap_top30_clean.png", width = 5.5, height = 7.5, units = "in", res = 300)
# grid.newpage()
# grid.draw(ph$gtable)
# dev.off()