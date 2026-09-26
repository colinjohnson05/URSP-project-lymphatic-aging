# Create output directories if they do not exist
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# 1. Generate synthetic cell atlas data
set.seed(42)
n_cells <- 600

mock_atlas <- data.frame(
  UMAP_1 = c(rnorm(200, -3), rnorm(200, 3), rnorm(200, 0)),
  UMAP_2 = c(rnorm(200, -3), rnorm(200, 3), rnorm(200, 4)),
  dataset = sample(c("Dataset_A", "Dataset_B", "Dataset_C"), n_cells, replace = TRUE),
  cell_type = rep(c("T-cell", "B-cell", "Monocyte"), each = 200)
)

# 2. Save mock summary table
summary_tbl <- aggregate(UMAP_1 ~ dataset + cell_type, data = mock_atlas, FUN = length)
colnames(summary_tbl)[3] <- "cell_count"
write.csv(summary_tbl, "output/tables/atlas_summary.csv", row.names = FALSE)

# 3. Generate and save mock UMAP plot (using base R graphics)
png("output/figures/umap_integrated.png", width = 800, height = 600, res = 100)

colors <- c("T-cell" = "#E41A1C", "B-cell" = "#377EB8", "Monocyte" = "#4DAF4A")
plot(
  mock_atlas$UMAP_1, mock_atlas$UMAP_2,
  col = colors[mock_atlas$cell_type],
  pch = 16, cex = 0.8,
  xlab = "UMAP 1", ylab = "UMAP 2",
  main = "Simulated Atlas Integration (UMAP)"
)
legend("topright", legend = names(colors), col = colors, pch = 16, title = "Cell Type")

dev.off()

cat("Pipeline run complete. Outputs saved to output/figures/ and output/tables/\n")
