library(decoupleR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Verify and Align Cell Barcodes
# -----------------------------------------------------------------------------
DefaultAssay(seurat_obj) <- "RNA"

# Get exact cell names from Seurat
seurat_cells <- colnames(seurat_obj)

# Ensure scale.data exists
if (!"scale.data" %in% slotNames(seurat_obj[["RNA"]]) || nrow(GetAssayData(seurat_obj, slot = "scale.data")) == 0) {
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
}

mat_scaled <- GetAssayData(seurat_obj, layer = "scale.data")

# -----------------------------------------------------------------------------
# 2. Rerun decoupleR ULM
# -----------------------------------------------------------------------------
tf_activities <- run_ulm(
  mat = mat_scaled,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# -----------------------------------------------------------------------------
# 3. Reconstruct Matrix with Explicit Barcode Alignment
# -----------------------------------------------------------------------------
tf_matrix <- tf_activities %>%
  filter(statistic == "ulm") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score, values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()

# Subset and reorder columns to match Seurat EXACTLY
common_cells <- intersect(seurat_cells, colnames(tf_matrix))
message("Matching cells: ", length(common_cells), " out of ", length(seurat_cells))

tf_matrix_aligned <- matrix(
  0, 
  nrow = nrow(tf_matrix), 
  ncol = length(seurat_cells),
  dimnames = list(rownames(tf_matrix), seurat_cells)
)

tf_matrix_aligned[, common_cells] <- tf_matrix[, common_cells]

# -----------------------------------------------------------------------------
# 4. Overwrite tf_ulm Assay
# -----------------------------------------------------------------------------
seurat_obj[["tf_ulm"]] <- CreateAssayObject(data = tf_matrix_aligned)
DefaultAssay(seurat_obj) <- "tf_ulm"

# -----------------------------------------------------------------------------
# 5. Plot with Midpoint Divergent Gradient
# -----------------------------------------------------------------------------
# Calculate symmetric score limits for clean plotting
max_val <- max(abs(range(tf_matrix_aligned[c("MEF2C", "GRHL2", "PROX1", "EHF"), ])))

plot_tf_umap <- FeaturePlot(
  seurat_obj, 
  features = c("MEF2C", "GRHL2", "PROX1", "EHF"),
  reduction = "umap.integrated",
  cols = c("#2b83ba", "#f7f7f7", "#d7191c"), # Blue (Repressed/Low), White (Neutral), Red (High)
  keep.scale = "all",
  ncol = 2
) & 
  theme_classic() &
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 8)
  )

print(plot_tf_umap)