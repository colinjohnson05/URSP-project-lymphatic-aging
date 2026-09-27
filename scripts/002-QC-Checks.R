# This script plots the overall UMAP to make sure that cells aren't clustering
# by donor or dataset

# ==============================================================================
# Diagnostic Script: Check Dataset & Donor Batch Effects (PCA & UMAP)
# ==============================================================================

library(qs2)
library(Seurat)
library(ggplot2)
library(patchwork)

# 1. Load overall merged dataset
data_list <- qs_read("data/overall_merged_data.qs2")
counts    <- data_list$counts
meta      <- data_list$meta

# 2. Create Seurat object and run standard uncorrected pipeline
seu <- CreateSeuratObject(counts = counts, meta.data = meta)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20)


# ==============================================================================
# 3. Generate Visual Diagnostic Plots
# ==============================================================================

# --- PCA Diagnostics (Checking primary axes of variance) ---
p_pca_dataset <- DimPlot(seu, reduction = "pca", dims = c(1, 2), group.by = "dataset") +
  ggtitle("PCA (PC1 vs PC2) — Dataset") +
  theme_minimal()

p_pca_donor <- DimPlot(seu, reduction = "pca", dims = c(1, 2), group.by = "donor_unique") +
  ggtitle("PCA (PC1 vs PC2) — Donor") +
  theme_minimal() +
  theme(legend.position = "none") # Suppress legend if donor count is large

# --- UMAP Diagnostics (Checking overall cluster mixing) ---
p_umap_dataset <- DimPlot(seu, reduction = "umap", group.by = "dataset") +
  ggtitle("UMAP — Colored by Dataset") +
  theme_minimal()

p_umap_donor <- DimPlot(seu, reduction = "umap", group.by = "donor_unique") +
  ggtitle("UMAP — Colored by Donor") +
  theme_minimal() +
  theme(legend.position = "none") # Suppress legend if donor count is large

p_umap_celltype <- DimPlot(seu, reduction = "umap", group.by = "cell_type_annotation") +
  ggtitle("UMAP — Colored by Cell Type") +
  theme_minimal()


# ==============================================================================
# 4. Assemble and Save Diagnostic Figure
# ==============================================================================

# Layout: PCA row on top, UMAP row on bottom
qc_figure <- (p_pca_dataset | p_pca_donor) / 
  (p_umap_dataset | p_umap_donor | p_umap_celltype)

ggsave("output/figures/batch_qc_pca_umap.png", plot = qc_figure, width = 16, height = 10, dpi = 300)

print("QC figure saved as 'output/figures/batch_qc_pca_umap.png'")