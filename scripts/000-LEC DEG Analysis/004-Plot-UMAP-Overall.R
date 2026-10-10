# ==============================================================================
# Script: Harmony Integration and UMAP Visualization (Young vs. Old)
# ==============================================================================

library(qs2)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# 1. Load age-filtered overall dataset
data_list <- qs_read("data/overall_aged_young_data.qs2")
counts    <- data_list$counts
meta      <- data_list$meta

# Ensure age_group is a factor with explicit level order
meta$age_group <- factor(meta$age_group, levels = c("Young", "Old"))

# 2. Create Seurat object and run standard preprocessing
seu <- CreateSeuratObject(counts = counts, meta.data = meta)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# 3. Batch correction across donors using Harmony
seu <- RunHarmony(
  seu, 
  group.by.vars  = "donor_unique", 
  reduction.save = "harmony", 
  verbose        = FALSE
)

# 4. Generate UMAP embedding on Harmony dimensions
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, reduction.name = "umap_harmony")

# ==============================================================================
# 5. Generate Diagnostic Visualizations
# ==============================================================================

# Custom palette for clear age contrast
age_colors <- c("Young" = "#3182bd", "Old" = "#de2d26")

# UMAP 1: Age Group (Young vs Old)
p_age <- DimPlot(
  seu, 
  reduction = "umap_harmony", 
  group.by  = "age_group", 
  cols      = age_colors, 
  pt.size   = 0.3
) +
  ggtitle("Endothelial Landscape — Young vs. Old") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# UMAP 2: Cell Type Annotations
p_celltype <- DimPlot(
  seu, 
  reduction = "umap_harmony", 
  group.by  = "cell_type_annotation", 
  pt.size   = 0.3
) +
  ggtitle("Cell Type Annotations") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# UMAP 3: Dataset Origin
p_dataset <- DimPlot(
  seu, 
  reduction = "umap_harmony", 
  group.by  = "dataset", 
  pt.size   = 0.3
) +
  ggtitle("Dataset Origin (Lung vs. Lymph)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

# Combine into a clean multi-panel figure
final_umap_panel <- (p_age | p_celltype) / p_dataset

# Save figure
ggsave("output/figures/overall_harmony_umap_age_group.png", plot = final_umap_panel, width = 14, height = 11, dpi = 300)

cat("Integrated UMAP figure saved as 'output/figures/overall_harmony_umap_age_group.png'\n")