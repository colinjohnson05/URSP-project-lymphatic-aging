# This script loads data, merges atlases, and outputs merged atlas
# Load necessary packages
library(qs2)
library(Seurat)
library(anndata)
library(dplyr)
library(ggplot2)
library(patchwork)
library(devil)
library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(SummarizedExperiment)

# Load datasets
lymph <- qs_read("data/Lymphatics_processed.qs2")
lung <- readH5AD("data/lungMergedEndo10.h5ad")
lymph
lymph_counts <- LayerData(lymph, assay = "RNA", layer = "counts")
lymph_meta <- lymph@meta.data
lung
lung_counts <- assay(lung, "X")
lung_meta <- as.data.frame(colData(lung))

# Figure out how cells are labeled
# Print unique cell annotations from each dataset metadata
cat("Lymphatics_processed unique cell type annotations:\n ", unique(lymph_meta$celltype_level3_fullname))
cat("lungMergedEndo10 unique cell type annotations:\n", unique(lung_meta$my_subtype_annotations))
# Convert "my_subtype_annotations" in lung_meta to plaintext characters
lung_meta$my_subtype_annotations_char <- as.character(lung_meta$my_subtype_annotations)
cat("lungMergedEndo10 unique cell type annotations:\n", unique(lung_meta$my_subtype_annotations_char))

### Merge datasets
# Define max mitochondrial threshold
max_pmt <- 10.0

# 1A. Standardize Lymph Dataset
colnames(lymph_counts) <- paste0("lymph_", colnames(lymph_counts))
rownames(lymph_meta)   <- paste0("lymph_", rownames(lymph_meta))

lymph_meta$dataset               <- "lymph"
lymph_meta$cell_type_annotation  <- lymph_meta$celltype_level3_fullname
lymph_meta$is_lymphatic          <- lymph_meta$celltype_level3_fullname == "Lymphatic endothelial cell"
lymph_meta$donor_unique          <- lymph_meta$DonorID
lymph_meta$age                   <- lymph_meta$Age
lymph_meta$percent_mito          <- lymph_meta$pMT

# Apply QC filter to Lymph
keep_lymph_qc <- lymph_meta$percent_mito <= max_pmt & !is.na(lymph_meta$percent_mito)
lymph_counts  <- lymph_counts[, keep_lymph_qc, drop = FALSE]
lymph_meta    <- lymph_meta[keep_lymph_qc, ]

# 1B. Standardize Lung Dataset
colnames(lung_counts) <- paste0("lung_", colnames(lung_counts))
rownames(lung_meta)   <- paste0("lung_", rownames(lung_meta))

lung_meta$dataset               <- "lung"
lung_meta$cell_type_annotation  <- lung_meta$my_subtype_annotations_char
lung_meta$is_lymphatic          <- lung_meta$my_subtype_annotations_char == "Lymphatic"
lung_meta$donor_unique          <- lung_meta$donorID
lung_meta$percent_mito          <- lung_meta$perMito

# Apply QC filter to Lung
keep_lung_qc <- lung_meta$percent_mito <= max_pmt & !is.na(lung_meta$percent_mito)
lung_counts  <- lung_counts[, keep_lung_qc, drop = FALSE]
lung_meta    <- lung_meta[keep_lung_qc, ]

cat(sprintf("QC Filter (pMT <= %.1f%%) — Kept Lymph: %d cells | Kept Lung: %d cells\n",
            max_pmt, ncol(lymph_counts), ncol(lung_counts)))

# Find genes present in both matrices
common_genes <- intersect(rownames(lymph_counts), rownames(lung_counts))
cat("Common genes across datasets:", length(common_genes), "\n")

# Subset count matrices to common genes and combine (cbind)
overall_counts <- cbind(
  lymph_counts[common_genes, ],
  lung_counts[common_genes, ]
)

# Select key metadata columns needed for downstream DEVIL & Harmony analysis
# Note: Ensure 'age' and 'donor_unique' match your actual column names in lymph_meta/lung_meta
target_meta_cols <- c("dataset", "age", "donor_unique", "cell_type_annotation", "is_lymphatic")

# Helper function to align metadata columns smoothly
align_metadata <- function(df, cols) {
  missing_cols <- setdiff(cols, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  return(df[, cols, drop = FALSE])
}

overall_meta <- rbind(
  align_metadata(lymph_meta, target_meta_cols),
  align_metadata(lung_meta, target_meta_cols)
)

overall_data <- list(
  counts = overall_counts,
  meta   = overall_meta
)

# Save merged
qs_save(overall_data, "data/overall_merged_data.qs2")
cat("Saved overall merged dataset — Genes:", nrow(overall_counts), "Cells:", ncol(overall_counts), "\n")
