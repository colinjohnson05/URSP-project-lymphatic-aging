# Overall Integrated Filtering by Age, UMAP, & Lymphatic DEVIL Input Filtering

library(qs2)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(Matrix)


# ==============================================================================
# 1. Define Age Cutoffs and Add Categorical Age Group
# ==============================================================================

young_cutoff <- 40
old_cutoff   <- 60

overall_meta <- overall_meta %>%
  mutate(
    age_group = case_when(
      age <= young_cutoff ~ "Young",
      age >= old_cutoff   ~ "Old",
      TRUE                ~ NA_character_ # Middle-aged donors and missing ages become NA
    )
  )


# ==============================================================================
# 2. Subset and Save Overall Dataset (Filtered by Age: Young & Old Only)
# ==============================================================================

# Keep cells that have a valid age group (drops NAs and middle-aged cells)
overall_keep_idx <- !is.na(overall_meta$age_group)

overall_age_counts <- overall_counts[, overall_keep_idx, drop = FALSE]
overall_age_meta   <- overall_meta[overall_keep_idx, ]

overall_age_data <- list(
  counts = overall_age_counts,
  meta   = overall_age_meta
)

qs_save(overall_age_data, "data/overall_aged_young_data.qs2")
cat(sprintf("Saved age-filtered overall dataset — Genes: %d, Cells: %d\n", 
            nrow(overall_age_counts), ncol(overall_age_counts)))


# ==============================================================================
# 3. Subset and Save Lymphatics Only (Lymphatic + Valid Age Group)
# ==============================================================================

lymphatic_keep_idx <- overall_meta$is_lymphatic & !is.na(overall_meta$age_group)

lymphatic_counts <- overall_counts[, lymphatic_keep_idx, drop = FALSE]
lymphatic_meta   <- overall_meta[lymphatic_keep_idx, ]

lymphatic_data <- list(
  counts = lymphatic_counts,
  meta   = lymphatic_meta
)

qs_save(lymphatic_data, "data/lymphatics_merged_data.qs2")
cat(sprintf("Saved lymphatics merged dataset — Genes: %d, Cells: %d\n", 
            nrow(lymphatic_counts), ncol(lymphatic_counts)))

# ==============================================================================
# PART 2: Lymphatics 1% Gene Prevalence Filtering for DEVIL
# ==============================================================================

cat("--- 2. Filtering Lymphatic Subset for DEVIL ---\n")

# 1. Load raw lymphatic dataset
lymph_data <- qs_read("data/lymphatics_merged_data.qs2")
counts_raw <- lymph_data$counts
meta_lymph <- lymph_data$meta

# 2. Calculate 1% cell detection threshold (count >= 1)
min_cells <- max(1, floor(0.01 * ncol(counts_raw)))

# 3. Identify genes expressed in >= 1% of lymphatic cells
keep_genes <- Matrix::rowSums(counts_raw >= 1) >= min_cells
counts_devil_ready <- counts_raw[keep_genes, , drop = FALSE]

cat(sprintf("Prevalence Filtering Results (expressed in >= %d cells / 1%%):\n", min_cells))
cat(sprintf("  • Initial genes : %d\n", nrow(counts_raw)))
cat(sprintf("  • Retained genes: %d\n", nrow(counts_devil_ready)))
cat(sprintf("  • Filtered genes: %d\n\n", nrow(counts_raw) - nrow(counts_devil_ready)))

# 4. Save analysis-ready input file for DEVIL
devil_input <- list(
  counts = counts_devil_ready, # Sparse raw integer counts matrix
  meta   = meta_lymph          # Unmodified metadata containing 'age' and 'donor_unique'
)

qs_save(devil_input, "data/lymphatics_for_devil.qs2")
cat("Saved ready-for-DEVIL object to 'data/lymphatics_for_devil.qs2'\n")