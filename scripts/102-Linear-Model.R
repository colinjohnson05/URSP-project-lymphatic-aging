### Use linear model with interaction terms to identify differentially
### expressed genes in aging lymphatics

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(tibble)
library(devil)
library(Matrix)

# Load data (if not in workspace already)
atlas_filtered <- readRDS("data/HAECA_filtered.RDS")

# Extract metadata
set.seed(123)
meta <- atlas_filtered@meta.data

# Set age bracket and lineage_broad as factors
meta$Age_brackets <- factor(
  meta$Age_brackets,
  levels = c("young", "aged")
  )

meta$lineage_broad <- factor(
  meta$lineage_broad,
  levels = c("Capillary", "Lymphatic")
  )

# Set up model design, using an interaction model
design <- model.matrix(
  ~ Age_brackets + lineage_broad + Age_brackets:lineage_broad,
  data = meta
)

head(design)
colnames(design)

# ---- Counts matrix -----------------------------------------------------------
# Join the two split layers and take the counts in OBJECT cell order.
#
# Do NOT use as.matrix() here. At ~17,000 genes x ~149,000 cells a dense double
# matrix is ~19 GB, and cbind() of two of them peaks near 38 GB on a 32 GB
# machine. devil accepts a sparse dgCMatrix directly (~2.5-4 GB).
#
# Joining also fixes a silent correctness bug: `design` and `meta$Donor` are in
# object cell order, but cbind(counts.1, counts.2) is in layer order. Those two
# orders are not guaranteed to agree, and if they disagree every gene's result
# is meaningless. LayerData() after JoinLayers() returns object order, and the
# assertions below make any mismatch fail loudly instead of silently.
atlas_filtered[["RNA"]] <- JoinLayers(atlas_filtered[["RNA"]])
Y_combined <- SeuratObject::LayerData(atlas_filtered, assay = "RNA", layer = "counts")

stopifnot(
  identical(colnames(Y_combined), rownames(meta)),  # cells line up with design
  nrow(design) == ncol(Y_combined)                  # model.matrix dropped no NA rows
)

# ---- Gene filter -------------------------------------------------------------
# Every gene carried into the test costs multiple-testing budget: BH compares
# ranked p-values against (i/m)*alpha, so thousands of barely-detected genes
# directly shrink the number of discoveries at a given FDR. Dropping genes that
# cannot support a meaningful estimate is a power gain, not a cherry-pick,
# provided the rule ignores the group labels -- this one does.
MIN_CELLS_DETECTED <- 10
keep_gene <- Matrix::rowSums(Y_combined > 0) >= MIN_CELLS_DETECTED
cat(sprintf("genes: %d -> %d retained (detected in >= %d cells)\n",
            length(keep_gene), sum(keep_gene), MIN_CELLS_DETECTED))
Y_combined <- Y_combined[keep_gene, ]

grouped_data <- devil::group_data(
  Y_combined,
  design,
  meta$Donor
  )

# # Fit model
devil_fit <- fit_devil(
  x              = grouped_data$input_matrix,
  design_matrix  = grouped_data$design_matrix,
  clusters       = grouped_data$clusters,
  verbose        = TRUE,
  init_beta_rough = TRUE,
  size_factors   = "normed_sum",
  overdispersion = "MOM",
  BPPARAM        = BiocParallel::SnowParam(workers = 4, progressbar = TRUE)
)

# ...
# ------------------------------------------------------------------------------
# 1. Define Contrast Vectors for Hypothesis Testing
# ------------------------------------------------------------------------------
# Extract column names from the grouped design matrix
design_cols <- colnames(grouped_data$design_matrix)
cat("Design matrix columns:\n", paste("-", design_cols, collapse = "\n"), "\n")

# Helper function to construct single-term contrast vectors
make_term_contrast <- function(design_matrix, term_name) {
  stopifnot(term_name %in% colnames(design_matrix))
  c_vec <- rep(0, ncol(design_matrix))
  names(c_vec) <- colnames(design_matrix)
  c_vec[term_name] <- 1
  return(c_vec)
}

# Contrast 1: Interaction Term (Age x Lineage resilience factor)
contrast_interaction <- make_term_contrast(
  grouped_data$design_matrix,
  "Age_bracketsaged:lineage_broadLymphatic"
)

# Contrast 2: Baseline Capillary Aging (Main Age Effect in Reference Lineage)
contrast_capillary_aging <- make_term_contrast(
  grouped_data$design_matrix,
  "Age_bracketsaged"
)

# ------------------------------------------------------------------------------
# 2. Execute DE Testing with DEVIL
# ------------------------------------------------------------------------------
cat("Testing Interaction term (Lineage Resilience)...\n")
res_interaction_raw <- devil::test_de(
  devil_fit,
  contrast = contrast_interaction
)

cat("Testing Baseline Capillary Aging term...\n")
res_capillary_raw <- devil::test_de(
  devil_fit,
  contrast = contrast_capillary_aging
)

# ------------------------------------------------------------------------------
# 1. Format Interaction Results
# ------------------------------------------------------------------------------
res_interaction <- as.data.frame(res_interaction_raw) %>%
  dplyr::rename(
    gene = name,
    beta_interaction = lfc,
    pval_interaction = pval,
    padj_interaction = adj_pval
  )

# ------------------------------------------------------------------------------
# 2. Format Baseline Capillary Aging Results
# ------------------------------------------------------------------------------
res_capillary <- as.data.frame(res_capillary_raw) %>%
  dplyr::rename(
    gene = name,
    LFC_Capillary = lfc,
    pval_Capillary = pval,
    padj_Capillary = adj_pval
  )

# ------------------------------------------------------------------------------
# 3. Merge Trajectories & Classify Resilience Mechanisms
# ------------------------------------------------------------------------------
res_combined <- res_interaction %>%
  inner_join(
    res_capillary %>% select(gene, LFC_Capillary, pval_Capillary, padj_Capillary),
    by = "gene"
  ) %>%
  mutate(
    # Compute derived aging trajectory in LECs (Capillary Baseline + Delta)
    LFC_LEC = LFC_Capillary + beta_interaction,
    
    # Categorize biological resilience shifts
    resilience_type = case_when(
      beta_interaction >  0.5 & LFC_Capillary < -0.2 ~ "Preserved/Buffered in LECs (Protected from Loss)",
      beta_interaction < -0.5 & LFC_Capillary >  0.2 ~ "Blunted/Protected in LECs (Protected from Surge)",
      beta_interaction >  0.5 & LFC_Capillary >  0.2 ~ "Hyper-induced in LECs",
      beta_interaction < -0.5 & LFC_Capillary < -0.2 ~ "Hyper-suppressed in LECs",
      TRUE ~ "Moderate / Unclassified Shift"
    )
  ) %>%
  arrange(padj_interaction)

# ------------------------------------------------------------------------------
# 4. Preview Top Lineage Resilience Candidates
# ------------------------------------------------------------------------------
head(
  res_combined %>% 
    select(gene, LFC_Capillary, beta_interaction, LFC_LEC, resilience_type, padj_interaction),
  20
)

write.csv(res_combined, file = "output/tables/deg_interaction_20261009.csv", row.names = TRUE)




