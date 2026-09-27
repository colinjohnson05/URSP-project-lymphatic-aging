# Run DEVIL pipeline for DEG analysis
library(devil)
library(scRNAseq)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load dataset
sce <- qs_read("data/lymphatics_for_devil.qs2")
sce

### initialize DEVIL model
set.seed(123)
meta_df <- sce[[2]]

meta_df$age_group <- factor(meta_df$age_group, levels = c("Young", "Old"))
meta_df$donor_unique <- factor(meta_df$donor_unique)

# Set up model design
design <- model.matrix(
  ~ 0 + age_group,
  data = meta_df
)

head(design)
colnames(design)

# Group data by patients
Y <- as.matrix(sce[[1]])
grouped_data = devil::group_data(Y, design, meta_df$donor_unique)

# Fit model
devil_fit <- fit_devil(
  x = grouped_data$input_matrix,
  design_matrix = grouped_data$design_matrix,
  clusters = grouped_data$clusters,
  verbose = TRUE,
  init_beta_rough = TRUE,
  size_factors = "normed_sum",
  overdispersion = "MOM"
)

# Test DEGs
make_contrast <- function(design, from, to) {
  stopifnot(from %in% colnames(design), to %in% colnames(design))
  c <- rep(0, ncol(design))
  names(c) <- colnames(design)
  c[from] <- 1
  c[to] <- -1
  as.numeric(c)
}

# Make contrast for DEG testing, +1 = up in old
contrast <- make_contrast(design, from = "age_groupOld", to = "age_groupYoung")
contrast

test <- devil::test_de(
  devil_fit,
  contrast = contrast,
  max_lfc  = 10 # Cap extreme fold changes
)

test %>%
  dplyr::arrange(adj_pval, desc(abs(lfc))) %>%
  dplyr::select(name, lfc, pval, adj_pval) %>%
  head(20)

# Volcano plot
devil::plot_volcano(
  test,
  lfc_cut = 1,
  pval_cut = 0.05,
  labels = TRUE,
  point_size = 1.8,
  title = "young vs aged"
)