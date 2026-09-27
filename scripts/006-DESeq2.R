# Pseudobulk DESeq2 Analysis Pipeline for Aging Lymphatics
library(qs2)
library(DESeq2)
library(Matrix)
library(dplyr)

set.seed(123)

# 1. Load dataset (list containing counts and metadata)
sce <- qs_read("data/lymphatics_for_devil.qs2")

# Extract counts and metadata matching your DEVIL object structure
counts  <- if (!is.null(sce$counts)) sce$counts else sce[[1]]
meta_df <- if (!is.null(sce$meta)) sce$meta else sce[[2]]

# Set factor levels (Young as baseline reference)
meta_df$age_group    <- factor(meta_df$age_group, levels = c("Young", "Old"))
meta_df$donor_unique <- factor(meta_df$donor_unique)

# 2. Construct Pseudobulk Samples by Donor
# Create a unique sample identifier for each donor
meta_df$sample_id <- paste(meta_df$donor_unique, meta_df$age_group, sep = "___")

# Filter out donors with fewer than 10 total cells to reduce low-coverage noise
cell_counts   <- table(meta_df$sample_id)
valid_samples <- names(cell_counts[cell_counts >= 10])

cells_to_keep <- meta_df$sample_id %in% valid_samples
counts_sub    <- counts[, cells_to_keep, drop = FALSE]
meta_sub      <- meta_df[cells_to_keep, ]

# 3. Aggregate Single-Cell Counts to Pseudobulk Counts (Donor Level)
sample_factor <- factor(meta_sub$sample_id)
group_mat     <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
colnames(group_mat) <- levels(sample_factor)

# Fast sparse matrix aggregation (Genes x Donors)
pb_counts <- as.matrix(counts_sub %*% group_mat)
colnames(pb_counts) <- levels(sample_factor)

# Build sample-level metadata
pb_meta <- unique(meta_sub[, c("sample_id", "donor_unique", "age_group")])
rownames(pb_meta) <- pb_meta$sample_id
pb_meta <- pb_meta[colnames(pb_counts), ]

# 4. Run DESeq2 Pipeline
# Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = pb_meta,
  design    = ~ age_group
)

# Pre-filter low-expression genes (>= 10 counts in at least 3 donor samples)
keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_genes, ]

# Fit DESeq2 model
dds <- DESeq(dds)

# Extract results: Old vs Young (Positive log2FC = Upregulated in Old)
res <- results(dds, contrast = c("age_group", "Old", "Young"))

# 5. Format and Save Results
deg_deseq2 <- as.data.frame(res) %>%
  mutate(gene = rownames(.)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  select(gene, log2FoldChange, pvalue, padj, everything())

# Save CSV and binary outputs
write.csv(deg_deseq2, "data/lymphatic_aging_degs_deseq2.csv", row.names = FALSE)
qs_save(deg_deseq2, "data/lymphatic_aging_degs_deseq2.qs2")

cat("\n--- DESeq2 Analysis Complete ---\n")
cat("Total pseudobulk samples (donors):", nrow(pb_meta), "\n")
cat("Significant DEGs (FDR < 0.05):", sum(deg_deseq2$padj < 0.05), "\n")
head(deg_deseq2, 10)