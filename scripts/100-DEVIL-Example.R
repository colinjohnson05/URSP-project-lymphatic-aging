library(devil)
library(scRNAseq)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(dplyr)

sce <- scRNAseq::BaronPancreasData()
sce

counts <- SummarizedExperiment::assay(sce, "counts")
meta <- as.data.frame(SummarizedExperiment::colData(sce))

cat("Genes:", nrow(counts), "\nCells:", ncol(counts), "\n")
stopifnot("label" %in% colnames(meta))
head(meta$label)

# keep 3 largest cell types
top3 <- names(sort(table(meta$label), decreasing = TRUE))[1:3]
keep_cells <- meta$label %in% top3
counts <- counts[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]

# gene filter: expressed (>=1 UMI) in >= 1% of kept cells
min_cells <- max(1, floor(0.01 * ncol(counts)))
keep_genes <- Matrix::rowSums(counts >= 1) >= min_cells
counts <- counts[keep_genes, , drop = FALSE]

cat("After filtering — Genes:", nrow(counts), "Cells:", ncol(counts), "\n")
table(meta$label)

# demo mode: top 500 genes by total counts
if (nrow(counts) > 500) {
  ord <- order(Matrix::rowSums(counts), decreasing = TRUE)
  counts <- counts[ord[seq_len(500)], ]
}

meta$label <- droplevels(factor(meta$label))
design <- model.matrix(~ 0 + label, data = meta)
colnames(design) <- gsub("^label", "", colnames(design))
colnames(design)

clusters <- NULL
if ("donor" %in% names(meta)) clusters <- factor(meta$donor)
if (is.null(clusters) && "patient" %in% names(meta)) clusters <- factor(meta$patient)

fit <- devil::fit_devil(
  x = as.matrix(counts),
  design_matrix = design,
  clusters = as.factor(clusters),
  overdispersion = "MOM",
  offset = 1e-6,
  init_overdispersion = NULL,
  size_factors = "normed_sum",
  
  verbose = TRUE,
  max_iter = 200,
  tolerance = 1e-4
)

make_contrast <- function(design, from, to) {
  stopifnot(from %in% colnames(design), to %in% colnames(design))
  c <- rep(0, ncol(design))
  names(c) <- colnames(design)
  c[from] <- 1
  c[to] <- -1
  as.numeric(c)
}

contrast <- make_contrast(design, from = "beta", to = "ductal")
contrast

test <- devil::test_de(
  fit,
  contrast = contrast,
  max_lfc  = 10 # Cap extreme fold changes
)

# Add gene names if missing
if (!("name" %in% colnames(test))) {
  if (!is.null(rownames(counts))) {
    test$name <- rownames(counts)
  } else {
    test$name <- as.character(seq_len(nrow(test)))
  }
}

test %>%
  dplyr::arrange(adj_pval, desc(abs(lfc))) %>%
  dplyr::select(name, lfc, pval, adj_pval) %>%
  head(10)

devil::plot_volcano(
  test,
  lfc_cut = 1,
  pval_cut = 0.05,
  labels = TRUE,
  point_size = 1.8,
  title = "beta vs ductal"
)

