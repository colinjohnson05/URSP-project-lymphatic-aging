# This script loads data
# Load necessary packages
library(qs2)
library(Seurat)
library(anndata)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# Load datasets
lymph <- qs_read("data/Lymphatics_processed.qs2")

lung <- read_h5ad("data/lungMergedEndo10.h5ad")