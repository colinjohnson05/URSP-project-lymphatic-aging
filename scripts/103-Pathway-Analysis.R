library(dplyr)
library(fgsea)
library(msigdbr)
library(data.table)

# 1. Load interaction results from CSV
res_combined <- read.csv("output/tables/deg_interaction_20261009.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(pval_interaction) & !is.na(beta_interaction))

# 2. Build ranked gene vector and break ties deterministically
set.seed(42)
ranked_genes <- res_combined %>%
  mutate(
    # Signed log-pvalue rank metric
    rank_metric = sign(beta_interaction) * -log10(pval_interaction) + rnorm(n(), 0, 1e-6)
  ) %>%
  filter(!is.infinite(rank_metric)) %>%
  arrange(desc(rank_metric))

gene_ranks <- setNames(ranked_genes$rank_metric, ranked_genes$gene)

# 3. Retrieve Gene Sets (Uppercase collection names matching msigdbr_collections)
# GO Biological Process (Collection "C5", Subcollection "GO:BP")
go_bp_pathways <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# MSigDB Hallmark (Collection "H")
hallmark_pathways <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# 4. Run fgsea on GO Biological Process
fgsea_res_go <- fgsea(
  pathways    = go_bp_pathways, 
  stats       = gene_ranks,
  scoreType   = "std",        # Standard two-tailed GSEA for signed ranks
  nPermSimple = 10000,        # Resolves unbalanced permutation warnings
  minSize     = 15,
  maxSize     = 500,
  nproc       = 1             # Prevents parallel process serialization warnings
) %>%
  arrange(padj)

# View top enriched GO Biological Processes
head(fgsea_res_go %>% select(pathway, NES, pval, padj), 15)

# Save GSEA table
fwrite(fgsea_res_go, "output/tables/gsea_go_bp_interaction.csv")
