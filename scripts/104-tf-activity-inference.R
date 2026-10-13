library(decoupleR)
library(dplyr)
library(tibble)

# 1. Load interaction results from CSV
res_combined <- read.csv("output/tables/deg_interaction_20261009.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(pval_interaction) & !is.na(beta_interaction))

# 2. Format signed statistic matrix for decoupleR
mat <- res_combined %>%
  mutate(stat = sign(beta_interaction) * -log10(pval_interaction)) %>%
  select(gene, stat) %>%
  filter(!is.na(stat) & !is.infinite(stat)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 3. Download CollecTRI directly from the OmniPath rescued mirror
collectri_url <- "https://rescued.omnipathdb.org/CollecTRI.csv"
dest_file <- "collectri_human.csv"

if (!file.exists(dest_file)) {
  download.file(collectri_url, destfile = dest_file, method = "libcurl")
}

# 4. Read local CollecTRI CSV and prepare network format
net <- read.csv(dest_file, stringsAsFactors = FALSE) %>%
  select(source, target, mor = weight) %>%
  distinct()

# 5. Infer TF activities using Univariate Linear Model (ULM)
tf_ulm <- run_ulm(
  mat = mat, 
  net = net, 
  .source = "source", 
  .target = "target", 
  .mor = "mor", 
  minsize = 5
)

# 6. Format and adjust p-values
tf_results <- tf_ulm %>%
  select(source, score, p_value) %>%
  rename(
    TF = source, 
    Activity_Score = score, 
    pval = p_value
  ) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  arrange(padj)

# 7. Save TF activity table to output
write.csv(tf_results, "output/tables/decoupler_tf_activities.csv", row.names = FALSE)

# Inspect top 15 inferred TFs
head(tf_results, 15)