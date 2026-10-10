library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

### Load endothelial aging atlas
atlas <- readRDS("data/HAECA_global_cleaned_corrected.RDS")

### View cell type statistics
meta <- atlas@meta.data

capillary_labels <- c("capillary", 
                     "capillary_arterial", 
                     "capillary_venous", 
                     "capillary_angiogenic")


meta$lineage_broad <- ifelse(
  meta$clustering_global == "lymphatic", "Lymphatic",
  ifelse(
    meta$clustering_global %in% capillary_labels, "Capillary", NA
  )
)

atlas@meta.data <- meta

donor_summary <- meta %>%
  group_by(Donor, Age_brackets, lineage_broad) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = lineage_broad, values_from = n_cells, values_fill = 0) %>%
  mutate(
    total_lec_cap = Capillary + Lymphatic,
    pct_lec  = (Lymphatic / total_lec_cap) * 100
  )

# Print donor summary table to console
print(donor_summary)

valid_donors <- donor_summary %>%
  filter(Lymphatic >= 15 & Capillary >= 15) %>%
  pull(Donor)

# Subset Seurat object to high-quality donors only
atlas_filtered <- subset(atlas, subset = Donor %in% valid_donors)
cells_keep <- colnames(atlas_filtered)[!is.na(atlas_filtered@meta.data$lineage_broad)]
atlas_filtered <- subset(atlas_filtered, cells = cells_keep)

# Save lung EC only
cells_keep <- colnames(atlas_filtered)[atlas_filtered@meta.data$Tissue == "Lung"]
atlas_filtered <- subset(atlas_filtered, cells = cells_keep)

cat(sprintf("Retained %d out of %d donors with sufficient cell counts for interaction modeling.\n", 
            length(valid_donors), nrow(donor_summary)))

saveRDS(atlas_filtered, file = "data/HAECA_filtered.RDS")

rm(atlas)
