library(dplyr)
library(ggplot2)

# 1. Load TF activity results from CSV
tf_results <- read.csv("output/tables/decoupler_tf_activities.csv", stringsAsFactors = FALSE)

# 2. Select top 20 TFs for plot
top_tfs <- tf_results %>%
  filter(!is.na(Activity_Score)) %>%
  slice_max(order_by = abs(Activity_Score), n = 20) %>%
  mutate(
    Direction = ifelse(Activity_Score > 0, "Preserved / Higher in LEC Aging", "Blunted / Lower in LEC Aging")
  )

# 3. Render Bar Plot
plot_tf_activity <- ggplot(top_tfs, aes(x = reorder(TF, Activity_Score), y = Activity_Score, fill = Direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Preserved / Higher in LEC Aging" = "#2b83ba",
    "Blunted / Lower in LEC Aging"   = "#d7191c"
  )) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
  theme_classic() +
  labs(
    title = "Inferred Master Transcription Factor Activities",
    subtitle = "CollectTRI Regulon Activity on Interaction Statistics (decoupleR ULM)",
    x = "Transcription Factor",
    y = "ULM Activity Score",
    fill = "Regulon Trajectory"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15, unit = "pt")
  )

print(plot_tf_activity)

# Export plot
ggsave("output/figures/decoupler_tf_activity_top20.svg", plot = plot_tf_activity, width = 7.0, height = 5.5, units = "in")
