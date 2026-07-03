library(ggplot2)

# 需要umap_df
p_time <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = time)) +
  geom_point(size = 0.4, alpha = 0.6) +
  coord_fixed() +
  theme_classic() +
  labs(x = "UMAP1", y = "UMAP2", color = "Time (h)") +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 3, alpha = 1)
    )
  )

print(p_time)

ggsave("figures/umap_by_time.pdf", p_time, width = 4.8, height = 4.5)
