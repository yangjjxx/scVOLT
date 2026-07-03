library(ggplot2)
library(scales)

# =========================
# 画 UMAP，颜色用百分位数
# =========================
p_score <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = drug_percentile)) +
  geom_point(size = 0.6, alpha = 0.75) +
  scale_color_gradientn(
    colors = c("#7DB7E8", "#BFDDF6", "#F9B56D", "#F47F6B", "#D94B3D"),
    limits = c(10, 90),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    oob = scales::squish
  ) +
  theme_classic() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = ""
  ) +
  theme(
    legend.title = element_text(size = 14, margin = margin(b = 20)),
    legend.text  = element_text(size = 12),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  ) +
  guides(
    color = guide_colorbar(
      ticks = TRUE,
      barheight = unit(3.2, "cm"),
      barwidth = unit(0.45, "cm")
    )
  )

print(p_score)

ggsave("figures/umap_by_drug_score.pdf", p_score, width = 5, height = 4)