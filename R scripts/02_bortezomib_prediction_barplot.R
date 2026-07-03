library(ggplot2)
library(dplyr)

# =========================
# 定义 Resistant / Sensitive
# =========================
bar_df <- plot_df %>%
  mutate(
    time = as.character(time),
    time = gsub("t", "", time),
    time = factor(time, levels = c("0", "12", "48", "96")),
    prediction = ifelse(score_centered > 0, "Resistant", "Sensitive"),
    prediction = factor(prediction, levels = c("Sensitive", "Resistant"))
  ) %>%
  group_by(time, prediction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time) %>%
  mutate(
    percent = n / sum(n) * 100
  ) %>%
  ungroup()

print(bar_df)

# =========================
# 配色
# =========================
pred_colors <- c(
  "Resistant" = "#C23B17",
  "Sensitive" = "#8BB7DD"
)

# =========================
# 画堆叠柱状图
# =========================
p_bar <- ggplot(bar_df, aes(x = time, y = percent, fill = prediction)) +
  geom_col(
    width = 0.28,
    color = NA
  ) +
  scale_fill_manual(
    values = pred_colors,
    breaks = c("Resistant", "Sensitive")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Time point",
    y = "% of cells",
    fill = "Bortezomib\nscVolt prediction"
  ) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    legend.direction = "horizontal",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  )

print(p_bar)

ggsave("figures/bortezomib_prediction_barplot.pdf", p_bar, width = 5, height = 4.2)