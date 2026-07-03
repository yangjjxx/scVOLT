library(ggplot2)
library(dplyr)

# =========================
# 准备数据
# =========================
# 准备umap_df ，包含：
# time: 0 / 12 / 48 / 96 
# drug_score: 药效分数，越高越耐药

plot_df <- umap_df %>%
  mutate(
    time = as.character(time),
    time = gsub("t", "", time),
    time = factor(time, levels = c("0", "12", "48", "96"))
)

# =========================
# 以 t0 的中位数为 baseline
# =========================
baseline_t0 <- median(plot_df$drug_score[plot_df$time == "0"], na.rm = TRUE)

plot_df <- plot_df %>%
  mutate(
    score_centered = drug_score - baseline_t0
  )

# =========================
# 计算 25% 75%分位数
# =========================
boot_median_ci <- function(x, n_boot = 1000, conf = 0.95) {
  x <- x[!is.na(x)]
  boots <- replicate(n_boot, median(sample(x, replace = TRUE)))
  alpha <- (1 - conf) / 2
  quantile(boots, probs = c(alpha, 1 - alpha))
}

summary_df <- plot_df %>%
  group_by(time) %>%
  summarise(
    median_score = median(score_centered, na.rm = TRUE),
    ci_low = boot_median_ci(score_centered)[1],
    ci_high = boot_median_ci(score_centered)[2],
    n_cells = n(),
    .groups = "drop"
  )

summary_df <- summary_df %>%
  mutate(
    time_plot = as.numeric(time)  # factor: 0,12,48,96 -> 1,2,3,4
  )

print(summary_df)

# =========================
# 配色
# =========================
time_colors <- c(
  "0"  = "#F8766D",
  "12" = "#7CAE00",
  "48" = "#00BFC4",
  "96" = "#C77CFF"
)
summary_df <- plot_df %>%
  group_by(time) %>%
  summarise(
    median_score = median(score_centered, na.rm = TRUE),
    q25 = quantile(score_centered, 0.25, na.rm = TRUE),
    q75 = quantile(score_centered, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(time_plot = as.numeric(time))


# =========================
# 画图
# =========================
p <- ggplot(summary_df, aes(x = time_plot, y = median_score, color = time)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#D94B3D",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.15,
    color = "black",
    linewidth = 0.9
  ) +
  geom_point(size = 10) +
  scale_color_manual(values = time_colors) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("0", "12", "48", "96")
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Time (h)",
    y = "Centered resistance score"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p)

ggsave("figures/centered_median_q25_q75_drug_score.pdf", p, width = 4.2, height = 3.8)

