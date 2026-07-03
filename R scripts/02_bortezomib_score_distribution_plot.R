library(ggplot2)
library(dplyr)
library(tidyr)

# =========================
# 保证顺序
# =========================
plot_df$time <- factor(plot_df$time, levels = c("0", "12", "48", "96"))
plot_df$state <- factor(plot_df$state, levels = c("sensitive-like", "resistant-like"))

# =========================
# 设置 bin
# =========================
n_bins <- 35

x_min <- min(plot_df$score_centered, na.rm = TRUE)
x_max <- max(plot_df$score_centered, na.rm = TRUE)

binwidth <- (x_max - x_min) / n_bins

breaks <- seq(x_min, x_max + binwidth, by = binwidth)

# =========================
# 手动统计每个 time × state × bin 的细胞比例
# =========================
plot_df_bin <- plot_df %>%
  mutate(
    bin = cut(
      score_centered,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  )

bin_df <- plot_df_bin %>%
  group_by(time, state, bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(time, state, bin, fill = list(n = 0)) %>%
  group_by(time) %>%
  mutate(
    n_time = sum(n),
    frac = n / n_time
  ) %>%
  ungroup()
# =========================
# 给每个 bin 计算中心位置
# =========================
bin_levels <- levels(plot_df_bin$bin)

bin_info <- data.frame(
  bin = bin_levels,
  bin_id = seq_along(bin_levels),
  x = breaks[-length(breaks)] + binwidth / 2
)

bin_df <- bin_df %>%
  left_join(bin_info, by = "bin")

# =========================
# 画图 bar+折线
# =========================
p <- ggplot(bin_df, aes(x = x, y = frac, fill = state, color = state)) +
  geom_col(
    width = binwidth,
    alpha = 0.55,
    position = "identity"
  ) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  facet_grid(time ~ ., scales = "fixed") +
  coord_cartesian(
    xlim = c(-0.01, 0.010),
    ylim = c(0, 0.24)
  ) +
  scale_fill_manual(values = c(
    "sensitive-like" = "#7DB7E8",
    "resistant-like" = "#F47F6B"
  )) +
  scale_color_manual(values = c(
    "sensitive-like" = "#4C78A8",
    "resistant-like" = "#D94B3D"
  )) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Fraction"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text.y = element_text(size = 16),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    panel.spacing.y = unit(0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p)

ggsave("figures/bortezomib_score_distribution.pdf", p, width = 4.5, height = 7.2)