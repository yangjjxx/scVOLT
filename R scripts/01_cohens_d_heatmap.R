# R 4.3.3
# Figure: Cohen's d heatmap across methods and drugs

required_packages <- c("ggplot2", "grid")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it first.", pkg), call. = FALSE)
  }
}

library(ggplot2)
library(grid)

dir.create("figures", showWarnings = FALSE)

df <- data.frame(
  Drug = rep(c("PLX4720_A375", "PLX4720_451Lu", "Sorafenib", "Docetaxel", "Paclitaxel"), each = 6),
  Method = rep(c("CaDRReS-Sc", "SCAD", "scIDUC", "Precily", "BeyondCell", "scVOLT"), times = 5),
  Value = c(
    # PLX4720_A375
    0.2119, 0.9494, 2.0546, -1.8252, 0.8998, 3.0362,
    # PLX4720_451Lu
    0.2563, 0.0111, 0.3214, -0.5298, 0.3816, 0.9118,
    # Sorafenib
    0.5442, 0.5365, 0.1041, 0.5640, 0.9147, 2.0882,
    # Docetaxel
    -0.1201, 0.6682, 1.0250, 0.4918, 0.3923, 1.6124,
    # Paclitaxel
    -0.9133, 2.2467, -0.0335, -0.1545, 1.2978, 1.3999
  )
)

# Drug order from top to bottom in the heatmap.
df$Drug <- factor(
  df$Drug,
  levels = rev(c("PLX4720_A375", "Sorafenib", "PLX4720_451Lu", "Paclitaxel", "Docetaxel"))
)

# Method order from left to right.
df$Method <- factor(
  df$Method,
  levels = c("CaDRReS-Sc", "SCAD", "scIDUC", "Precily", "BeyondCell", "scVOLT")
)

p <- ggplot(df, aes(x = Method, y = Drug, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(
      label = sprintf("%.2f", Value),
      color = abs(Value) > 1.5
    ),
    size = 6
  ) +
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-3.2, 3.2),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    name = "Cohen's d"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 18,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.title.align = 0.5,
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.margin = margin(5, 10, 5, 5)
  )

print(p)
ggsave("figures/cohen_heatmap.pdf", plot = p, width = 7, height = 5)
