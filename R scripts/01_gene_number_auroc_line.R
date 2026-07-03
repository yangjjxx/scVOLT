# R 4.3.3
# Figure: AUROC sensitivity to the number of genes

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
  Drug = rep(
    c("Docetaxel", "PLX4720_A375", "Paclitaxel", "Sorafenib", "PLX4720_451Lu"),
    each = 6
  ),
  Gene = rep(c(600, 800, 1000, 1200, 1400, 1600), times = 5),
  AUROC = c(
    # Docetaxel_data.h5ad
    0.787836, 0.763171, 0.744835, 0.736958, 0.867188, 0.690664,
    # A375_PLX-4720_data.h5ad
    0.995316, 0.974606, 0.966716, 0.919132, 0.885602, 0.965976,
    # Paclitaxel_data.h5ad
    0.314918, 0.445314, 0.531885, 0.626111, 0.567959, 0.695722,
    # Sorafenib_data.h5ad
    0.604167, 0.255833, 0.366389, 0.799167, 0.933611, 0.856111,
    # 451Lu_PLX-4720_data.h5ad
    0.379546, 0.621955, 0.790290, 0.766934, 0.726393, 0.464631
  )
)

df$Drug <- factor(
  df$Drug,
  levels = c("PLX4720_A375", "Sorafenib", "PLX4720_451Lu", "Paclitaxel", "Docetaxel")
)

p <- ggplot(df, aes(x = Gene, y = AUROC, color = Drug, group = Drug)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = c(
      "PLX4720_A375" = "#D55E00",
      "Sorafenib" = "#0072B2",
      "PLX4720_451Lu" = "#009E73",
      "Paclitaxel" = "#CC79A7",
      "Docetaxel" = "#E69F00"
    )
  ) +
  scale_x_continuous(
    limits = c(600, 1600),
    breaks = c(600, 800, 1000, 1200, 1400, 1600),
    expand = c(mult = 0.05, add = 0)
  ) +
  scale_y_continuous(
    limits = c(0.0, 1.02),
    breaks = seq(0.0, 1.0, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = "The number of genes",
    y = "AUROC",
    color = NULL
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0.2, "cm"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.35, "cm")
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE)
  )

print(p)
ggsave("figures/gene_line_auroc.pdf", plot = p, width = 4, height = 2.5)
