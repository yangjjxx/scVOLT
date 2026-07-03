# R 4.3.3
# Figure: AUROC sensitivity to latent dimension K

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
  K = rep(c(10, 20, 30, 40, 50, 60), times = 5),
  AUROC = c(
    # Docetaxel
    0.739992, 0.873063, 0.736958, 0.782606, 0.692794, 0.651343,
    # PLX4720_A375
    0.932199, 0.997041, 0.919132, 0.508876, 0.744083, 0.550789,
    # Paclitaxel
    0.844234, 0.645744, 0.626111, 0.640195, 0.647942, 0.583141,
    # Sorafenib
    0.209444, 0.659444, 0.799167, 0.454444, 0.549722, 0.536944,
    # PLX4720_451Lu
    0.361361, 0.429096, 0.766934, 0.795629, 0.519686, 0.697531
  )
)

df$Drug <- factor(
  df$Drug,
  levels = c("PLX4720_A375", "Sorafenib", "PLX4720_451Lu", "Paclitaxel", "Docetaxel")
)

p <- ggplot(df, aes(x = K, y = AUROC, color = Drug, group = Drug)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
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
    limits = c(10, 60),
    breaks = c(10, 20, 30, 40, 50, 60),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = expression("The latent dimension " * italic(K)),
    y = "AUROC",
    color = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0.25, "cm"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.4, "cm")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE)
  )

print(p)
ggsave("figures/K_line_auroc.pdf", plot = p, width = 4, height = 2.5, dpi = 600)
