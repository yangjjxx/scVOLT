# R 4.3.3
# Figure: AUPRC sensitivity to latent dimension K

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
    0.807595, 	0.883495, 	0.722756, 	0.809568, 	0.714890, 	0.620924, 
    0.910357, 	0.995561, 	0.903734, 	0.453577, 	0.639540, 	0.474297, 
    0.933591, 	0.843371, 	0.833779, 	0.834078, 	0.837849, 	0.782340, 
    0.350952, 	0.745810, 	0.773284, 	0.480385, 	0.575769, 	0.507785, 
    0.513501, 	0.516547, 	0.799114, 	0.782053, 	0.593100, 	0.658553
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
    y = "AUPRC",
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
ggsave("figures/K_line_auprc.pdf", plot = p, width = 4, height = 2.5, dpi = 600)
