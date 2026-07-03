# R 4.3.3
# Figure: Cohen's d sensitivity to latent dimension K

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
    0.966974 	,1.612407 	,0.859420 	,1.035440 	,0.720526 	,0.553748,
    1.270916 	,3.036222 	,1.772068 	,-0.041130 ,	0.894778 ,	0.242574,
    1.399910 	,0.511417 	,0.461513 	,0.490321 	,0.510176 	,0.233179,
    -0.959328 ,	0.093751 ,	1.111813 ,	-0.246022, 	0.275194, 	0.199020,
    -0.237548 ,	-0.138203, 	0.911774, 	1.091356, 	0.090012, 	0.584465
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
    limits = c(-2, 3.1),
    breaks = seq(-2, 3.1, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = expression("The latent dimension " * italic(K)),
    y = "Cohen's d",
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
ggsave("figures/K_line_cohens_d.pdf", plot = p, width = 4, height = 2.5, dpi = 600)
