# R 4.3.3
# Figure: AUROC heatmap across methods and drugs

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
    0.5562, 0.8632, 0.9292, 0.9066, 0.7393, 0.9970,
    # PLX4720_451Lu
    0.5179, 0.3625, 0.4862, 0.6374, 0.6077, 0.7669,
    # Sorafenib
    0.6356, 0.7569, 0.5244, 0.6533, 0.7490, 0.7992,
    # Docetaxel
    0.4588, 0.6506, 0.7816, 0.6411, 0.6131, 0.8731,
    # Paclitaxel
    0.2500, 0.9936, 0.5017, 0.5439, 0.8178, 0.8442
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
      color = Value > 0.75
    ),
    size = 5
  ) +
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#F7FBFF",
      "#DEEBF7",
      "#9ECAE1",
      "#4292C6",
      "#08519C"
    ),
    values = c(0, 0.45, 0.65, 0.82, 1),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    name = "AUROC"
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
      size = 15,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black"
    ),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title.align = 0.5,
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.margin = margin(5, 10, 5, 5)
  )

print(p)
ggsave("figures/auroc_heatmap.pdf", plot = p, width = 6, height = 4)
