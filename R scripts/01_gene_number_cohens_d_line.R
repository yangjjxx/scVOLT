# R 4.3.3
# Figure: Cohen's d sensitivity to the number of genes

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
  Cohensd = c(
    1.110281 ,	0.933517 	,0.820331 	,0.859420 	,1.601147 	,0.724015,
    3.958652 ,	2.947101 	,2.460195 	,1.772068 	,1.476798 	,2.275467 ,
    -0.588494, 	-0.133126, 	0.145138, 	0.461513, 	0.242651, 	0.704980, 
    0.365314 ,	-0.792304 ,	-0.479780, 	1.111813, 	2.088164, 	1.262805, 
    -0.411548, 	0.363245 ,	1.119096 ,	0.911774 ,	0.653325 ,	-0.068105 
  )
)



df$Drug <- factor(
  df$Drug,
  levels = c("PLX4720_A375", "Sorafenib", "PLX4720_451Lu", "Paclitaxel", "Docetaxel")
)

p <- ggplot(df, aes(x = Gene, y = Cohensd, color = Drug, group = Drug)) +
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
    limits = c(-4, 4),
    breaks = seq(-4, 4.0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "The number of genes",
    y = "Cohen's d",
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
ggsave("figures/gene_line_cohens_d.pdf", plot = p, width = 4, height = 2.5)
