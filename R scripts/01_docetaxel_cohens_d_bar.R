# R 4.3.3
# Figure: Docetaxel Cohen's d bar plot

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
  Drug = "Docetaxel",
  Method = factor(
    c("CaDRReS-Sc", "SCAD", "scIDUC", "Precily", "BeyondCell", "scVOLT"),
    levels = c("CaDRReS-Sc", "SCAD", "scIDUC", "Precily", "BeyondCell", "scVOLT")
  ),
  Value = c(-0.9133 ,2.2467 ,-0.0335 ,-0.1545 ,1.2978 ,1.3999 )
)

 

p <- ggplot(df, aes(x = Drug, y = Value, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "CaDRReS-Sc" = "#9EC4BE",
      "SCAD" = "#ABD0F1",
      "scIDUC" = "#DCE9F4",
      "Precily" = "#F19685",
      "BeyondCell" = "#FAC795",
      "scVOLT" = "#E56F5E"
    )
  ) +
  scale_y_continuous(
    limits = c(-2.0, 3.1),
    breaks = seq(-2.0, 3.0, 1.0),
    labels = function(x) sprintf("%.1f", x),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    expand = expansion(add = c(0.4, 0.3))
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.8, "cm"),
    legend.key.height = unit(0.15, "cm"),
    legend.spacing.x = unit(0.25, "cm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

print(p)
ggsave("figures/docetaxel_cohens_d.pdf", plot = p, width = 2.5, height = 2)
