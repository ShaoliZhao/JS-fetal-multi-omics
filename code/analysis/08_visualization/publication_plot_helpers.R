source("analysis/00_setup/project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(ComplexHeatmap)
})

js_palette <- c(
  Control = "#3565A8",
  JS = "#B5453C",
  Cerebellum = "#0D6B64",
  Kidney = "#8A5A20",
  Spatial = "#6B4E9B"
)

plot_grouped_fraction <- function(df, x, y, fill, filename) {
  p <- ggplot(df, aes({{ x }}, {{ y }}, fill = {{ fill }})) +
    geom_col(width = 0.8, color = "white", linewidth = 0.1) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_js_atlas()
  save_panel(p, filename, width = 5, height = 3.5)
}

plot_module_heatmap <- function(mat, annotation_col = NULL, filename = "module_heatmap.pdf") {
  pdf(file.path(paths$figures, filename), width = 6, height = 5)
  pheatmap(
    mat,
    annotation_col = annotation_col,
    scale = "row",
    border_color = NA,
    fontsize = 8
  )
  dev.off()
}

plot_dot_summary <- function(df, x, y, size, color, filename) {
  p <- ggplot(df, aes({{ x }}, {{ y }}, size = {{ size }}, color = {{ color }})) +
    geom_point(alpha = 0.85) +
    theme_js_atlas() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_panel(p, filename, width = 6, height = 4)
}
