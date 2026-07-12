# Project configuration template for the fetal JS multi-omics atlas.

project_root <- Sys.getenv("JS_MULTIOMICS_ROOT", unset = ".")

paths <- list(
  data_processed = file.path(project_root, "data", "processed"),
  data_external = file.path(project_root, "data", "external_controlled"),
  results = file.path(project_root, "results"),
  figures = file.path(project_root, "figures"),
  metadata = file.path(project_root, "metadata")
)

dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$metadata, recursive = TRUE, showWarnings = FALSE)

theme_js_atlas <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
    )
}

save_panel <- function(plot, filename, width = 5, height = 4) {
  ggplot2::ggsave(
    filename = file.path(paths$figures, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
}
