# Public spatial plotting helpers used by the cell-bin analysis.
#
# This file is a compact public version of the original STimport.R. It keeps the
# colors and spatial plotting helpers used in the manuscript figures, while
# removing server-specific parallel settings and exploratory functions.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

my10cols <- c(
  "#EEC2E5", "#F9B27C", "#B6E2DC", "#ff9fe2", "#8FDBF3", "#CEBAF0",
  "#C6C3E1", "#A2C4F1", "#F9B29C", "#C6DCB9", "#F9C89B", "#FBDF9D",
  "#E9E4AF", "#FFCFD1", "#FBE3CD", "#87CEFA", "#cd9deb", "#a1efbb"
)

my25cols <- c(
  "#EEC2E5", "#F9B27C", "#B6E2DC", "#ff9fe2", "#8FDBF3", "#CEBAF0",
  "#C6C3E1", "#A2C4F1", "#F9B29C", "#C6DCB9", "#F9C89B", "#FBDF9D",
  "#E9E4AF", "#FFCFD1", "#FBE3CD", "#87CEFA", "#cd9deb", "#a1efbb",
  "#FF7F50", "#6B8E23", "#FFD700", "#40E0D0", "#DC143C", "#7FFFD4",
  "#BA55D3", "#708090"
)

color.module <- c(
  "#CCCCCC", "#ecb888", "#af88bb", "#a032cb", "#efbed6", "#fc496a",
  "#b6d37f", "#589336", "#7fd68e", "#52c465", "#3372e0", "#84d7f6",
  "#5394c3", "#6376b3", "#7f6cd7", "#c4ceff", "#fc9d40", "#5c95e0",
  "#cd7560", "#ff70e4", "#ff8738", "#ffcead", "#1cbf8b", "#b76d38",
  "#1584ff", "#7f006d", "#ffd35f", "#E66F73", "#F57F20", "#1DBB95",
  "#9CB79F", "#F0B8D2", "#A0485E", "#A0688E", "#C7E1DF", "#51B1DF",
  "#6D97D7", "#5D6193", "#CEC3E0", "#A9917E", "#7C7D80", "#F4E192",
  "#ADD666"
)

sDimplot <- function(object, color = NULL, group = "celltype", bgcolor = "#eeeeee",
                     pt.size = 0.2) {
  stopifnot(all(c("Spatial_1", "Spatial_2") %in% colnames(object@meta.data)))
  df <- data.frame(
    x = object$Spatial_1,
    y = object$Spatial_2,
    cluster = factor(object@meta.data[[group]])
  )
  if (is.null(color)) color <- hue_pal()(nlevels(df$cluster))

  ggplot(df, aes(x = x, y = y, color = cluster)) +
    geom_point(shape = 19, size = pt.size) +
    scale_color_manual(values = color) +
    coord_fixed() +
    labs(x = NULL, y = NULL, color = NULL) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = bgcolor),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )
}

smultiDimplot <- function(object, color = NULL, group = "celltype",
                          sample = "orig.ident", pt.size = 0.2,
                          ncol = 2, bgcolor = "white") {
  stopifnot(all(c("Spatial_1", "Spatial_2") %in% colnames(object@meta.data)))
  df <- data.frame(
    x = object$Spatial_1,
    y = object$Spatial_2,
    cluster = factor(object@meta.data[[group]]),
    sample = object@meta.data[[sample]]
  )
  if (is.null(color)) color <- hue_pal()(nlevels(df$cluster))

  plots <- lapply(split(df, df$sample), function(subdf) {
    ggplot(subdf, aes(x = x, y = y, color = cluster)) +
      geom_point(shape = 19, size = pt.size) +
      scale_color_manual(values = color, drop = FALSE) +
      coord_fixed() +
      labs(x = NULL, y = NULL, color = NULL, title = unique(subdf$sample)) +
      theme_bw() +
      theme(
        panel.background = element_rect(fill = bgcolor),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
  })
  wrap_plots(plots, ncol = ncol)
}

smultiFeaturePlot <- function(object, features, sample = "orig.ident",
                              assay = NULL, slot = "data", pt.size = 0.2,
                              ncol = 2, low.color = "#eeeeee",
                              high.color = "navy", same.scale = FALSE) {
  stopifnot(all(c("Spatial_1", "Spatial_2") %in% colnames(object@meta.data)))
  if (is.null(assay)) assay <- DefaultAssay(object)
  expr <- GetAssayData(object, assay = assay, slot = slot)
  features <- intersect(features, rownames(expr))

  plots <- lapply(features, function(gene) {
    df <- data.frame(
      x = object$Spatial_1,
      y = object$Spatial_2,
      expr = as.numeric(expr[gene, ]),
      sample = object@meta.data[[sample]]
    )
    if (same.scale) {
      scale_layer <- scale_color_gradient(limits = range(df$expr), low = low.color, high = high.color)
    } else {
      scale_layer <- scale_color_gradient(low = low.color, high = high.color)
    }
    subplots <- lapply(split(df, df$sample), function(subdf) {
      ggplot(subdf[order(subdf$expr), ], aes(x = x, y = y, color = expr)) +
        geom_point(shape = 19, size = pt.size) +
        scale_layer +
        coord_fixed() +
        labs(x = NULL, y = NULL, color = "Expr", title = paste(gene, unique(subdf$sample))) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
    })
    wrap_plots(subplots, ncol = ncol)
  })
  wrap_plots(plots, ncol = 1)
}
