# Stereo-seq/cell-bin label transfer and spatial plotting workflow.
#
# Example:
# Rscript analysis/04_spatial_transcriptomics/stereo_label_transfer_spatial_plots.R \
#   --reference results/cerebellum_atlas_public_workflow.rds \
#   --spatial results/stereo_cellbin_object.rds \
#   --prefix ofd1_cellbin \
#   --celltype-col celltype \
#   --x-col spatial_x \
#   --y-col spatial_y \
#   --genes NRN1,BARHL1,S100B,CALB1

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  keys <- grep("^--", args)
  for (i in keys) {
    key <- sub("^--", "", args[[i]])
    value <- if (i + 1 <= length(args) && !grepl("^--", args[[i + 1]])) args[[i + 1]] else TRUE
    out[[key]] <- value
  }
  out
}

resolve_path <- function(path) {
  if (is.null(path) || is.na(path) || path == "") return(path)
  if (grepl("^/", path)) return(path)
  file.path(project_root, path)
}

require_file <- function(path, label) {
  path <- resolve_path(path)
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
  path
}

choose_column <- function(obj, requested, candidates) {
  if (!is.null(requested) && requested %in% colnames(obj@meta.data)) return(requested)
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) > 0) return(hit[[1]])
  stop("Could not find metadata column. Tried: ", paste(c(requested, candidates), collapse = ", "), call. = FALSE)
}

detect_coordinate_columns <- function(obj, x_col, y_col) {
  meta <- colnames(obj@meta.data)
  x <- choose_column(obj, x_col, c("spatial_x", "x", "coord_x", "imagecol", "pxl_col_in_fullres"))
  y <- choose_column(obj, y_col, c("spatial_y", "y", "coord_y", "imagerow", "pxl_row_in_fullres"))
  if (!(x %in% meta && y %in% meta)) stop("Spatial coordinate columns not found.", call. = FALSE)
  c(x = x, y = y)
}

run_label_transfer <- function(reference, spatial, celltype_col, dims = 1:30) {
  DefaultAssay(reference) <- "RNA"
  DefaultAssay(spatial) <- "RNA"
  anchors <- FindTransferAnchors(
    reference = reference,
    query = spatial,
    normalization.method = "LogNormalize",
    dims = dims
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference@meta.data[[celltype_col]],
    dims = dims
  )
  AddMetaData(spatial, metadata = predictions)
}

plot_spatial_feature <- function(obj, feature, coord_cols, prefix, point_size = 0.08) {
  if (!feature %in% rownames(obj)) {
    warning("Skipping missing feature: ", feature)
    return(invisible(NULL))
  }
  df <- obj@meta.data
  df$expr <- FetchData(obj, vars = feature)[, 1]
  p <- ggplot(df, aes(.data[[coord_cols["x"]]], .data[[coord_cols["y"]]], color = expr)) +
    geom_point(size = point_size, alpha = 0.95) +
    scale_color_viridis_c(option = "magma") +
    coord_equal() +
    theme_void() +
    labs(color = feature)
  save_panel(p, paste0(prefix, "_feature_", feature, ".pdf"), width = 5, height = 5)
}

plot_spatial_category <- function(obj, category_col, coord_cols, prefix, point_size = 0.08) {
  df <- obj@meta.data
  p <- ggplot(df, aes(.data[[coord_cols["x"]]], .data[[coord_cols["y"]]], color = .data[[category_col]])) +
    geom_point(size = point_size, alpha = 0.95) +
    coord_equal() +
    theme_void() +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(color = category_col)
  save_panel(p, paste0(prefix, "_", category_col, ".pdf"), width = 5.6, height = 5)
}

add_program_score <- function(obj, gene_string, score_name) {
  genes <- intersect(strsplit(gene_string, ",")[[1]], rownames(obj))
  if (length(genes) < 2) {
    warning("Program ", score_name, " has fewer than 2 genes present; skipping.")
    return(obj)
  }
  obj <- AddModuleScore(obj, features = list(genes), name = paste0(score_name, "_"))
  score_col <- paste0(score_name, "_1")
  colnames(obj@meta.data)[colnames(obj@meta.data) == score_col] <- score_name
  obj
}

export_spatial_metadata <- function(obj, prefix) {
  meta <- obj@meta.data
  meta$cellbin_id <- rownames(meta)
  write_csv(meta, file.path(paths$results, paste0(prefix, "_spatial_metadata_with_predictions.csv")))
}

main <- function() {
  opts <- parse_args(list(
    reference = "results/cerebellum_atlas_public_workflow.rds",
    spatial = "results/stereo_cellbin_object.rds",
    prefix = "stereo_cellbin",
    celltype_col = "celltype",
    x_col = "",
    y_col = "",
    genes = "NRN1,BARHL1,S100B,CALB1",
    scaffold_genes = "HES1,SOX2,VIM,FABP7,HOPX",
    gc_spec_genes = "ATOH1,EOMES,NEUROD1,NHLH1,NRN1",
    gc_wire_genes = "BARHL1,CALB1,GRID2,PCP4"
  ))

  reference <- readRDS(require_file(opts$reference, "Reference Seurat object"))
  spatial <- readRDS(require_file(opts$spatial, "Spatial Seurat object"))
  celltype_col <- choose_column(reference, opts$celltype_col, c("celltype", "cell_type", "ct", "seurat_clusters"))
  coord_cols <- detect_coordinate_columns(spatial, opts$x_col, opts$y_col)

  spatial <- run_label_transfer(reference, spatial, celltype_col)
  spatial$predicted_cell_type <- spatial$predicted.id

  spatial <- add_program_score(spatial, opts$scaffold_genes, "scaffold_score")
  spatial <- add_program_score(spatial, opts$gc_spec_genes, "gc_spec_score")
  spatial <- add_program_score(spatial, opts$gc_wire_genes, "gc_wire_score")

  export_spatial_metadata(spatial, opts$prefix)
  plot_spatial_category(spatial, "predicted_cell_type", coord_cols, opts$prefix)

  genes <- trimws(strsplit(opts$genes, ",")[[1]])
  for (gene in genes) plot_spatial_feature(spatial, gene, coord_cols, opts$prefix)
  for (score in intersect(c("scaffold_score", "gc_spec_score", "gc_wire_score"), colnames(spatial@meta.data))) {
    df <- spatial@meta.data
    df$expr <- df[[score]]
    p <- ggplot(df, aes(.data[[coord_cols["x"]]], .data[[coord_cols["y"]]], color = expr)) +
      geom_point(size = 0.08, alpha = 0.95) +
      scale_color_viridis_c(option = "magma") +
      coord_equal() +
      theme_void() +
      labs(color = score)
    save_panel(p, paste0(opts$prefix, "_program_", score, ".pdf"), width = 5, height = 5)
  }

  saveRDS(spatial, file.path(paths$results, paste0(opts$prefix, "_with_label_transfer.rds")))
  invisible(spatial)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
