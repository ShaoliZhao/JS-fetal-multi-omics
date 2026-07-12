# Monocle2 lineage reconstruction from a Seurat subset.
#
# Example:
# Rscript analysis/06_lineage_velocity/monocle_lineage_template.R \
#   --object results/kidney_atlas_public_workflow.rds \
#   --prefix kidney_npc_nephron \
#   --subset-col celltype \
#   --subset-values NPC,Podocyte,PT,LOH,LOH_DTL,DCT \
#   --group-col group \
#   --celltype-col celltype

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
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

choose_column <- function(obj, requested, candidates) {
  if (!is.null(requested) && requested %in% colnames(obj@meta.data)) return(requested)
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) > 0) return(hit[[1]])
  stop("Metadata column not found. Tried: ", paste(c(requested, candidates), collapse = ", "), call. = FALSE)
}

subset_for_lineage <- function(obj, subset_col, subset_values) {
  if (subset_col == "" || subset_values == "") return(obj)
  if (!subset_col %in% colnames(obj@meta.data)) stop("Subset column not found: ", subset_col, call. = FALSE)
  values <- trimws(strsplit(subset_values, ",")[[1]])
  keep <- colnames(obj)[obj@meta.data[[subset_col]] %in% values]
  if (length(keep) < 100) warning("Lineage subset has fewer than 100 cells.")
  subset(obj, cells = keep)
}

seurat_to_monocle2 <- function(obj) {
  counts <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
  pd <- new("AnnotatedDataFrame", data = obj@meta.data)
  fd <- new(
    "AnnotatedDataFrame",
    data = data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  )
  newCellDataSet(
    as(as.matrix(counts), "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size()
  )
}

choose_ordering_genes <- function(obj, cds, method = "variable", group_col = "group", max_genes = 2000) {
  if (method == "variable" && length(VariableFeatures(obj)) > 0) {
    return(head(VariableFeatures(obj), max_genes))
  }

  if (method == "dispersion") {
    disp <- dispersionTable(cds)
    genes <- disp %>%
      filter(mean_expression >= 0.1) %>%
      arrange(desc(dispersion_empirical)) %>%
      pull(gene_id)
    return(head(genes, max_genes))
  }

  if (method == "group_de" && group_col %in% colnames(pData(cds))) {
    formula <- as.formula(paste0("~", group_col))
    deg <- differentialGeneTest(cds, fullModelFormulaStr = deparse(formula), cores = 1)
    return(head(rownames(deg[order(deg$qval), ]), max_genes))
  }

  head(rownames(obj), max_genes)
}

run_monocle2 <- function(obj, prefix, group_col, celltype_col, ordering_method) {
  cds <- seurat_to_monocle2(obj)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  ordering_genes <- choose_ordering_genes(obj, cds, method = ordering_method, group_col = group_col)
  cds <- setOrderingFilter(cds, ordering_genes)
  write_csv(tibble(gene = ordering_genes), file.path(paths$results, paste0(prefix, "_ordering_genes.csv")))

  cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
  cds <- orderCells(cds)

  p_group <- plot_cell_trajectory(cds, color_by = group_col)
  p_celltype <- plot_cell_trajectory(cds, color_by = celltype_col)
  p_pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime")
  save_panel(p_group, paste0(prefix, "_trajectory_by_group.pdf"), width = 5, height = 4)
  save_panel(p_celltype, paste0(prefix, "_trajectory_by_celltype.pdf"), width = 5.5, height = 4)
  save_panel(p_pseudotime, paste0(prefix, "_trajectory_by_pseudotime.pdf"), width = 5, height = 4)

  meta <- pData(cds) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id")
  write_csv(meta, file.path(paths$results, paste0(prefix, "_monocle_cell_metadata.csv")))
  saveRDS(cds, file.path(paths$results, paste0(prefix, "_monocle2_cds.rds")))
  cds
}

main <- function() {
  opts <- parse_args(list(
    object = "results/kidney_atlas_public_workflow.rds",
    prefix = "lineage",
    subset_col = "",
    subset_values = "",
    group_col = "group",
    celltype_col = "celltype",
    ordering_method = "variable"
  ))

  object_path <- resolve_path(opts$object)
  if (!file.exists(object_path)) stop("Seurat object not found: ", object_path, call. = FALSE)
  obj <- readRDS(object_path)
  obj <- subset_for_lineage(obj, opts$subset_col, opts$subset_values)
  group_col <- choose_column(obj, opts$group_col, c("group", "disease_status"))
  celltype_col <- choose_column(obj, opts$celltype_col, c("celltype", "cell_type", "ct", "seurat_clusters"))
  run_monocle2(obj, opts$prefix, group_col, celltype_col, opts$ordering_method)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
