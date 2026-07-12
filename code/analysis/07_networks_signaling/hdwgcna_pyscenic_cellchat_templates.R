# Network and signaling workflows: CellChat, pySCENIC export and hdWGCNA.
#
# Examples:
# Rscript analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R \
#   --mode cellchat \
#   --object results/kidney_atlas_public_workflow.rds \
#   --prefix kidney \
#   --celltype-col celltype \
#   --group-col group
#
# Rscript analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R \
#   --mode pyscenic_export \
#   --object results/kidney_atlas_public_workflow.rds \
#   --prefix kidney

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
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

subset_object <- function(obj, subset_col = "", subset_values = "") {
  if (subset_col == "" || subset_values == "") return(obj)
  if (!subset_col %in% colnames(obj@meta.data)) stop("Subset column not found: ", subset_col, call. = FALSE)
  values <- trimws(strsplit(subset_values, ",")[[1]])
  subset(obj, cells = colnames(obj)[obj@meta.data[[subset_col]] %in% values])
}

run_cellchat_comparison <- function(obj, prefix, group_col, celltype_col, min_cells = 10) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat is required for --mode cellchat.", call. = FALSE)
  }
  library(CellChat)

  object_list <- SplitObject(obj, split.by = group_col)
  chats <- lapply(object_list, function(x) {
    chat <- createCellChat(x, group.by = celltype_col)
    chat@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
    chat <- subsetData(chat)
    chat <- identifyOverExpressedGenes(chat)
    chat <- identifyOverExpressedInteractions(chat)
    chat <- computeCommunProb(chat)
    chat <- filterCommunication(chat, min.cells = min_cells)
    chat <- computeCommunProbPathway(chat)
    chat <- aggregateNet(chat)
    netAnalysis_computeCentrality(chat)
  })

  merged <- mergeCellChat(chats, add.names = names(chats))
  saveRDS(merged, file.path(paths$results, paste0(prefix, "_cellchat_merged.rds")))

  interactions <- subsetCommunication(merged)
  write_csv(interactions, file.path(paths$results, paste0(prefix, "_cellchat_interactions.csv")))

  p_count <- compareInteractions(merged, show.legend = FALSE, group = seq_along(chats), measure = "count")
  p_weight <- compareInteractions(merged, show.legend = FALSE, group = seq_along(chats), measure = "weight")
  save_panel(p_count + p_weight, paste0(prefix, "_cellchat_interaction_summary.pdf"), width = 6, height = 3)
  merged
}

export_pyscenic_inputs <- function(obj, prefix) {
  out_dir <- file.path(paths$results, paste0(prefix, "_pyscenic_input"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  counts <- GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
  Matrix::writeMM(counts, file.path(out_dir, "matrix.mtx"))
  write_lines(rownames(counts), file.path(out_dir, "genes.tsv"))
  write_lines(colnames(counts), file.path(out_dir, "cells.tsv"))

  meta <- obj@meta.data %>%
    tibble::rownames_to_column("cell_id")
  write_csv(meta, file.path(out_dir, "cell_metadata.csv"))

  command_template <- c(
    "# Fill these paths with local pySCENIC database files before running.",
    "pyscenic grn matrix.loom hs_hgnc_tfs.txt -o adjacencies.tsv --num_workers 20",
    "pyscenic ctx adjacencies.tsv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname matrix.loom --mode dask_multiprocessing --output regulons.csv --num_workers 20",
    "pyscenic aucell matrix.loom regulons.csv --output auc_mtx.loom --num_workers 20"
  )
  write_lines(command_template, file.path(out_dir, "pyscenic_commands_template.sh"))
  out_dir
}

run_hdwgcna_if_available <- function(obj, prefix, group_col, celltype_col, subset_name = "network") {
  if (!requireNamespace("hdWGCNA", quietly = TRUE)) {
    message("hdWGCNA is not installed; writing metacell input summary only.")
    meta <- obj@meta.data %>%
      count(.data[[group_col]], .data[[celltype_col]], name = "n_cells")
    write_csv(meta, file.path(paths$results, paste0(prefix, "_hdwgcna_group_celltype_counts.csv")))
    return(invisible(NULL))
  }
  library(hdWGCNA)

  obj <- SetupForWGCNA(
    obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = subset_name
  )
  obj <- MetacellsByGroups(
    obj,
    group.by = c(celltype_col, group_col),
    reduction = "pca",
    k = 25,
    max_shared = 10,
    ident.group = celltype_col
  )
  obj <- NormalizeMetacells(obj)
  obj <- SetDatExpr(obj, group_name = unique(obj@meta.data[[celltype_col]])[1], group.by = celltype_col)
  obj <- TestSoftPowers(obj, networkType = "signed")
  saveRDS(obj, file.path(paths$results, paste0(prefix, "_hdwgcna_setup.rds")))
  obj
}

main <- function() {
  opts <- parse_args(list(
    mode = "cellchat",
    object = "results/kidney_atlas_public_workflow.rds",
    prefix = "network",
    group_col = "group",
    celltype_col = "celltype",
    subset_col = "",
    subset_values = "",
    min_cells = "10"
  ))

  object_path <- resolve_path(opts$object)
  if (!file.exists(object_path)) stop("Seurat object not found: ", object_path, call. = FALSE)
  obj <- readRDS(object_path)
  obj <- subset_object(obj, opts$subset_col, opts$subset_values)
  group_col <- choose_column(obj, opts$group_col, c("group", "disease_status"))
  celltype_col <- choose_column(obj, opts$celltype_col, c("celltype", "cell_type", "ct", "seurat_clusters"))

  if (opts$mode == "cellchat") {
    run_cellchat_comparison(obj, opts$prefix, group_col, celltype_col, min_cells = as.integer(opts$min_cells))
  } else if (opts$mode == "pyscenic_export") {
    export_pyscenic_inputs(obj, opts$prefix)
  } else if (opts$mode == "hdwgcna") {
    run_hdwgcna_if_available(obj, opts$prefix, group_col, celltype_col)
  } else {
    stop("Unknown --mode: ", opts$mode, ". Use cellchat, pyscenic_export or hdwgcna.", call. = FALSE)
  }
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
