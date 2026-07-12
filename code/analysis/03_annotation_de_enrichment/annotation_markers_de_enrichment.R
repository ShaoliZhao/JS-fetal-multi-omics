# Cell-type annotation, marker detection, disease DE and GO enrichment.
#
# Example:
# Rscript analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R \
#   --object results/cerebellum_atlas_public_workflow.rds \
#   --prefix cerebellum \
#   --celltype-col celltype \
#   --group-col group \
#   --ident1 JS \
#   --ident2 ctrl

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  if (length(args) == 0) return(out)
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
  stop("None of these metadata columns were found: ", paste(c(requested, candidates), collapse = ", "), call. = FALSE)
}

marker_panels <- list(
  cerebellum = c("SOX2", "HES1", "ATOH1", "EOMES", "NRN1", "CALB1", "PAX2", "SOX14", "AQP4", "FOXJ1", "CLDN5", "MBP"),
  kidney = c("SIX2", "NPHS1", "NPHS2", "LRP2", "CUBN", "UMOD", "SLC12A3", "GATA3", "AQP2", "COL1A1", "PECAM1", "C1QA")
)

plot_marker_dotplot <- function(obj, genes, group_col, prefix) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(invisible(NULL))
  p <- DotPlot(obj, features = genes, group.by = group_col) +
    RotatedAxis() +
    theme_js_atlas() +
    ggtitle(paste(prefix, "marker genes"))
  save_panel(p, paste0(prefix, "_marker_dotplot.pdf"), width = 9, height = 4)
}

run_markers <- function(obj, celltype_col, prefix) {
  Idents(obj) <- celltype_col
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  write_csv(markers, file.path(paths$results, paste0(prefix, "_celltype_markers.csv")))
  markers
}

run_celltype_de <- function(obj, celltype_col, group_col, ident1, ident2, prefix) {
  celltypes <- sort(unique(as.character(obj@meta.data[[celltype_col]])))
  de_dir <- file.path(paths$results, paste0(prefix, "_celltype_de"))
  dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

  de_tables <- lapply(celltypes, function(ct) {
    sub_obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[celltype_col]] == ct])
    if (ncol(sub_obj) < 20 || length(unique(sub_obj@meta.data[[group_col]])) < 2) return(NULL)
    Idents(sub_obj) <- group_col
    de <- FindMarkers(
      sub_obj,
      ident.1 = ident1,
      ident.2 = ident2,
      min.pct = 0.1,
      logfc.threshold = 0.25
    )
    de$gene <- rownames(de)
    de$celltype <- ct
    safe_ct <- gsub("[^A-Za-z0-9_]+", "_", ct)
    write_csv(de, file.path(de_dir, paste0(prefix, "_", safe_ct, "_", ident1, "_vs_", ident2, ".csv")))
    de
  })

  bind_rows(de_tables)
}

run_go <- function(de_table, prefix, lfc_col = "avg_log2FC") {
  if (nrow(de_table) == 0) return(invisible(NULL))
  go_dir <- file.path(paths$results, paste0(prefix, "_go"))
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

  split_table <- de_table %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05, abs(.data[[lfc_col]]) > 0.25) %>%
    mutate(direction = ifelse(.data[[lfc_col]] > 0, "up", "down")) %>%
    group_by(celltype, direction) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")

  for (i in seq_len(nrow(split_table))) {
    genes <- split_table$genes[[i]]
    if (length(genes) < 5) next
    entrez <- suppressMessages(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    if (nrow(entrez) < 5) next
    ego <- enrichGO(
      gene = unique(entrez$ENTREZID),
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      readable = TRUE
    )
    safe_ct <- gsub("[^A-Za-z0-9_]+", "_", split_table$celltype[[i]])
    out <- file.path(go_dir, paste0(prefix, "_", safe_ct, "_", split_table$direction[[i]], "_GO_BP.csv"))
    write_csv(as.data.frame(ego), out)
  }
}

plot_composition_from_object <- function(obj, celltype_col, group_col, prefix) {
  comp <- obj@meta.data %>%
    count(.data[[group_col]], .data[[celltype_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup()
  colnames(comp)[1:2] <- c("group", "celltype")
  write_csv(comp, file.path(paths$results, paste0(prefix, "_celltype_composition.csv")))
  p <- ggplot(comp, aes(celltype, percent, fill = group)) +
    geom_col(position = "dodge", width = 0.75) +
    coord_flip() +
    theme_js_atlas() +
    labs(x = NULL, y = "Percent of nuclei")
  save_panel(p, paste0(prefix, "_celltype_composition.pdf"), width = 5, height = 5)
}

main <- function() {
  opts <- parse_args(list(
    object = "results/cerebellum_atlas_public_workflow.rds",
    prefix = "cerebellum",
    celltype_col = "",
    group_col = "group",
    ident1 = "JS",
    ident2 = "ctrl"
  ))

  obj <- readRDS(require_file(opts$object, "Seurat object"))
  celltype_col <- choose_column(obj, opts$celltype_col, c("celltype", "cell_type", "ct", "seurat_clusters"))
  group_col <- choose_column(obj, opts$group_col, c("group", "disease_status"))

  message("Using cell type column: ", celltype_col)
  message("Using group column: ", group_col)

  panel_key <- ifelse(grepl("kidney", opts$prefix, ignore.case = TRUE), "kidney", "cerebellum")
  plot_marker_dotplot(obj, marker_panels[[panel_key]], celltype_col, opts$prefix)
  markers <- run_markers(obj, celltype_col, opts$prefix)
  de <- run_celltype_de(obj, celltype_col, group_col, opts$ident1, opts$ident2, opts$prefix)
  write_csv(de, file.path(paths$results, paste0(opts$prefix, "_all_celltype_de.csv")))
  run_go(de, opts$prefix)
  plot_composition_from_object(obj, celltype_col, group_col, opts$prefix)

  invisible(list(markers = markers, de = de))
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
