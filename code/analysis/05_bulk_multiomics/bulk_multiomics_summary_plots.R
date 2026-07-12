# Bulk RNA/proteome/metabolome summary workflow.
#
# This script expects public-safe summary tables rather than raw human data.
# It can run PCA from a feature x sample matrix, visualize differential feature
# tables, run GO enrichment for gene/protein symbols and summarize cross-omic
# overlaps.
#
# Example:
# Rscript analysis/05_bulk_multiomics/bulk_multiomics_summary_plots.R \
#   --sample-meta metadata/bulk_sample_metadata_public_template.csv \
#   --matrix data/processed/bulk_rna_normalized_matrix.csv \
#   --de-dir data/processed/bulk_de_tables \
#   --prefix bulk_multiomics

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
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

read_feature_matrix <- function(path) {
  path <- resolve_path(path)
  if (!file.exists(path)) stop("Feature matrix not found: ", path, call. = FALSE)
  mat <- read_csv(path, show_col_types = FALSE)
  feature_col <- colnames(mat)[1]
  mat <- as.data.frame(mat)
  rownames(mat) <- mat[[feature_col]]
  mat[[feature_col]] <- NULL
  as.matrix(mat)
}

run_pca_from_matrix <- function(mat, sample_meta, prefix, modality = "bulk") {
  common <- intersect(colnames(mat), sample_meta$sample_id)
  if (length(common) < 3) stop("Need at least three overlapping samples for PCA.", call. = FALSE)

  mat <- mat[, common, drop = FALSE]
  sample_meta <- sample_meta[match(common, sample_meta$sample_id), , drop = FALSE]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$sample_id <- rownames(scores)
  scores <- left_join(scores, sample_meta, by = "sample_id")
  scores$modality <- modality
  write_csv(scores, file.path(paths$results, paste0(prefix, "_", modality, "_pca_scores.csv")))

  p <- ggplot(scores, aes(PC1, PC2, color = group, shape = tissue, label = sample_id)) +
    geom_point(size = 3) +
    theme_js_atlas() +
    labs(title = paste(modality, "PCA"))
  save_panel(p, paste0(prefix, "_", modality, "_pca.pdf"), width = 5, height = 4)
  pca
}

standardize_de_table <- function(tbl) {
  names(tbl) <- gsub("^avg_log2FC$", "log2FC", names(tbl))
  names(tbl) <- gsub("^p_val_adj$|^FDR$|^adj.P.Val$", "padj", names(tbl))
  if (!"gene" %in% names(tbl)) names(tbl)[1] <- "gene"
  if (!"log2FC" %in% names(tbl)) stop("DE table needs a log2FC or avg_log2FC column.", call. = FALSE)
  if (!"padj" %in% names(tbl)) stop("DE table needs a padj/FDR/p_val_adj column.", call. = FALSE)
  tbl
}

plot_volcano <- function(tbl, prefix, label_top = 12) {
  tbl <- standardize_de_table(tbl) %>%
    mutate(
      direction = case_when(
        padj < 0.05 & log2FC > 0.3 ~ "Up",
        padj < 0.05 & log2FC < -0.3 ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10 = -log10(pmax(padj, 1e-300))
    )

  write_csv(tbl, file.path(paths$results, paste0(prefix, "_de_standardized.csv")))
  p <- ggplot(tbl, aes(log2FC, neglog10, color = direction)) +
    geom_point(size = 0.8, alpha = 0.75) +
    scale_color_manual(values = c(Down = "#3565A8", NS = "grey82", Up = "#B5453C")) +
    theme_js_atlas() +
    labs(x = "log2 fold change", y = "-log10 adjusted P")
  save_panel(p, paste0(prefix, "_volcano.pdf"), width = 4.5, height = 4)
  tbl
}

run_go_from_de <- function(tbl, prefix) {
  tbl <- standardize_de_table(tbl)
  for (direction in c("up", "down")) {
    genes <- tbl %>%
      filter(padj < 0.05, if (direction == "up") log2FC > 0.3 else log2FC < -0.3) %>%
      pull(gene) %>%
      unique()
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
    write_csv(as.data.frame(ego), file.path(paths$results, paste0(prefix, "_", direction, "_GO_BP.csv")))
  }
}

summarize_de_directory <- function(de_dir, prefix) {
  de_dir <- resolve_path(de_dir)
  if (!dir.exists(de_dir)) {
    warning("DE directory not found; skipping: ", de_dir)
    return(invisible(NULL))
  }
  files <- list.files(de_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)
  if (length(files) == 0) return(invisible(NULL))

  all_sets <- lapply(files, function(path) {
    tbl <- if (grepl("\\.tsv$", path)) read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)
    tbl <- standardize_de_table(tbl)
    tibble(
      contrast = tools::file_path_sans_ext(basename(path)),
      up = list(unique(tbl$gene[tbl$padj < 0.05 & tbl$log2FC > 0.3])),
      down = list(unique(tbl$gene[tbl$padj < 0.05 & tbl$log2FC < -0.3]))
    )
  }) %>% bind_rows()

  overlap <- all_sets %>%
    transmute(contrast, n_up = lengths(up), n_down = lengths(down))
  write_csv(overlap, file.path(paths$results, paste0(prefix, "_de_set_sizes.csv")))
  all_sets
}

main <- function() {
  opts <- parse_args(list(
    sample_meta = "metadata/bulk_sample_metadata_public_template.csv",
    matrix = "",
    de_table = "",
    de_dir = "",
    prefix = "bulk_multiomics",
    modality = "bulk"
  ))

  if (opts$sample_meta != "" && opts$matrix != "") {
    sample_meta <- read_csv(resolve_path(opts$sample_meta), show_col_types = FALSE)
    mat <- read_feature_matrix(opts$matrix)
    run_pca_from_matrix(mat, sample_meta, opts$prefix, opts$modality)
  } else {
    message("Skipping PCA because --sample-meta or --matrix was not provided.")
  }

  if (opts$de_table != "") {
    de_path <- resolve_path(opts$de_table)
    de <- if (grepl("\\.tsv$", de_path)) read_tsv(de_path, show_col_types = FALSE) else read_csv(de_path, show_col_types = FALSE)
    de <- plot_volcano(de, opts$prefix)
    run_go_from_de(de, opts$prefix)
  }

  if (opts$de_dir != "") {
    summarize_de_directory(opts$de_dir, opts$prefix)
  }
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
