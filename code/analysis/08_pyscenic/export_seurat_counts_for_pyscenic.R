# Export Seurat RNA counts for pySCENIC.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

export_counts <- function(seurat_file, prefix) {
  obj <- readRDS(seurat_file)
  DefaultAssay(obj) <- "RNA"
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  out_dir <- file.path(paths$results, "pyscenic", paste0(prefix, "_matrix"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeMM(counts, file.path(out_dir, "matrix.mtx"))
  write.table(rownames(counts), file.path(out_dir, "genes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(counts), file.path(out_dir, "cells.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

cerebellum_file <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
if (file.exists(cerebellum_file)) export_counts(cerebellum_file, "cerebellum")

kidney_file <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
if (file.exists(kidney_file)) export_counts(kidney_file, "kidney")
