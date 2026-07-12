# Export the spatial cell-bin object to pySCENIC-compatible Matrix Market files.
#
# Source provenance: code/ST_code_backup/Rpro/ana/cellbin/0127/pyscenic/changest.R.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

brain_file <- file.path(paths$results, "spatial_cellbin/cellbin_brain_merge_annotated.rds")
out_dir <- file.path(paths$results, "pyscenic/spatial_cellbin_matrix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

brain_merge <- readRDS(brain_file)
DefaultAssay(brain_merge) <- "Spatial"
counts <- GetAssayData(brain_merge, assay = "Spatial", slot = "counts")

writeMM(counts, file.path(out_dir, "matrix.mtx"))
write.table(rownames(counts), file.path(out_dir, "genes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(counts), file.path(out_dir, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

brain_merge$ctgroup <- paste0(brain_merge$celltype, "_", brain_merge$orig.ident)
set.seed(12345)
target_total <- 30000
ct_count <- as.data.frame(table(brain_merge$ctgroup))
colnames(ct_count) <- c("ctgroup", "n_total")
ct_count$weight <- sqrt(ct_count$n_total)
ct_count$n_sample <- round(target_total * ct_count$weight / sum(ct_count$weight))
ct_count$n_sample <- pmin(ct_count$n_sample, ct_count$n_total)

sampled_cells <- unlist(lapply(seq_len(nrow(ct_count)), function(i) {
  cells <- colnames(brain_merge)[brain_merge$ctgroup == ct_count$ctgroup[i]]
  sample(cells, ct_count$n_sample[i])
}))
brain_scenic_sub <- subset(brain_merge, cells = sampled_cells)
saveRDS(brain_scenic_sub, file.path(out_dir, "spatial_cellbin_balanced_for_pyscenic.rds"))
