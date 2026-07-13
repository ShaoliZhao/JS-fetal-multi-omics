# Monocle2 DDRTree lineage analysis.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(monocle)
  library(ggridges)
  library(reshape2)
  library(VGAM)
})

out_dir <- file.path(paths$results, "06_lineage_velocity", "monocle2")
fig_dir <- file.path(paths$figures, "06_lineage_velocity", "monocle2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(paths$data_processed, "seurat", "kidney_final.rds")
lineage_name <- "kidney_nephron"
celltype_col <- "celltype"
group_col <- "group"
sample_col <- "orig.ident"
lineage_celltypes <- c("NPC", "Podocyte", "PT", "LOH", "PEC", "DCT")
max_cells_per_sample <- 500
ordering_gene_number <- 500
root_state <- 3

obj <- readRDS(object_file)
obj <- subset(obj, subset = celltype %in% lineage_celltypes)
obj[[celltype_col]][, 1] <- droplevels(factor(obj[[celltype_col]][, 1], levels = lineage_celltypes))

set.seed(12345)
cells_by_sample <- split(colnames(obj), obj[[sample_col]][, 1])
selected_cells <- unlist(lapply(cells_by_sample, function(x) {
  if (length(x) > max_cells_per_sample) sample(x, max_cells_per_sample) else x
}))
obj <- subset(obj, cells = selected_cells)
obj <- JoinLayers(obj)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- as(as.matrix(counts), "sparseMatrix")

cell_meta <- obj@meta.data
pd <- new("AnnotatedDataFrame", data = cell_meta)

gene_meta <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts)
)
fd <- new("AnnotatedDataFrame", data = gene_meta)

cds <- newCellDataSet(
  counts,
  phenoData = pd,
  featureData = fd,
  expressionFamily = VGAM::negbinomial.size(),
  lowerDetectionLimit = 1
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))
pData(cds)[[celltype_col]] <- as.factor(pData(cds)[[celltype_col]])

diff_test_res <- differentialGeneTest(
  cds[expressed_genes, ],
  fullModelFormulaStr = paste0("~", celltype_col)
)
write.csv(diff_test_res, file.path(out_dir, paste0(lineage_name, "_ordering_gene_test.csv")))

ordering_genes <- rownames(diff_test_res)[order(diff_test_res$qval)][seq_len(min(ordering_gene_number, nrow(diff_test_res)))]
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds, max_components = 2, num_dim = 15, method = "DDRTree")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = root_state)

pdf(file.path(fig_dir, paste0(lineage_name, "_trajectory_celltype_state_pseudotime.pdf")), width = 8, height = 10)
print(plot_cell_trajectory(cds, color_by = celltype_col))
print(plot_cell_trajectory(cds, color_by = "State"))
print(plot_cell_trajectory(cds, color_by = "Pseudotime"))
if (group_col %in% colnames(pData(cds))) {
  print(plot_cell_trajectory(cds, color_by = "Pseudotime") + facet_wrap(as.formula(paste("~", group_col))))
}
dev.off()

if (group_col %in% colnames(pData(cds))) {
  pdata <- pData(cds)

  p <- ggplot(pdata, aes(x = Pseudotime, fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    scale_fill_manual(values = c("ctrl" = "#4DBBD5FF", "JS" = "#E64B35FF")) +
    scale_color_manual(values = c("ctrl" = "#4DBBD5FF", "JS" = "#E64B35FF")) +
    labs(x = "Pseudotime", y = "Cell density") +
    theme_js_atlas(9)
  ggsave(file.path(fig_dir, paste0(lineage_name, "_pseudotime_density_JS_ctrl.pdf")), p, width = 4.5, height = 3.5)

  state_prop <- pdata %>%
    count(.data[[group_col]], State, name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()
  write.csv(state_prop, file.path(out_dir, paste0(lineage_name, "_state_proportion.csv")), row.names = FALSE)
}

target_genes <- c(
  "WT1", "NPHS2", "SLC12A1", "DCDC2", "YAP1", "GLS", "EPB41L5",
  "MAGI2", "PKP4", "AFDN", "ROBO2", "PBX1", "TCF7L2", "MAML2"
)
valid_genes <- rownames(fData(cds))[fData(cds)$gene_short_name %in% target_genes]
if (length(valid_genes) > 0) {
  pdf(file.path(fig_dir, paste0(lineage_name, "_selected_genes_pseudotime.pdf")), width = 8, height = 7)
  print(plot_genes_in_pseudotime(cds[valid_genes, ], color_by = group_col, ncol = 3, cell_size = 0))
  dev.off()
}

saveRDS(cds, file.path(out_dir, paste0(lineage_name, "_monocle2_cds.rds")))
