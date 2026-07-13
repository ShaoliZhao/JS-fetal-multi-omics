# Monocle3 graph-based pseudotime analysis.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(monocle3)
  library(igraph)
  library(pheatmap)
  library(viridis)
})

out_dir <- file.path(paths$results, "06_lineage_velocity", "monocle3")
fig_dir <- file.path(paths$figures, "06_lineage_velocity", "monocle3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(paths$data_processed, "seurat", "kidney_final.rds")
lineage_name <- "kidney_nephron_monocle3"
celltype_col <- "celltype"
group_col <- "group"
lineage_celltypes <- c("NPC", "Podocyte", "PT", "LOH", "PEC", "DCT", "UB_CD", "Cycling")
root_cluster <- "5"

obj <- readRDS(object_file)
obj <- subset(obj, subset = celltype %in% lineage_celltypes)
obj[[celltype_col]][, 1] <- droplevels(obj[[celltype_col]][, 1])
obj <- JoinLayers(obj)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- counts[VariableFeatures(obj), ]

gene_annotation <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts)
)
cell_metadata <- obj@meta.data[colnames(counts), , drop = FALSE]
cell_metadata$barcode <- rownames(cell_metadata)

cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

if ("umap" %in% names(obj@reductions)) {
  reducedDims(cds)[["UMAP"]] <- obj@reductions$umap@cell.embeddings[colnames(cds), ]
  cds@clusters$UMAP$clusters <- obj$seurat_clusters[colnames(cds)]
  cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)), levels = 1)
  names(cds@clusters$UMAP$partitions) <- colnames(cds)
} else {
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
}

cds <- learn_graph(cds)
colData(cds)$cluster <- cds@clusters$UMAP$clusters

cell_ids <- which(colData(cds)$cluster == root_cluster)
if (length(cell_ids) > 0) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
} else {
  cds <- order_cells(cds)
}

pdf(file.path(fig_dir, paste0(lineage_name, "_cells_pseudotime_celltype.pdf")), width = 10, height = 5)
print(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE))
print(plot_cells(cds, color_cells_by = celltype_col, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE))
dev.off()

pt_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
pt_res <- pt_res[order(pt_res$morans_I, decreasing = TRUE), ]
write.csv(pt_res, file.path(out_dir, paste0(lineage_name, "_principal_graph_genes.csv")))

trajectory_genes <- rownames(subset(pt_res, q_value < 0.01))
if (length(trajectory_genes) > 0) {
  pdf(file.path(fig_dir, paste0(lineage_name, "_top_genes_pseudotime.pdf")), width = 8, height = 6)
  print(plot_genes_in_pseudotime(cds[head(trajectory_genes, 12), ], min_expr = 0.1, ncol = 3))
  dev.off()
}

if (length(trajectory_genes) > 20 && group_col %in% colnames(colData(cds))) {
  pt <- pseudotime(cds)
  obj$pseudotime <- pt[colnames(obj)]
  obj$ptime_bin <- floor(obj$pseudotime / 2.5)

  avg_exp <- AverageExpression(
    obj,
    assays = "RNA",
    features = head(trajectory_genes, 200),
    group.by = "ptime_bin"
  )[[1]]
  avg_exp <- avg_exp[apply(avg_exp, 1, max, na.rm = TRUE) > 0, , drop = FALSE]
  avg_exp <- t(apply(avg_exp, 1, function(x) x / max(x, na.rm = TRUE)))

  pdf(file.path(fig_dir, paste0(lineage_name, "_trajectory_gene_heatmap.pdf")), width = 5, height = 8)
  pheatmap(avg_exp, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA, color = viridis::magma(50))
  dev.off()
}

saveRDS(cds, file.path(out_dir, paste0(lineage_name, "_monocle3_cds.rds")))
