# pySCENIC downstream analysis for regulon activity, RSS and TF-target plots.
#
# Source scripts:
# - code/Rpr0/cebJS/pyscenic/ana928/scenicceb.R
# - code/Rpr0/kid0326/pyscenic/rana/pyscenic.kid.downstream.R

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggraph)
  library(tidygraph)
})

out_dir <- file.path(paths$results, "pyscenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_scenic_downstream <- function(seurat_file, loom_file, adj_tsv, prefix,
                                  subset_celltypes = NULL,
                                  rss_group_cols = c("celltype", "group", "genotype"),
                                  rss_z = c(celltype = 2, group = 1, genotype = 1)) {
  sce <- readRDS(seurat_file)
  grn <- read_tsv(adj_tsv, show_col_types = FALSE)
  loom <- open_loom(loom_file)
  on.exit(close_loom(loom), add = TRUE)

  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

  if (!is.null(subset_celltypes)) {
    sce <- subset(sce, celltype %in% subset_celltypes)
  }
  regulonAUC <- regulonAUC[, match(colnames(sce), colnames(regulonAUC))]
  regulon_names <- regulonAUC@NAMES
  sce@meta.data <- cbind(sce@meta.data, t(assay(regulonAUC[regulon_names, ])))

  for (group_col in rss_group_cols) {
    if (!group_col %in% colnames(sce@meta.data)) next
    cellinfo <- data.frame(group = sce@meta.data[[group_col]], row.names = colnames(sce))
    sub_auc <- regulonAUC[, colnames(regulonAUC) %in% rownames(cellinfo)]
    rss <- calcRSS(AUC = getAUC(sub_auc), cellAnnotation = cellinfo[colnames(sub_auc), "group"])
    cutoff <- rss_z[[group_col]]
    if (is.null(cutoff)) cutoff <- 1
    rssPlot <- plotRSS(
      rss,
      zThreshold = cutoff,
      cluster_columns = FALSE,
      order_rows = TRUE,
      thr = 0.1,
      varName = group_col,
      col.low = "#330066",
      col.mid = "#66CC66",
      col.high = "#FFCC33"
    )
    write.csv(rssPlot$df, file.path(out_dir, paste0(prefix, "_", group_col, "_RSS.csv")), row.names = FALSE)
    pdf(file.path(out_dir, paste0(prefix, "_", group_col, "_RSS.pdf")), width = 5, height = 8)
    print(rssPlot$plot)
    dev.off()

    cellsPerGroup <- split(rownames(cellinfo), cellinfo$group)
    activity_by_group <- sapply(cellsPerGroup, function(cells) rowMeans(getAUC(sub_auc)[, cells]))
    activity_scaled <- na.omit(t(scale(t(activity_by_group))))
    pdf(file.path(out_dir, paste0(prefix, "_", group_col, "_regulon_activity_heatmap.pdf")), width = 8, height = 12)
    print(draw(Heatmap(activity_scaled, name = "Regulon activity",
                       row_names_gp = grid::gpar(fontsize = 6))))
    dev.off()
  }

  # GO and network plots for the RSS-selected regulons from the last RSS table.
  rss_genes <- unique(gsub("\\(\\+\\)", "", rownames(rss)))
  pdf(file.path(out_dir, paste0(prefix, "_TF_GO_top300_targets.pdf")), width = 10, height = 5)
  for (tf in rss_genes) {
    sub_grn <- grn %>% filter(TF == tf) %>% group_by(TF) %>% top_n(300, importance) %>% ungroup()
    entrez <- bitr(sub_grn$target, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (nrow(entrez) > 5) {
      ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
      print(dotplot(ego, showCategory = 10, label_format = 60) + ggtitle(tf))
    }
  }
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_TF_target_network_top50.pdf")), width = 8, height = 7)
  for (tf in rss_genes) {
    sub_grn <- grn %>% filter(TF == tf) %>% group_by(TF) %>% top_n(50, importance) %>% ungroup()
    if (nrow(sub_grn) == 0) next
    nodes <- data.frame(node = unique(c(sub_grn$TF, sub_grn$target)))
    nodes$node.size <- ifelse(nodes$node == tf, 2, 1.5)
    edges <- data.frame(from = sub_grn$TF, to = sub_grn$target, importance = sub_grn$importance)
    graph_data <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    p <- ggraph(graph_data, layout = "stress", circular = TRUE) +
      geom_edge_arc(aes(edge_colour = importance, edge_width = importance)) +
      geom_node_point(aes(size = node.size), color = "black") +
      geom_node_label(aes(label = node), size = 3, repel = TRUE) +
      theme_void() +
      ggtitle(tf)
    print(p)
  }
  dev.off()

  sce
}

# Fill these files with local pySCENIC outputs before running.
cerebellum_seurat <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
cerebellum_loom <- file.path(paths$data, "controlled/pyscenic/cerebellum.loom")
cerebellum_adj <- file.path(paths$data, "controlled/pyscenic/cerebellum_adj.sample.tsv")
if (file.exists(cerebellum_seurat) && file.exists(cerebellum_loom) && file.exists(cerebellum_adj)) {
  run_scenic_downstream(
    cerebellum_seurat,
    cerebellum_loom,
    cerebellum_adj,
    prefix = "cerebellum",
    subset_celltypes = c("Cellcycle", "PKCs", "UBCs", "VZP", "GCs", "INs"),
    rss_group_cols = c("celltype", "group", "name"),
    rss_z = c(celltype = 2, group = 1, name = 1)
  )
}

kidney_seurat <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
kidney_loom <- file.path(paths$data, "controlled/pyscenic/kidney.loom")
kidney_adj <- file.path(paths$data, "controlled/pyscenic/kidney_adj.sample.tsv")
if (file.exists(kidney_seurat) && file.exists(kidney_loom) && file.exists(kidney_adj)) {
  run_scenic_downstream(
    kidney_seurat,
    kidney_loom,
    kidney_adj,
    prefix = "kidney",
    subset_celltypes = c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT", "UB_CD", "Stromal", "Endo"),
    rss_group_cols = c("celltype", "group", "geneotype"),
    rss_z = c(celltype = 2.5, group = 1, geneotype = 1)
  )
}
