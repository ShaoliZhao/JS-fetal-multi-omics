# pySCENIC / SCENIC downstream analysis.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggraph)
  library(tidygraph)
})

out_dir <- file.path(paths$results, "pyscenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

add_regulon_auc <- function(sce, loom_file) {
  loom <- open_loom(loom_file)
  on.exit(close_loom(loom), add = TRUE)
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  regulonAUC <- regulonAUC[, match(colnames(sce), colnames(regulonAUC))]
  sce@meta.data <- cbind(sce@meta.data, t(assay(regulonAUC[regulonAUC@NAMES, ])))
  list(object = sce, regulons = regulons, regulonAUC = regulonAUC)
}

run_rss_block <- function(regulonAUC, metadata, group_col, prefix,
                          z_threshold = 1, width = 5, height = 8) {
  cellTypes <- data.frame(group = metadata[[group_col]], row.names = rownames(metadata))
  sub_regulonAUC <- regulonAUC[, colnames(regulonAUC) %in% rownames(cellTypes)]
  rss <- calcRSS(AUC = getAUC(sub_regulonAUC),
                 cellAnnotation = cellTypes[colnames(sub_regulonAUC), "group"])
  rssPlot <- plotRSS(
    rss,
    zThreshold = z_threshold,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr = 0.1,
    varName = group_col,
    col.low = "#330066",
    col.mid = "#66CC66",
    col.high = "#FFCC33"
  )
  write.csv(rssPlot$df, file.path(out_dir, paste0(prefix, "_", group_col, "_RSS.csv")),
            row.names = FALSE)
  pdf(file.path(out_dir, paste0(prefix, "_", group_col, "_RSS.pdf")), width = width, height = height)
  print(rssPlot$plot)
  dev.off()

  cellsPerGroup <- split(rownames(cellTypes), cellTypes$group)
  activity_by_group <- sapply(cellsPerGroup, function(cells) rowMeans(getAUC(sub_regulonAUC)[, cells]))
  write.csv(activity_by_group,
            file.path(out_dir, paste0(prefix, "_", group_col, "_regulon_activity_mean.csv")))
  pdf(file.path(out_dir, paste0(prefix, "_", group_col, "_activity_heatmap.pdf")),
      width = width, height = height)
  pheatmap(activity_by_group,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("white", "red"))(50),
           main = paste(prefix, group_col, "regulon activity"),
           fontsize = 10,
           scale = "none")
  dev.off()
  rss
}

plot_tf_go_and_network <- function(grn, regulons, rss, prefix, top_target = 300) {
  tfs <- gsub("\\(\\+\\)", "", rownames(rss))
  sub_grn_all <- grn %>%
    filter(TF %in% tfs) %>%
    group_by(TF) %>%
    slice_max(order_by = importance, n = top_target, with_ties = FALSE) %>%
    ungroup()
  write.csv(sub_grn_all, file.path(out_dir, paste0(prefix, "_RSS_TF_top", top_target, "_targets.csv")),
            row.names = FALSE)

  pdf(file.path(out_dir, paste0(prefix, "_TF_GO_top", top_target, "_targets.pdf")),
      width = 10, height = 5)
  for (tf in unique(sub_grn_all$TF)) {
    entrez <- suppressMessages(bitr(
      sub_grn_all$target[sub_grn_all$TF == tf],
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ))
    if (nrow(entrez) > 5) {
      ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
      print(dotplot(ego, showCategory = 10, label_format = 60) + ggtitle(tf))
    }
  }
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_TF_network_top50.pdf")), width = 8, height = 7)
  for (tf in unique(sub_grn_all$TF)) {
    sub_grn <- sub_grn_all %>% filter(TF == tf) %>% slice_max(importance, n = 50)
    if (nrow(sub_grn) == 0) next
    nodes <- data.frame(node = unique(c(sub_grn$TF, sub_grn$target)))
    edges <- data.frame(from = sub_grn$TF, to = sub_grn$target, importance = sub_grn$importance)
    graph_data <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    p <- ggraph(graph_data, layout = "stress", circular = TRUE) +
      geom_edge_arc(aes(edge_colour = importance, edge_width = importance)) +
      geom_node_point(size = 1.5, color = "black") +
      geom_node_label(aes(label = node), size = 3, repel = TRUE) +
      theme_void() +
      ggtitle(tf)
    print(p)
  }
  dev.off()
}

run_scenic_downstream <- function(seurat_file, loom_file, adj_tsv, prefix,
                                  subset_celltypes,
                                  rss_settings,
                                  name_recode = NULL) {
  sce <- readRDS(seurat_file)
  grn <- read_tsv(adj_tsv, show_col_types = FALSE)

  if (!is.null(name_recode) && "name" %in% colnames(sce@meta.data)) {
    for (nm in names(name_recode)) sce$name[sce$name %in% name_recode[[nm]]] <- nm
  }
  if (!is.null(subset_celltypes)) sce <- subset(sce, celltype %in% subset_celltypes)
  sce$grct <- paste(sce$celltype, sce$group, sep = "_")

  scenic <- add_regulon_auc(sce, loom_file)
  sce <- scenic$object
  regulonAUC <- scenic$regulonAUC
  regulons <- scenic$regulons

  last_rss <- NULL
  for (setting in rss_settings) {
    if (!setting$group_col %in% colnames(sce@meta.data)) next
    last_rss <- run_rss_block(
      regulonAUC,
      sce@meta.data,
      group_col = setting$group_col,
      prefix = prefix,
      z_threshold = setting$z,
      width = setting$width,
      height = setting$height
    )
  }
  if (!is.null(last_rss)) plot_tf_go_and_network(grn, regulons, last_rss, prefix)
  saveRDS(sce, file.path(out_dir, paste0(prefix, "_Seurat_with_regulonAUC.rds")))
  sce
}

cerebellum_seurat <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
cerebellum_loom <- file.path(paths$data, "controlled/pyscenic/cerebellum_cc14all928.loom")
cerebellum_adj <- file.path(paths$data, "controlled/pyscenic/cerebellum_adj.sample.tsv")
if (file.exists(cerebellum_seurat) && file.exists(cerebellum_loom) && file.exists(cerebellum_adj)) {
  run_scenic_downstream(
    cerebellum_seurat,
    cerebellum_loom,
    cerebellum_adj,
    prefix = "cerebellum",
    subset_celltypes = c("Cellcycle", "PKCs", "UBCs", "VZP", "GCs", "INs"),
    name_recode = list(OFD1 = c("KIAA0586")),
    rss_settings = list(
      list(group_col = "celltype", z = 2, width = 5, height = 9),
      list(group_col = "group", z = 1, width = 4, height = 6),
      list(group_col = "name", z = 1.2, width = 4, height = 6)
    )
  )
}

kidney_seurat <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
kidney_loom <- file.path(paths$data, "controlled/pyscenic/kidney330_or_kd918.loom")
kidney_adj <- file.path(paths$data, "controlled/pyscenic/kidney_adj.sample.tsv")
if (file.exists(kidney_seurat) && file.exists(kidney_loom) && file.exists(kidney_adj)) {
  run_scenic_downstream(
    kidney_seurat,
    kidney_loom,
    kidney_adj,
    prefix = "kidney",
    subset_celltypes = c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT",
                         "UB_CD", "Stromal", "Endo", "CD", "Prolif", "TAL",
                         "PODO", "IC", "M"),
    rss_settings = list(
      list(group_col = "celltype", z = 2.5, width = 5, height = 9),
      list(group_col = "group", z = 1, width = 3, height = 5),
      list(group_col = "geneotype", z = 1.2, width = 4, height = 8)
    )
  )
}
