# Annotation, marker, cell-type DE and GO enrichment.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

out_dir <- file.path(paths$results, "03_annotation_de_enrichment")
fig_dir <- file.path(paths$figures, "03_annotation_de_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

objects <- data.frame(
  dataset = c("cerebellum", "kidney"),
  rds = c(
    file.path(paths$data_processed, "seurat", "cerebellum_final.rds"),
    file.path(paths$data_processed, "seurat", "kidney_final.rds")
  ),
  celltype_col = c("celltype", "celltype"),
  group_col = c("group", "group"),
  stringsAsFactors = FALSE
)

cerebellum_marker_genes <- c(
  "ATOH1", "EOMES", "PAX6", "GABRA6", "CA8", "PCP2", "SOX2", "NES",
  "DCX", "GAD1", "OLIG1", "PDGFRA", "AQP4", "P2RY12", "CLDN5"
)

kidney_marker_genes <- c(
  "SIX2", "CITED1", "WT1", "NPHS1", "NPHS2", "CUBN", "LRP2", "SLC12A1",
  "UMOD", "SLC12A3", "AQP2", "KRT8", "COL1A1", "PECAM1", "TOP2A"
)

for (i in seq_len(nrow(objects))) {
  dataset <- objects$dataset[i]
  obj_file <- objects$rds[i]
  celltype_col <- objects$celltype_col[i]
  group_col <- objects$group_col[i]

  if (!file.exists(obj_file)) {
    message("Skip ", dataset, ": object not found at ", obj_file)
    next
  }

  obj <- readRDS(obj_file)
  DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
  obj[[celltype_col]][, 1] <- droplevels(as.factor(obj[[celltype_col]][, 1]))
  Idents(obj) <- celltype_col

  marker_genes <- if (dataset == "cerebellum") cerebellum_marker_genes else kidney_marker_genes
  marker_genes <- intersect(marker_genes, rownames(obj))

  if (length(marker_genes) > 0) {
    p <- DotPlot(obj, features = marker_genes, group.by = celltype_col) +
      RotatedAxis() +
      scale_color_gradient2(low = "#3A6EA5", mid = "grey92", high = "#B23A48") +
      theme_js_atlas(8)
    ggsave(file.path(fig_dir, paste0(dataset, "_marker_dotplot.pdf")), p, width = 8, height = 4)
  }

  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  write.csv(markers, file.path(out_dir, paste0(dataset, "_all_positive_markers.csv")), row.names = FALSE)

  top_markers <- markers %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>%
    slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
    ungroup()
  write.csv(top_markers, file.path(out_dir, paste0(dataset, "_top20_markers_by_celltype.csv")), row.names = FALSE)

  if (group_col %in% colnames(obj@meta.data)) {
    obj[[group_col]][, 1] <- factor(obj[[group_col]][, 1], levels = c("ctrl", "JS"))

    composition <- obj@meta.data %>%
      count(.data[[group_col]], .data[[celltype_col]], name = "n") %>%
      group_by(.data[[group_col]]) %>%
      mutate(percent = n / sum(n) * 100) %>%
      ungroup()
    write.csv(composition, file.path(out_dir, paste0(dataset, "_celltype_composition.csv")), row.names = FALSE)

    p <- ggplot(composition, aes(x = .data[[group_col]], y = percent, fill = .data[[celltype_col]])) +
      geom_col(width = 0.7, color = "white", linewidth = 0.1) +
      labs(x = NULL, y = "Nuclei (%)", fill = "Cell type") +
      theme_js_atlas(8)
    ggsave(file.path(fig_dir, paste0(dataset, "_celltype_composition.pdf")), p, width = 5, height = 4)

    for (ct in levels(obj[[celltype_col]][, 1])) {
      sub <- subset(obj, subset = .data[[celltype_col]] == ct)
      if (!all(c("ctrl", "JS") %in% sub[[group_col]][, 1])) next
      if (min(table(sub[[group_col]][, 1])) < 20) next

      Idents(sub) <- group_col
      de <- FindMarkers(
        sub,
        ident.1 = "JS",
        ident.2 = "ctrl",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        test.use = "wilcox"
      )
      de$gene <- rownames(de)
      ct_name <- gsub("[^A-Za-z0-9]+", "_", ct)
      write.csv(de, file.path(out_dir, paste0(dataset, "_", ct_name, "_JS_vs_ctrl_DE.csv")), row.names = FALSE)

      de_sig <- de %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
      if (nrow(de_sig) == 0) next

      for (direction in c("up_in_JS", "down_in_JS")) {
        genes <- if (direction == "up_in_JS") {
          de_sig %>% filter(avg_log2FC > 0.25) %>% pull(gene)
        } else {
          de_sig %>% filter(avg_log2FC < -0.25) %>% pull(gene)
        }
        genes <- unique(genes)
        if (length(genes) < 10) next

        entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        if (nrow(entrez) == 0) next

        ego <- enrichGO(
          gene = entrez$ENTREZID,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          readable = TRUE,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
        write.csv(as.data.frame(ego), file.path(out_dir, paste0(dataset, "_", ct_name, "_", direction, "_GO_BP.csv")), row.names = FALSE)

        if (nrow(as.data.frame(ego)) > 0) {
          p <- barplot(ego, showCategory = 12, label_format = 60) +
            ggtitle(paste(dataset, ct, direction)) +
            theme_js_atlas(8)
          ggsave(file.path(fig_dir, paste0(dataset, "_", ct_name, "_", direction, "_GO_BP.pdf")), p, width = 6, height = 4)
        }
      }
    }
  }
}
