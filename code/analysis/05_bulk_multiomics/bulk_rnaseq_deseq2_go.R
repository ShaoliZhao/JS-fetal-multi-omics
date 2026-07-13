# Bulk RNA-seq DESeq2 and GO workflow.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggVennDiagram)
  library(pheatmap)
  library(RColorBrewer)
  library(readxl)
})

out_dir <- file.path(paths$results, "05_bulk_multiomics", "bulk_rnaseq")
fig_dir <- file.path(paths$figures, "05_bulk_multiomics", "bulk_rnaseq")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sample_sheet_file <- file.path(paths$metadata, "bulk_rnaseq_sample_sheet_template.csv")
bulk_samples <- read.csv(sample_sheet_file, stringsAsFactors = FALSE)
bulk_samples <- bulk_samples[bulk_samples$include == TRUE, ]

logfc_cutoff <- 0.5
padj_cutoff <- 0.05

analyses <- unique(bulk_samples$analysis)
bulk_deg_sets <- list()

for (analysis in analyses) {
  sample_info <- bulk_samples[bulk_samples$analysis == analysis, ]
  sample_info <- sample_info[match(sample_info$sample_id, sample_info$sample_id), ]
  sample_info$count_path <- file.path(project_root, sample_info$count_file)

  if (!all(file.exists(sample_info$count_path))) {
    missing_files <- sample_info$count_path[!file.exists(sample_info$count_path)]
    message("Skip DESeq2 rerun for ", analysis, ": missing featureCounts files. First missing file: ", missing_files[[1]])
    next
  }

  count_list <- list()
  gene_length <- NULL
  gene_id <- NULL

  for (i in seq_len(nrow(sample_info))) {
    tab <- read.table(sample_info$count_path[i], header = TRUE, check.names = FALSE)
    count_col <- if ("count_column" %in% colnames(sample_info) && !is.na(sample_info$count_column[i])) {
      as.integer(sample_info$count_column[i])
    } else {
      7
    }

    if (is.null(gene_id)) {
      gene_id <- tab$Geneid
      gene_length <- tab$Length
    }

    count_list[[sample_info$sample_id[i]]] <- tab[, count_col]
  }

  count_matrix <- do.call(cbind, count_list)
  rownames(count_matrix) <- gene_id
  count_matrix <- round(as.matrix(count_matrix))
  count_matrix <- count_matrix[rowSums(count_matrix) > 10, , drop = FALSE]

  col_data <- data.frame(
    row.names = sample_info$sample_id,
    group = factor(sample_info$group, levels = c("ctrl", "JS"))
  )

  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~group
  )
  dds <- DESeq(dds, fitType = "mean")
  res <- results(dds, contrast = c("group", "JS", "ctrl"))
  res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res$ENSEMBL <- rownames(res)
  res$length <- gene_length[match(res$ENSEMBL, gene_id)]
  res <- cbind(res, as.data.frame(counts(dds, normalized = TRUE, replaced = FALSE)))

  nametable <- bitr(
    res$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )
  diffbox <- merge(res, nametable, by = "ENSEMBL", all.x = TRUE)
  diffbox <- diffbox[order(diffbox$padj), ]
  write.csv(diffbox, file.path(out_dir, paste0(analysis, "_diffbox.csv")), row.names = FALSE)

  normalized_counts <- counts(dds, normalized = TRUE)
  write.csv(normalized_counts, file.path(out_dir, paste0(analysis, "_normalized_counts.csv")))

  vsd <- vst(dds, blind = FALSE)
  pdf(file.path(fig_dir, paste0(analysis, "_PCA_and_sample_distance.pdf")), width = 8, height = 4)
  print(plotPCA(vsd, intgroup = "group") + ggtitle(analysis))
  sample_dists <- dist(t(assay(vsd)))
  pheatmap(as.matrix(sample_dists), clustering_distance_rows = sample_dists,
           clustering_distance_cols = sample_dists, main = paste(analysis, "sample distance"))
  dev.off()

  js_up <- diffbox %>% filter(padj < padj_cutoff & log2FoldChange > logfc_cutoff)
  ctrl_up <- diffbox %>% filter(padj < padj_cutoff & log2FoldChange < -logfc_cutoff)
  write.csv(js_up, file.path(out_dir, paste0(analysis, "_JS_up_logFC0.5_padj0.05.csv")), row.names = FALSE)
  write.csv(ctrl_up, file.path(out_dir, paste0(analysis, "_ctrl_up_logFC0.5_padj0.05.csv")), row.names = FALSE)

  bulk_deg_sets[[paste0(analysis, "_JS_up")]] <- js_up$ENSEMBL
  bulk_deg_sets[[paste0(analysis, "_ctrl_up")]] <- ctrl_up$ENSEMBL

  for (direction in c("JS_up", "ctrl_up")) {
    genes <- if (direction == "JS_up") js_up$ENSEMBL else ctrl_up$ENSEMBL
    genes <- unique(genes[!is.na(genes)])
    if (length(genes) < 10) next

    entrez <- bitr(genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
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
    write.csv(as.data.frame(ego), file.path(out_dir, paste0(analysis, "_", direction, "_GO_BP.csv")), row.names = FALSE)

    if (nrow(as.data.frame(ego)) > 0) {
      p <- barplot(ego, showCategory = 15, label_format = 60) +
        labs(title = paste(analysis, direction)) +
        theme_js_atlas(8)
      ggsave(file.path(fig_dir, paste0(analysis, "_", direction, "_GO_BP.pdf")), p, width = 6, height = 4.5)
    }
  }
}

if (length(bulk_deg_sets) >= 2) {
  pdf(file.path(fig_dir, "bulk_DEG_venn_sets.pdf"), width = 7, height = 6)
  print(ggVennDiagram(bulk_deg_sets, label_alpha = 0))
  dev.off()
}

bulk_supp_xlsx <- file.path(paths$data_processed, "supplementary_tables", "Supplementary_Table_2_bulk_transcriptomics_and_celltype_composition.xlsx")
if (file.exists(bulk_supp_xlsx)) {
  supp_deg_sheets <- c(
    kidney = "Bulk_kidneydeg",
    cerebellum = "Bulk_cebdeg",
    frontal_lobe = "Bulk_FCdeg",
    occipital_lobe = "Bulk_OLdeg"
  )

  supp_deg_sets <- list()
  for (nm in names(supp_deg_sheets)) {
    tab <- read_xlsx(bulk_supp_xlsx, sheet = supp_deg_sheets[[nm]])
    js_up <- tab %>% filter(padj < padj_cutoff & log2FoldChange > logfc_cutoff)
    ctrl_up <- tab %>% filter(padj < padj_cutoff & log2FoldChange < -logfc_cutoff)
    supp_deg_sets[[paste0(nm, "_JS_up")]] <- js_up$ENSEMBL
    supp_deg_sets[[paste0(nm, "_ctrl_up")]] <- ctrl_up$ENSEMBL
  }

  if (length(supp_deg_sets) >= 2) {
    pdf(file.path(fig_dir, "bulk_DEG_venn_sets_from_supplementary_table.pdf"), width = 7, height = 6)
    print(ggVennDiagram(supp_deg_sets, label_alpha = 0))
    dev.off()
  }

  go_sheets <- data.frame(
    tissue = c("cerebellum", "cerebellum", "frontal_lobe", "frontal_lobe", "kidney", "kidney", "occipital_lobe", "occipital_lobe"),
    direction = c("ctrl", "JS", "ctrl", "JS", "ctrl", "JS", "ctrl", "JS"),
    sheet = c("ceb.ctrl.0.5", "ceb.JS.0.5", "FC.ctrl.0.5", "FC.JS.0.5", "kd.ctrl.0.5", "kd.JS.0.5", "OL.ctrl.0.5", "OL.JS.0.5"),
    stringsAsFactors = FALSE
  )
  go_summary <- data.frame()

  for (i in seq_len(nrow(go_sheets))) {
    tab <- read_xlsx(bulk_supp_xlsx, sheet = go_sheets$sheet[i])
    tab$tissue <- go_sheets$tissue[i]
    tab$direction <- go_sheets$direction[i]
    go_summary <- bind_rows(go_summary, tab)
  }

  if (all(c("GeneRatio", "pvalue", "Description", "direction") %in% colnames(go_summary))) {
    go_summary <- go_summary %>%
      mutate(
        GeneRatio_value = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
        signed_ratio = ifelse(direction == "ctrl", -GeneRatio_value, GeneRatio_value),
        signed_logp = ifelse(direction == "ctrl", log10(pvalue), -log10(pvalue))
      )

    p <- ggplot(go_summary, aes(x = signed_ratio, y = reorder(Description, signed_ratio))) +
      geom_segment(aes(x = 0, xend = signed_ratio, yend = Description), color = "grey75", linewidth = 0.35) +
      geom_point(aes(size = Count, color = signed_logp)) +
      facet_wrap(~tissue, scales = "free_y") +
      scale_color_gradient2(low = "#663E92", mid = "grey90", high = "#F6F39B", midpoint = 0) +
      labs(x = "Gene ratio (JS up positive, ctrl up negative)", y = NULL, color = "Signed -log10(P)") +
      theme_js_atlas(7)
    ggsave(file.path(fig_dir, "bulk_four_region_directional_GO_lollipop.pdf"), p, width = 9, height = 7)
  }

  cilia_gene_sets <- list(
    Axoneme = c("DNAH5", "DNAH9", "DNAI1", "DNAI2", "CCDC39", "CCDC40", "RSPH1", "RSPH4A", "RSPH9"),
    IFT = c("IFT20", "IFT27", "IFT43", "IFT46", "IFT52", "IFT57", "IFT74", "IFT80", "IFT81", "IFT88", "IFT122", "IFT140", "IFT172"),
    Transition_zone = c("CEP290", "OFD1", "TMEM67", "RPGRIP1L", "MKS1", "NPHP1", "NPHP4", "TCTN1", "TCTN2", "TCTN3"),
    Basal_body = c("CEP120", "CEP135", "CEP164", "CEP192", "CEP250", "CETN2", "CETN3", "ODF2", "PCM1", "PLK4")
  )

  cilia_genes <- unique(unlist(cilia_gene_sets))
  cilia_group <- stack(cilia_gene_sets)
  colnames(cilia_group) <- c("gene", "module")
  deg_sheets <- supp_deg_sheets

  lfc_mat <- matrix(0, nrow = length(cilia_genes), ncol = length(deg_sheets))
  rownames(lfc_mat) <- cilia_genes
  colnames(lfc_mat) <- names(deg_sheets)

  for (nm in names(deg_sheets)) {
    tab <- read_xlsx(bulk_supp_xlsx, sheet = deg_sheets[[nm]])
    gene_col <- if ("SYMBOL" %in% colnames(tab)) "SYMBOL" else if ("gene" %in% colnames(tab)) "gene" else colnames(tab)[1]
    lfc_col <- "log2FoldChange"
    tab <- tab %>% select(gene = all_of(gene_col), lfc = all_of(lfc_col))
    hit <- intersect(cilia_genes, tab$gene)
    lfc_mat[hit, nm] <- tab$lfc[match(hit, tab$gene)]
  }

  row_anno <- data.frame(Module = cilia_group$module[match(rownames(lfc_mat), cilia_group$gene)])
  rownames(row_anno) <- rownames(lfc_mat)

  pdf(file.path(fig_dir, "bulk_cilia_gene_log2FC_heatmap.pdf"), width = 5.5, height = 10)
  pheatmap(
    lfc_mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_row = row_anno,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
    breaks = seq(-2, 2, length.out = 102),
    border_color = NA
  )
  dev.off()
}
