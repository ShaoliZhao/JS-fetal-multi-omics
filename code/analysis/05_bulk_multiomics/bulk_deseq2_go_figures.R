# Bulk RNA-seq, proteomics and metabolomics figure workflows.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(pheatmap)
  library(RColorBrewer)
  library(ggVennDiagram)
})

out_dir <- file.path(paths$results, "05_bulk_multiomics")
fig_dir <- file.path(paths$figures, "05_bulk_multiomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

logfc_cutoff <- 0.5
padj_cutoff <- 0.05

bulk_diff_files <- data.frame(
  tissue = c("cerebellum", "kidney", "frontal_lobe", "occipital_lobe"),
  file = c(
    file.path(paths$data_processed, "bulk", "cerebellum", "diffbox.csv"),
    file.path(paths$data_processed, "bulk", "kidney", "diffbox.csv"),
    file.path(paths$data_processed, "bulk", "PFC", "FCdiffbox.csv"),
    file.path(paths$data_processed, "bulk", "PFC", "OLdiffbox.csv")
  ),
  id_type = c("ENSEMBL", "ENSEMBL", "ENSEMBL", "ENSEMBL"),
  stringsAsFactors = FALSE
)

bulk_go_sets <- list()

for (i in seq_len(nrow(bulk_diff_files))) {
  tissue <- bulk_diff_files$tissue[i]
  file <- bulk_diff_files$file[i]
  id_type <- bulk_diff_files$id_type[i]

  if (!file.exists(file)) {
    message("Skip ", tissue, ": DEG table not found at ", file)
    next
  }

  deg <- read.csv(file, check.names = FALSE)
  gene_col <- if ("gene" %in% colnames(deg)) "gene" else colnames(deg)[1]
  deg$gene_id <- deg[[gene_col]]

  js_up <- deg %>% filter(padj < padj_cutoff & log2FoldChange > logfc_cutoff)
  ctrl_up <- deg %>% filter(padj < padj_cutoff & log2FoldChange < -logfc_cutoff)

  write.csv(js_up, file.path(out_dir, paste0(tissue, "_JS_up_logFC0.5_padj0.05.csv")), row.names = FALSE)
  write.csv(ctrl_up, file.path(out_dir, paste0(tissue, "_ctrl_up_logFC0.5_padj0.05.csv")), row.names = FALSE)

  bulk_go_sets[[paste0(tissue, "_JS_up")]] <- js_up$gene_id
  bulk_go_sets[[paste0(tissue, "_ctrl_up")]] <- ctrl_up$gene_id

  for (direction in c("JS_up", "ctrl_up")) {
    genes <- if (direction == "JS_up") js_up$gene_id else ctrl_up$gene_id
    genes <- unique(genes[!is.na(genes)])
    if (length(genes) < 10) next

    entrez <- bitr(genes, fromType = id_type, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
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
    write.csv(as.data.frame(ego), file.path(out_dir, paste0(tissue, "_", direction, "_GO_BP.csv")), row.names = FALSE)

    if (nrow(as.data.frame(ego)) > 0) {
      p <- barplot(ego, showCategory = 15, label_format = 60) +
        ggtitle(paste(tissue, direction, "GO BP")) +
        theme_js_atlas(8)
      ggsave(file.path(fig_dir, paste0(tissue, "_", direction, "_GO_BP.pdf")), p, width = 6, height = 4.5)
    }
  }
}

if (length(bulk_go_sets) >= 2) {
  pdf(file.path(fig_dir, "bulk_DEG_venn_sets.pdf"), width = 7, height = 6)
  print(ggVennDiagram(bulk_go_sets, label_alpha = 0))
  dev.off()
}

bulk_summary_xlsx <- file.path(paths$data_processed, "bulk", "bulk4organs.xlsx")
if (file.exists(bulk_summary_xlsx)) {
  summary_sheets <- c("ceb", "kid", "FC", "OL")
  go_summary <- data.frame()

  for (sheet in summary_sheets) {
    tab <- read_xlsx(bulk_summary_xlsx, sheet = sheet)
    tab$sheet <- sheet
    go_summary <- bind_rows(go_summary, tab)
  }

  if (all(c("GeneRatio", "pvalue", "Description", "type") %in% colnames(go_summary))) {
    go_summary <- go_summary %>%
      mutate(
        GeneRatio_value = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
        signed_ratio = ifelse(type == "down", -GeneRatio_value, GeneRatio_value),
        signed_logp = ifelse(type == "down", log10(pvalue), -log10(pvalue))
      )

    p <- ggplot(go_summary, aes(x = signed_ratio, y = reorder(Description, signed_ratio))) +
      geom_segment(aes(x = 0, xend = signed_ratio, yend = Description), color = "grey75", linewidth = 0.35) +
      geom_point(aes(size = Count, color = signed_logp)) +
      facet_wrap(~sheet, scales = "free_y") +
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

  deg_sheets <- c(
    kidney = "kidneydeg",
    cerebellum = "cebdeg",
    frontal_lobe = "FCdeg",
    occipital_lobe = "OLdeg"
  )
  lfc_mat <- matrix(0, nrow = length(cilia_genes), ncol = length(deg_sheets))
  rownames(lfc_mat) <- cilia_genes
  colnames(lfc_mat) <- names(deg_sheets)

  for (nm in names(deg_sheets)) {
    tab <- read_xlsx(bulk_summary_xlsx, sheet = deg_sheets[[nm]])
    gene_col <- if ("gene" %in% colnames(tab)) "gene" else colnames(tab)[1]
    lfc_col <- if ("log2FoldChange" %in% colnames(tab)) "log2FoldChange" else "avg_log2FC"
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

proteomics_files <- data.frame(
  tissue = c("kidney", "cerebellum"),
  file = c(
    file.path(paths$data_processed, "proteomics", "JSKVSctrlK_diff_annotation.xlsx"),
    file.path(paths$data_processed, "proteomics", "JSVSctrl_diff_annotation.xlsx")
  ),
  stringsAsFactors = FALSE
)

protein_sets <- list()
for (i in seq_len(nrow(proteomics_files))) {
  tissue <- proteomics_files$tissue[i]
  file <- proteomics_files$file[i]
  if (!file.exists(file)) next

  tab <- read_xlsx(file)
  gene_col <- if ("Gene" %in% colnames(tab)) "Gene" else colnames(tab)[1]
  if (!"log2FC" %in% colnames(tab) && all(c("mean_JS", "mean_ctrl") %in% colnames(tab))) {
    tab$log2FC <- log2(tab$mean_JS / tab$mean_ctrl)
  }

  protein_up <- tab %>% filter(Significant == "yes" & log2FC > logfc_cutoff)
  protein_down <- tab %>% filter(Significant == "yes" & log2FC < -logfc_cutoff)
  protein_sets[[paste0(tissue, "_protein_JS_up")]] <- protein_up[[gene_col]]
  protein_sets[[paste0(tissue, "_protein_ctrl_up")]] <- protein_down[[gene_col]]

  write.csv(protein_up, file.path(out_dir, paste0(tissue, "_protein_JS_up.csv")), row.names = FALSE)
  write.csv(protein_down, file.path(out_dir, paste0(tissue, "_protein_ctrl_up.csv")), row.names = FALSE)
}

if (length(protein_sets) >= 2) {
  pdf(file.path(fig_dir, "proteomics_venn_sets.pdf"), width = 6, height = 5)
  print(ggVennDiagram(protein_sets, label_alpha = 0))
  dev.off()
}
