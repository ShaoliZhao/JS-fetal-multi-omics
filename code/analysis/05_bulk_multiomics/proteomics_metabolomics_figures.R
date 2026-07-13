# Proteomics, metabolomics and RNA-protein overlap figures.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(readxl)
  library(ggVennDiagram)
  library(pheatmap)
  library(RColorBrewer)
})

out_dir <- file.path(paths$results, "05_bulk_multiomics", "proteomics_metabolomics")
fig_dir <- file.path(paths$figures, "05_bulk_multiomics", "proteomics_metabolomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

logfc_cutoff <- 0.5
supp_dir <- file.path(paths$data_processed, "supplementary_tables")
protein_xlsx <- file.path(supp_dir, "Supplementary_Table_3_proteomics.xlsx")
metabolite_xlsx <- file.path(supp_dir, "Supplementary_Table_4_metabolomics.xlsx")
bulk_xlsx <- file.path(supp_dir, "Supplementary_Table_2_bulk_transcriptomics_and_celltype_composition.xlsx")

protein_diff_sheets <- data.frame(
  tissue = c("cerebellum", "kidney"),
  sheet = c("ceb JS_ctrl_diff_annotation", "kidney JS_ctrl_diff_anno"),
  go_up_sheet = c("cerebellum GO upregulated", "kidney GO upregulated"),
  go_down_sheet = c("cerebellum GO downregulated", "kidney GO downregulated"),
  stringsAsFactors = FALSE
)

protein_sets <- list()

if (file.exists(protein_xlsx)) {
  for (i in seq_len(nrow(protein_diff_sheets))) {
    tissue <- protein_diff_sheets$tissue[i]
    tab <- read_xlsx(protein_xlsx, sheet = protein_diff_sheets$sheet[i])
    gene_col <- "gene_name"
    if (!"log2FC" %in% colnames(tab)) {
      tab$log2FC <- log2(as.numeric(tab$FC))
    }

    protein_up <- tab %>% filter(Significant == "yes" & log2FC > logfc_cutoff)
    protein_down <- tab %>% filter(Significant == "yes" & log2FC < -logfc_cutoff)
    protein_sets[[paste0(tissue, "_protein_JS_up")]] <- protein_up[[gene_col]]
    protein_sets[[paste0(tissue, "_protein_ctrl_up")]] <- protein_down[[gene_col]]

    write.csv(protein_up, file.path(out_dir, paste0(tissue, "_protein_JS_up.csv")), row.names = FALSE)
    write.csv(protein_down, file.path(out_dir, paste0(tissue, "_protein_ctrl_up.csv")), row.names = FALSE)

    go_up <- read_xlsx(protein_xlsx, sheet = protein_diff_sheets$go_up_sheet[i])
    go_down <- read_xlsx(protein_xlsx, sheet = protein_diff_sheets$go_down_sheet[i])
    write.csv(go_up, file.path(out_dir, paste0(tissue, "_protein_JS_up_GO_BP.csv")), row.names = FALSE)
    write.csv(go_down, file.path(out_dir, paste0(tissue, "_protein_ctrl_up_GO_BP.csv")), row.names = FALSE)
  }

  if (length(protein_sets) >= 2) {
    pdf(file.path(fig_dir, "proteomics_venn_sets.pdf"), width = 6, height = 5)
    print(ggVennDiagram(protein_sets, label_alpha = 0))
    dev.off()
  }
}

if (file.exists(bulk_xlsx) && length(protein_sets) > 0) {
  rna_diff_sheets <- data.frame(
    tissue = c("kidney", "cerebellum"),
    sheet = c("Bulk_kidneydeg", "Bulk_cebdeg"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(rna_diff_sheets))) {
    tissue <- rna_diff_sheets$tissue[i]
    rna <- read_xlsx(bulk_xlsx, sheet = rna_diff_sheets$sheet[i])
    rna_up <- rna %>% filter(padj < 0.05 & log2FoldChange > 0.5) %>% pull(ENSEMBL)
    rna_down <- rna %>% filter(padj < 0.05 & log2FoldChange < -0.5) %>% pull(ENSEMBL)

    protein_up <- protein_sets[[paste0(tissue, "_protein_JS_up")]]
    protein_down <- protein_sets[[paste0(tissue, "_protein_ctrl_up")]]
    if (is.null(protein_up) || is.null(protein_down)) next

    protein_up_ensembl <- bitr(protein_up, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)$ENSEMBL
    protein_down_ensembl <- bitr(protein_down, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)$ENSEMBL

    gene_sets <- list(
      RNA_JS_up = rna_up,
      RNA_ctrl_up = rna_down,
      protein_JS_up = protein_up_ensembl,
      protein_ctrl_up = protein_down_ensembl
    )

    pdf(file.path(fig_dir, paste0(tissue, "_RNA_protein_overlap_venn.pdf")), width = 6, height = 5)
    print(ggVennDiagram(gene_sets, label_alpha = 0))
    dev.off()
  }
}

if (file.exists(metabolite_xlsx)) {
  metabolite_sheets <- data.frame(
    tissue = c("cerebellum", "kidney"),
    sheet = c("cerebellum JS_ctrl", "kidney JS_ctrl"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(metabolite_sheets))) {
    tissue <- metabolite_sheets$tissue[i]
    tab <- read_xlsx(metabolite_xlsx, sheet = metabolite_sheets$sheet[i])
    if (!"log2FC" %in% colnames(tab)) {
      tab$log2FC <- log2(as.numeric(tab$FC))
    }

    sig <- tab %>% filter(Significant == "yes")
    write.csv(sig, file.path(out_dir, paste0(tissue, "_significant_metabolites.csv")), row.names = FALSE)

    sample_cols <- grep("^(JS|ctrl)", colnames(tab), value = TRUE)
    if (length(sample_cols) == 0) next

    top_metabolites <- tab %>%
      arrange(q_value) %>%
      slice_head(n = min(30, n())) %>%
      as.data.frame()
    rownames(top_metabolites) <- make.unique(top_metabolites$Metabolites)
    metabo_mat <- as.matrix(top_metabolites[, sample_cols])
    storage.mode(metabo_mat) <- "numeric"

    pdf(file.path(fig_dir, paste0(tissue, "_metabolomics_top30_heatmap.pdf")), width = 6, height = 8)
    pheatmap(
      metabo_mat,
      scale = "row",
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
      border_color = NA
    )
    dev.off()
  }
}
