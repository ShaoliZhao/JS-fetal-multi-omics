# Bulk RNA-seq DESeq2 and GO analysis used for manuscript summary figures.
#
# Source scripts:
# - bulk0217.R: kidney bulk RNA-seq
# - 1010bulk.R: cerebellum bulk RNA-seq
# - cebPFC.R: forebrain/PFC overlap summaries
#
# This is intentionally written as a simple step-by-step script. Replace the
# featureCounts input folders below with local controlled-access files.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

out_dir <- file.path(paths$results, "bulk_rnaseq")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_featurecounts_folder <- function(folder) {
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  stopifnot(length(files) > 0)
  count_list <- lapply(files, function(f) {
    x <- read.delim(f, comment.char = "#", check.names = FALSE)
    data.frame(Geneid = x$Geneid, Length = x$Length, count = x[[7]], check.names = FALSE)
  })
  sample_ids <- tools::file_path_sans_ext(basename(files))
  exp <- count_list[[1]][, c("Geneid", "Length")]
  for (i in seq_along(count_list)) exp[[sample_ids[i]]] <- count_list[[i]]$count
  exp
}

run_deseq2 <- function(exp, sample_ids, groups, prefix, logfc_cutoff = 0.5,
                       p_col = "padj", p_cutoff = 0.05) {
  mtx <- as.matrix(exp[, sample_ids])
  rownames(mtx) <- exp$Geneid
  mtx <- mtx[rowSums(mtx) > 10, ]
  meta <- data.frame(sample = sample_ids, group = factor(groups, levels = c("ctrl", "JS")))
  rownames(meta) <- sample_ids

  dds <- DESeqDataSetFromMatrix(countData = round(mtx), colData = meta, design = ~ group)
  dds <- DESeq(dds, fitType = "mean")
  res <- as.data.frame(results(dds, contrast = c("group", "JS", "ctrl")))
  res$ENSEMBL <- rownames(res)
  write.csv(res, file.path(out_dir, paste0(prefix, "_DESeq2_all_genes.csv")), row.names = FALSE)

  pvalues <- if (p_col == "pvalue") res$pvalue else res$padj
  deg <- res %>%
    filter(!is.na(pvalues), pvalues < p_cutoff, abs(log2FoldChange) > logfc_cutoff)
  write.csv(deg, file.path(out_dir, paste0(prefix, "_DEG_", p_col, p_cutoff, "_logFC", logfc_cutoff, ".csv")), row.names = FALSE)

  gene_map <- bitr(deg$ENSEMBL, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  if (nrow(gene_map) > 5) {
    ego <- enrichGO(gene_map$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
    write.csv(as.data.frame(ego), file.path(out_dir, paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
    pdf(file.path(out_dir, paste0(prefix, "_GO_BP_barplot.pdf")), width = 9, height = 5)
    print(barplot(ego, showCategory = 10))
    dev.off()
  }

  res
}

# Cerebellum bulk RNA-seq: source 1010bulk.R
ceb_featurecounts <- file.path(paths$data, "controlled/bulk/cerebellum_featureCounts")
if (dir.exists(ceb_featurecounts)) {
  ceb_exp <- read_featurecounts_folder(ceb_featurecounts)
  write.csv(ceb_exp, file.path(out_dir, "cerebellum_featureCounts_matrix.csv"), row.names = FALSE)

  ceb_samples <- colnames(ceb_exp)[3:ncol(ceb_exp)]
  # Fill this vector to match the sample order in `ceb_samples`.
  # Original 1010bulk.R used 32 samples and compared JS vs control.
  ceb_groups <- c(
    rep("JS", 16),
    rep("ctrl", length(ceb_samples) - 16)
  )
  ceb_res <- run_deseq2(
    ceb_exp,
    sample_ids = ceb_samples,
    groups = ceb_groups,
    prefix = "cerebellum_JS_vs_ctrl",
    logfc_cutoff = 0.5,
    p_col = "padj",
    p_cutoff = 0.05
  )
}

# Kidney bulk RNA-seq: source bulk0217.R
kidney_featurecounts <- file.path(paths$data, "controlled/bulk/kidney_featureCounts")
if (dir.exists(kidney_featurecounts)) {
  kidney_exp <- read_featurecounts_folder(kidney_featurecounts)
  write.csv(kidney_exp, file.path(out_dir, "kidney_featureCounts_matrix.csv"), row.names = FALSE)

  kidney_samples <- colnames(kidney_exp)[3:ncol(kidney_exp)]
  kidney_groups <- c(rep("JS", 8), rep("ctrl", length(kidney_samples) - 8))

  kidney_all <- run_deseq2(
    kidney_exp,
    sample_ids = kidney_samples,
    groups = kidney_groups,
    prefix = "kidney_all_JS_vs_ctrl",
    logfc_cutoff = 0.5,
    p_col = "pvalue",
    p_cutoff = 0.05
  )

  # Genotype-specific contrasts from the original kidney script used
  # pvalue < 0.05 and abs(log2FC) > 0.5.
  genotype_sets <- list(
    CEP290 = list(samples = c("wxj9", "lyk10", "GSC1_4", "GSC1_5", "GSC2_1", "GSC2_2", "GSC3_1", "GSC3_2", "GSC4_1", "GSC4_2"),
                 groups = c(rep("JS", 2), rep("ctrl", 8))),
    KIAA0586 = list(samples = c("lyk10", "lyk12", "lyk9", "GSC1_4", "GSC1_5", "GSC2_1", "GSC2_2", "GSC3_1", "GSC3_2", "GSC4_1", "GSC4_2"),
                   groups = c(rep("JS", 3), rep("ctrl", 8))),
    TMEM67 = list(samples = c("yhk2_3", "lhk1", "lhk2", "GSC1_4", "GSC1_5", "GSC2_1", "GSC2_2", "GSC3_1", "GSC3_2", "GSC4_1", "GSC4_2"),
                 groups = c(rep("JS", 3), rep("ctrl", 8)))
  )

  for (nm in names(genotype_sets)) {
    use <- genotype_sets[[nm]]
    use_samples <- intersect(use$samples, colnames(kidney_exp))
    if (length(use_samples) == length(use$samples)) {
      run_deseq2(
        kidney_exp,
        sample_ids = use$samples,
        groups = use$groups,
        prefix = paste0("kidney_", nm, "_vs_ctrl"),
        logfc_cutoff = 0.5,
        p_col = "pvalue",
        p_cutoff = 0.05
      )
    }
  }
}

# PFC / forebrain overlaps: source cebPFC.R
# The original script tested several manually selected sample subsets and then
# used ggVennDiagram to compare cerebellum, kidney and cortical bulk DE genes.
# For the public repository, keep the final overlap step as table-driven code.
overlap_files <- list.files(out_dir, pattern = "_DEG_.*\\.csv$", full.names = TRUE)
if (length(overlap_files) > 1) {
  deg_sets <- lapply(overlap_files, function(f) {
    x <- read.csv(f, stringsAsFactors = FALSE)
    unique(x$ENSEMBL)
  })
  names(deg_sets) <- sub("_DEG_.*", "", basename(overlap_files))
  saveRDS(deg_sets, file.path(out_dir, "bulk_DEG_gene_sets_for_overlap.rds"))
}
