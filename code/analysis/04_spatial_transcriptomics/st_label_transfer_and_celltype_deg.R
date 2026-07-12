# Label transfer and cell-type-specific spatial DEG for cell-bin Stereo-seq.
#
# Source provenance: code/ST_code_backup/Rpro/ana/integrate/stintegrate.R.
# The public version starts from the annotated cell-bin object produced by
# cellbin_build_and_annotation.R and the final cerebellar snRNA-seq Seurat object.

source("analysis/00_setup/project_config.R")
source("analysis/04_spatial_transcriptomics/STimport_public.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

spatial_obj_file <- file.path(paths$results, "spatial_cellbin/cellbin_brain_merge_annotated.rds")
reference_file <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
out_dir <- file.path(paths$results, "spatial_label_transfer_deg")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

brain_merge <- readRDS(spatial_obj_file)
allen_reference <- readRDS(reference_file)

allen_reference$mystype <- as.character(allen_reference$seurat_clusters)
allen_reference$mystype[allen_reference$mystype %in% c("9", "12")] <- "UBCs"
allen_reference$mystype[allen_reference$mystype %in% c("0")] <- "pGC"
allen_reference$mystype[allen_reference$mystype %in% c("3")] <- "immatureGC"
allen_reference$mystype[allen_reference$mystype %in% c("13", "8")] <- "GCP"
allen_reference$mystype[allen_reference$mystype %in% c("6")] <- "aGC"
allen_reference$mystype[allen_reference$mystype %in% c("17")] <- "DN"
allen_reference$mystype[allen_reference$mystype %in% c("5")] <- "PKC"
allen_reference$mystype[allen_reference$mystype %in% c("4", "10", "11")] <- "AST"
allen_reference$mystype[allen_reference$mystype %in% c("1")] <- "BG"
allen_reference$mystype[allen_reference$mystype %in% c("2")] <- "VZP"
allen_reference$mystype[allen_reference$mystype %in% c("7")] <- "IN1"
allen_reference$mystype[allen_reference$mystype %in% c("16")] <- "IN2"
allen_reference$mystype[allen_reference$mystype %in% c("15")] <- "IN3"
allen_reference$mystype[allen_reference$mystype %in% c("20")] <- "Endo"
allen_reference$mystype[allen_reference$mystype %in% c("14")] <- "cellcycle"
allen_reference$mystype[allen_reference$mystype %in% c("19")] <- "MG"
allen_reference$mystype[allen_reference$mystype %in% c("21")] <- "ODC"
allen_reference$mystype[allen_reference$mystype %in% c("18")] <- "OPC"

DefaultAssay(allen_reference) <- "RNA"
allen_reference <- SCTransform(allen_reference, verbose = FALSE)
allen_reference <- RunPCA(allen_reference, verbose = FALSE)
allen_reference <- RunUMAP(allen_reference, dims = 1:30, verbose = FALSE)

transfer_one_sample <- function(query, reference, sample_name) {
  DefaultAssay(query) <- "Spatial"
  query <- SCTransform(query, assay = "Spatial", verbose = FALSE)
  query <- RunPCA(query, verbose = FALSE)
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "SCT",
    dims = 1:30
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$mystype,
    prediction.assay = TRUE,
    weight.reduction = query[["pca"]],
    dims = 1:30
  )
  query[["predictions"]] <- predictions
  pred_mat <- as.matrix(query[["predictions"]]@data)
  query$predicted_class <- rownames(pred_mat)[max.col(t(pred_mat), ties.method = "first")]
  query$orig.ident <- sample_name
  query
}

ofd1 <- subset(brain_merge, orig.ident == "OFD1")
ctrl <- subset(brain_merge, orig.ident == "ctrl")
ofd1 <- transfer_one_sample(ofd1, allen_reference, "OFD1")
ctrl <- transfer_one_sample(ctrl, allen_reference, "ctrl")
brain_transfer <- merge(ctrl, ofd1)

pdf(file.path(out_dir, "label_transfer_predicted_class.pdf"), width = 10, height = 5)
print(smultiDimplot(brain_transfer, color = my25cols, group = "predicted_class", ncol = 2, pt.size = 0.1))
dev.off()
saveRDS(brain_transfer, file.path(out_dir, "cellbin_label_transfer.rds"))

DefaultAssay(brain_transfer) <- "Spatial"
brain_transfer[["RNA"]] <- CreateAssayObject(counts = GetAssayData(brain_transfer, assay = "Spatial", slot = "counts"))
DefaultAssay(brain_transfer) <- "RNA"
brain_transfer <- NormalizeData(brain_transfer)
brain_transfer <- FindVariableFeatures(brain_transfer)

Idents(brain_transfer) <- "predicted_class"
celltypes <- sort(unique(brain_transfer$predicted_class))

for (ct in celltypes) {
  sub_obj <- subset(brain_transfer, predicted_class == ct)
  if (length(unique(sub_obj$orig.ident)) < 2) next

  deg <- FindMarkers(
    sub_obj,
    ident.1 = "OFD1",
    ident.2 = "ctrl",
    group.by = "orig.ident",
    min.pct = 0,
    min.diff.pct = 0,
    logfc.threshold = 0.5
  )
  deg$gene <- rownames(deg)
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", ct)
  write.csv(deg, file.path(out_dir, paste0(safe_ct, "_OFD1_vs_ctrl_DEG_logFC0.5.csv")), row.names = FALSE)

  js_genes <- deg %>% filter(p_val < 0.05, avg_log2FC > 0.5) %>% pull(gene)
  ctrl_genes <- deg %>% filter(p_val < 0.05, avg_log2FC < -0.5) %>% pull(gene)
  for (direction in c("JS_high", "ctrl_high")) {
    genes <- if (direction == "JS_high") js_genes else ctrl_genes
    entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (nrow(entrez) > 5) {
      ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
      write.csv(as.data.frame(ego), file.path(out_dir, paste0(safe_ct, "_", direction, "_GO_BP.csv")), row.names = FALSE)
    }
  }
}

wnt_genes <- grep("^WNT", rownames(brain_transfer), value = TRUE)
cilia_genes <- intersect(c("OFD1", "CEP290", "TMEM67", "KIAA0586", "IFT88", "ARL13B", "FOXJ1"), rownames(brain_transfer))
pdf(file.path(out_dir, "WNT_and_cilia_spatial_features.pdf"), width = 10, height = 8)
print(smultiFeaturePlot(brain_transfer, c(wnt_genes, cilia_genes), pt.size = 0.1, ncol = 2))
dev.off()
