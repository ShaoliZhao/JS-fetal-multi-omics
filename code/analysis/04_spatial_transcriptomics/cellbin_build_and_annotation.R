# Build and annotate JS/control cerebellar cell-bin Stereo-seq objects.
#
# Source provenance: code/ST_code_backup/Rpro/ana/cellbin/trycellbin.R.
# The script keeps the manuscript workflow and thresholds, but replaces private
# server paths with local variables. Supply the two Stereo-seq cell-bin Seurat
# objects below before running.

source("analysis/00_setup/project_config.R")
source("analysis/04_spatial_transcriptomics/STimport_public.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
})

ofd1_rds <- file.path(paths$data, "controlled/spatial/OFD1_cellbin_raw.rds")
ctrl_rds <- file.path(paths$data, "controlled/spatial/control_w18_cellbin_raw.rds")
out_dir <- file.path(paths$results, "spatial_cellbin")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ofd1_raw <- readRDS(ofd1_rds)
ctrl_raw <- readRDS(ctrl_rds)

build_cellbin_object <- function(raw_obj, project) {
  raw_counts <- GetAssayData(raw_obj, assay = "Spatial", slot = "counts")
  coords <- data.frame(Spatial_1 = raw_obj$x, Spatial_2 = raw_obj$y)
  rownames(coords) <- colnames(raw_counts)

  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(raw_counts),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  keep_ensembl <- is.na(gene_symbols) | duplicated(gene_symbols)
  gene_symbols[keep_ensembl] <- rownames(raw_counts)[keep_ensembl]
  rownames(raw_counts) <- gene_symbols

  obj <- CreateSeuratObject(counts = raw_counts, min.features = 0,
                            assay = "Spatial", project = project)
  obj <- AddMetaData(obj, coords)
  obj[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(coords[colnames(obj), ]),
    key = "Spatial_",
    assay = "Spatial"
  )

  DefaultAssay(obj) <- "Spatial"
  obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj)
  obj <- RunUMAP(obj, dims = 1:30)
  obj
}

ofd1_cellbin <- build_cellbin_object(ofd1_raw, "OFD1")
ctrl_cellbin <- build_cellbin_object(ctrl_raw, "ctrl")
saveRDS(ofd1_cellbin, file.path(out_dir, "OFD1_cellbin_sct.rds"))
saveRDS(ctrl_cellbin, file.path(out_dir, "ctrl_cellbin_sct.rds"))

pdf(file.path(out_dir, "cellbin_clusters_before_merge.pdf"), width = 10, height = 5)
print(sDimplot(ofd1_cellbin, my25cols, "seurat_clusters") + sDimplot(ctrl_cellbin, my25cols, "seurat_clusters"))
dev.off()

brain_merge <- merge(ctrl_cellbin, ofd1_cellbin)
DefaultAssay(brain_merge) <- "SCT"
VariableFeatures(brain_merge) <- unique(c(VariableFeatures(ctrl_cellbin), VariableFeatures(ofd1_cellbin)))
brain_merge <- RunPCA(brain_merge, verbose = FALSE)
brain_merge <- FindNeighbors(brain_merge, dims = 1:30)
brain_merge <- FindClusters(brain_merge)
brain_merge <- RunUMAP(brain_merge, dims = 1:30)
brain_merge <- PrepSCTFindMarkers(brain_merge)

Idents(brain_merge) <- "seurat_clusters"
markers <- FindAllMarkers(brain_merge, only.pos = TRUE, logfc.threshold = 0.5)
markers_top30 <- markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
write.csv(markers, file.path(out_dir, "cellbin_cluster_markers_logFC0.5.csv"), row.names = FALSE)
write.csv(markers_top30, file.path(out_dir, "cellbin_cluster_markers_top30.csv"), row.names = FALSE)

brain_merge$celltype <- as.character(brain_merge$seurat_clusters)
brain_merge$celltype[brain_merge$celltype %in% c("13")] <- "RL"
brain_merge$celltype[brain_merge$celltype %in% c("0", "3", "7", "9", "20", "21", "11")] <- "iGL"
brain_merge$celltype[brain_merge$celltype %in% c("22")] <- "PKCm"
brain_merge$celltype[brain_merge$celltype %in% c("6")] <- "PKCi"
brain_merge$celltype[brain_merge$celltype %in% c("5")] <- "brainstem"
brain_merge$celltype[brain_merge$celltype %in% c("19")] <- "ODC"
brain_merge$celltype[brain_merge$celltype %in% c("18")] <- "MG"
brain_merge$celltype[brain_merge$celltype %in% c("8")] <- "Ependymal"
brain_merge$celltype[brain_merge$celltype %in% c("14", "17")] <- "Endo"
brain_merge$celltype[brain_merge$celltype %in% c("10")] <- "ERY"
brain_merge$celltype[brain_merge$celltype %in% c("4")] <- "eGL"
brain_merge$celltype[brain_merge$celltype %in% c("16")] <- "Per"
brain_merge$celltype[brain_merge$celltype %in% c("1")] <- "BG/AST"
brain_merge$celltype[brain_merge$celltype %in% c("2")] <- "STR"
brain_merge$celltype[brain_merge$celltype %in% c("12")] <- "Cycling"
brain_merge$celltype[brain_merge$celltype %in% c("15")] <- "CFAP298_hi"

# Clusters 24-26 were excluded in the original analysis because they did not
# represent interpretable cerebellar cell-bin domains.
brain_merge <- subset(brain_merge, seurat_clusters %in% as.character(0:23))

pdf(file.path(out_dir, "cellbin_celltype_annotation.pdf"), width = 10, height = 5)
print(smultiDimplot(brain_merge, my25cols, group = "celltype", ncol = 2, pt.size = 0.1))
dev.off()

saveRDS(brain_merge, file.path(out_dir, "cellbin_brain_merge_annotated.rds"))

# Purkinje layer comparison used for inner/outer PKC interpretation.
Idents(brain_merge) <- "seurat_clusters"
pkc_deg <- FindMarkers(
  brain_merge,
  ident.1 = "6",
  ident.2 = "22",
  logfc.threshold = 0.3,
  min.diff.pct = 0.1,
  min.pct = 0
)
pkc_deg$gene <- rownames(pkc_deg)
pkc_marker <- FindMarkers(brain_merge, ident.1 = c("6", "22"), only.pos = TRUE, logfc.threshold = 1)
pkc_marker$gene <- rownames(pkc_marker)
pkc_core <- intersect(pkc_deg$gene, pkc_marker$gene)
pkc_js <- pkc_deg %>% filter(gene %in% pkc_core, p_val < 0.05, avg_log2FC > 0.3)
pkc_ctrl <- pkc_deg %>% filter(gene %in% pkc_core, p_val < 0.05, avg_log2FC < -0.3)
write.csv(pkc_deg, file.path(out_dir, "PKCi_vs_PKCm_FindMarkers_logFC0.3.csv"), row.names = FALSE)

for (nm in c("pkc_js", "pkc_ctrl")) {
  genes <- get(nm)$gene
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez) > 5) {
    ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
    write.csv(as.data.frame(ego), file.path(out_dir, paste0(nm, "_GO_BP.csv")), row.names = FALSE)
  }
}
