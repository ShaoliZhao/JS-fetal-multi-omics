# hdWGCNA module analysis following the single-cell and DME tutorials.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(hdWGCNA)
  library(WGCNA)
  library(UCell)
  library(cowplot)
  library(harmony)
  library(future)
})

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 64 * 1024^3)

analysis_name <- "kidney_nephron"
object_file <- file.path(paths$data_processed, "seurat", "kidney_final.rds")

celltypes_use <- c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT", "UB_CD", "Stromal", "Endo", "Cycling")
network_celltypes <- celltypes_use
celltype_col <- "celltype"
sample_col <- "orig.ident"
group_col <- "group"
wgcna_name <- "tutorial"

out_dir <- file.path(paths$results, "09_hdwgcna", analysis_name)
fig_dir <- file.path(paths$figures, "09_hdwgcna", analysis_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS(object_file)
seurat_obj <- subset(seurat_obj, subset = celltype %in% celltypes_use)
seurat_obj$cell_type <- droplevels(factor(seurat_obj[[celltype_col]][, 1], levels = celltypes_use))
seurat_obj$Sample <- seurat_obj[[sample_col]][, 1]

p <- DimPlot(seurat_obj, group.by = "cell_type", label = TRUE) +
  ggtitle(analysis_name) +
  NoLegend()
ggsave(file.path(fig_dir, "input_celltypes_umap.pdf"), p, width = 6, height = 5)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = wgcna_name
)

if (!"harmony" %in% names(seurat_obj@reductions)) {
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Sample")
}

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_type", "Sample"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "cell_type"
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- ScaleMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = "Sample")
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction = "harmony", dims = 1:15)

p1 <- DimPlotMetacells(seurat_obj, group.by = "cell_type") + umap_theme() + ggtitle("Cell type")
p2 <- DimPlotMetacells(seurat_obj, group.by = "Sample") + umap_theme() + ggtitle("Sample")
ggsave(file.path(fig_dir, "metacell_umap.pdf"), p1 | p2, width = 9, height = 4)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = NULL,
  group.by = NULL,
  assay = "RNA",
  layer = "data"
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed"
)

plot_list <- PlotSoftPowers(seurat_obj)
ggsave(file.path(fig_dir, "soft_power_diagnostics.pdf"), wrap_plots(plot_list, ncol = 2), width = 8, height = 6)

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(out_dir, "soft_power_table.csv"), row.names = FALSE)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = analysis_name
)

pdf(file.path(fig_dir, "module_dendrogram.pdf"), width = 10, height = 5)
PlotDendrogram(seurat_obj, main = paste0(analysis_name, " hdWGCNA dendrogram"))
dev.off()

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "Sample"
)

hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized = FALSE)
write.csv(hMEs, file.path(out_dir, "harmonized_module_eigengenes.csv"))
write.csv(MEs, file.path(out_dir, "module_eigengenes.csv"))

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = NULL,
  group_name = NULL
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(analysis_name, "-M")
)

p <- PlotKMEs(seurat_obj, ncol = 5)
ggsave(file.path(fig_dir, "kME_rankings.pdf"), p, width = 10, height = 8)

modules <- GetModules(seurat_obj) %>% subset(module != "grey")
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
write.csv(modules, file.path(out_dir, "module_assignments.csv"), row.names = FALSE)
write.csv(hub_df, file.path(out_dir, "top10_hub_genes.csv"), row.names = FALSE)

saveRDS(seurat_obj, file.path(out_dir, paste0(analysis_name, "_hdWGCNA_core.rds")))

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method = "UCell"
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = "hMEs",
  order = TRUE
)
ggsave(file.path(fig_dir, "module_eigengene_featureplots.pdf"), wrap_plots(plot_list, ncol = 6), width = 12, height = 10)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = "scores",
  order = "shuffle",
  ucell = TRUE
)
ggsave(file.path(fig_dir, "hub_score_featureplots.pdf"), wrap_plots(plot_list, ncol = 6), width = 12, height = 10)

pdf(file.path(fig_dir, "module_correlogram.pdf"), width = 8, height = 8)
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != "grey"]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs[colnames(seurat_obj), mods, drop = FALSE])

p <- DotPlot(seurat_obj, features = mods, group.by = "cell_type") +
  RotatedAxis() +
  scale_color_gradient2(high = "#E25857", mid = "grey95", low = "#3A6EA5")
ggsave(file.path(fig_dir, "module_dotplot_celltype.pdf"), p, width = 8, height = 4)

p <- DotPlot(seurat_obj, features = mods, group.by = group_col) +
  RotatedAxis() +
  scale_color_gradient2(high = "#E25857", mid = "grey95", low = "#3A6EA5")
ggsave(file.path(fig_dir, "module_dotplot_JS_ctrl.pdf"), p, width = 5, height = 4)

meta <- seurat_obj@meta.data
group1 <- rownames(meta)[meta$cell_type %in% network_celltypes & meta[[group_col]] == "JS"]
group2 <- rownames(meta)[meta$cell_type %in% network_celltypes & meta[[group_col]] == "ctrl"]

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  wgcna_name = wgcna_name
)
write.csv(DMEs, file.path(out_dir, "DME_JS_vs_ctrl.csv"))

p <- PlotDMEsLollipop(
  seurat_obj,
  DMEs,
  wgcna_name = wgcna_name,
  pvalue = "p_val_adj"
)
ggsave(file.path(fig_dir, "DME_JS_vs_ctrl_lollipop.pdf"), p, width = 7, height = 5)

p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = wgcna_name
)
ggsave(file.path(fig_dir, "DME_JS_vs_ctrl_volcano.pdf"), p, width = 5, height = 4)

saveRDS(seurat_obj, file.path(out_dir, paste0(analysis_name, "_hdWGCNA_with_DME.rds")))
