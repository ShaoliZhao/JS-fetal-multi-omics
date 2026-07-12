# hdWGCNA module analysis.

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

out_root <- file.path(paths$results, "09_hdwgcna")
fig_root <- file.path(paths$figures, "09_hdwgcna")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  cerebellum = list(
    object = file.path(paths$data_processed, "seurat", "cerebellum_final.rds"),
    celltypes = c("VZP", "IN", "PKC", "GC", "UBC", "Cycling"),
    wgcna_name = "cerebellum",
    fraction = 0.05
  ),
  kidney = list(
    object = file.path(paths$data_processed, "seurat", "kidney_final.rds"),
    celltypes = c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT", "UB_CD", "Stromal", "Endo", "Cycling"),
    wgcna_name = "kidney",
    fraction = 0.05
  )
)

for (dataset in names(datasets)) {
  cfg <- datasets[[dataset]]
  if (!file.exists(cfg$object)) {
    message("Skip ", dataset, ": object not found at ", cfg$object)
    next
  }

  out_dir <- file.path(out_root, dataset)
  fig_dir <- file.path(fig_root, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- readRDS(cfg$object)
  obj <- subset(obj, subset = celltype %in% cfg$celltypes)
  obj$celltype <- droplevels(factor(obj$celltype, levels = cfg$celltypes))
  obj$cell_type <- obj$celltype
  obj$Sample <- obj$orig.ident

  obj <- SetupForWGCNA(
    obj,
    gene_select = "fraction",
    fraction = cfg$fraction,
    wgcna_name = cfg$wgcna_name
  )

  obj <- RunHarmony(obj, group.by.vars = "Sample")
  obj <- MetacellsByGroups(
    seurat_obj = obj,
    group.by = c("Sample", "cell_type"),
    reduction = "harmony",
    ident.group = "cell_type"
  )
  obj <- NormalizeMetacells(obj)

  obj <- SetDatExpr(
    obj,
    group_name = NULL,
    group.by = NULL,
    assay = "RNA",
    layer = "data"
  )

  obj <- TestSoftPowers(obj, networkType = "signed")
  power_table <- GetPowerTable(obj)
  write.csv(power_table, file.path(out_dir, "soft_power_table.csv"), row.names = FALSE)

  p <- wrap_plots(PlotSoftPowers(obj), ncol = 2)
  ggsave(file.path(fig_dir, "soft_power_diagnostics.pdf"), p, width = 8, height = 6)

  obj <- ConstructNetwork(obj, tom_name = dataset)
  pdf(file.path(fig_dir, "module_dendrogram.pdf"), width = 10, height = 5)
  PlotDendrogram(obj, main = paste(dataset, "hdWGCNA dendrogram"))
  dev.off()

  obj <- ModuleEigengenes(obj, group.by.vars = "Sample")
  obj <- ModuleConnectivity(obj, group.by = NULL, group_name = NULL)

  modules <- GetModules(obj) %>% filter(module != "grey")
  hub_df <- GetHubGenes(obj, n_hubs = 10)
  write.csv(modules, file.path(out_dir, "module_assignments.csv"), row.names = FALSE)
  write.csv(hub_df, file.path(out_dir, "top10_hub_genes.csv"), row.names = FALSE)

  p <- PlotKMEs(obj, ncol = 3)
  ggsave(file.path(fig_dir, "module_kME_rankings.pdf"), p, width = 9, height = 9)

  obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell")

  p <- wrap_plots(ModuleFeaturePlot(obj, features = "hMEs", order = TRUE), ncol = 4)
  ggsave(file.path(fig_dir, "module_eigengene_featureplots.pdf"), p, width = 10, height = 12)

  p <- wrap_plots(ModuleFeaturePlot(obj, features = "scores", order = "shuffle", ucell = TRUE), ncol = 4)
  ggsave(file.path(fig_dir, "hub_gene_score_featureplots.pdf"), p, width = 10, height = 12)

  pdf(file.path(fig_dir, "module_correlogram.pdf"), width = 8, height = 8)
  ModuleCorrelogram(obj)
  dev.off()

  MEs <- GetMEs(obj, harmonized = TRUE)
  mods <- levels(GetModules(obj)$module)
  mods <- mods[mods != "grey"]
  obj@meta.data <- cbind(obj@meta.data, MEs[colnames(obj), mods, drop = FALSE])

  p <- DotPlot(obj, features = mods, group.by = "celltype") +
    RotatedAxis() +
    scale_color_gradient2(high = "#E25857", mid = "grey95", low = "#6EB65A") +
    theme_js_atlas(8)
  ggsave(file.path(fig_dir, "module_dotplot_celltype.pdf"), p, width = 7, height = 4)

  if ("group" %in% colnames(obj@meta.data)) {
    p <- DotPlot(obj, features = mods, group.by = "group") +
      RotatedAxis() +
      scale_color_gradient2(high = "#E25857", mid = "grey95", low = "#3A6EA5") +
      theme_js_atlas(8)
    ggsave(file.path(fig_dir, "module_dotplot_group.pdf"), p, width = 6, height = 4)
  }

  if ("geneotype" %in% colnames(obj@meta.data)) {
    p <- DotPlot(obj, features = mods, group.by = "geneotype") +
      RotatedAxis() +
      scale_color_gradient2(high = "#E25857", mid = "grey95", low = "#3A6EA5") +
      theme_js_atlas(8)
    ggsave(file.path(fig_dir, "module_dotplot_genotype.pdf"), p, width = 6, height = 4)
  }

  saveRDS(obj, file.path(out_dir, paste0(dataset, "_hdWGCNA.rds")))
}
