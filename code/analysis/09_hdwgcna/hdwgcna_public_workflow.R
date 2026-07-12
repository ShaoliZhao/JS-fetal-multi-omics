# hdWGCNA module analysis curated from code/Rpr0.
#
# Main source scripts:
# - code/Rpr0/cebJS/hdWGCNA/1021/finalcebhdwgcna.R
# - code/Rpr0/kid0326/hdwgcna/hdwgcna.kidney330.R
# - code/Rpr0/kidney/hdwgcna/wgcna0918.R

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
  library(enrichR)
  library(fgsea)
})

theme_set(theme_cowplot())
set.seed(12345)

out_dir <- file.path(paths$results, "hdwgcna")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)

prepare_cerebellum_for_hdwgcna <- function(obj) {
  if (!"name" %in% colnames(obj@meta.data)) obj$name <- obj$orig.ident
  obj$name[obj$name %in% c("pr", "g6", "g7", "g2", "lu", "wxs")] <- "ctrl"
  obj$name[obj$name %in% c("wxj", "wxjba")] <- "CEP290"
  obj$name[obj$name %in% c("zjj")] <- "CPLANE1"
  obj$name[obj$name %in% c("g5v", "g5h")] <- "TMEM67"
  obj$name[obj$name %in% c("lyh", "lyv", "ly1ceb", "KIAA0586")] <- "OFD1"

  obj$celltype[obj$celltype %in% c("ODC", "OPC", "OPCs&ODC")] <- "OPC&ODC"
  obj$celltype[obj$celltype %in% c("MICs")] <- "MG"
  subset(obj, celltype %in% c("Cellcycle", "GCs", "INs", "MG",
                              "OPC&ODC", "PKCs", "UBCs", "VZP"))
}

prepare_kidney_for_hdwgcna <- function(obj) {
  obj$celltype[obj$celltype %in% "Juxtaglomerular"] <- "IC"
  obj$celltype[obj$celltype %in% "Cycling"] <- "NPC"
  subset(obj, celltype %in% c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL",
                              "DCT", "UB_CD", "Stromal", "Endo", "Cycling",
                              "CD", "Prolif", "TAL", "PODO", "IC", "M",
                              "CD-PC", "Cycle", "EC", "LOHP", "SC", "MAC"))
}

run_hdwgcna_all_modules <- function(obj, prefix, min_cells = 50, k = 25,
                                    max_shared = 10, n_hubs = 30) {
  obj <- SetupForWGCNA(
    obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "tutorial"
  )
  obj$cell_type <- obj$celltype
  obj$Sample <- obj$orig.ident
  obj <- RunHarmony(obj, group.by.vars = "Sample")
  obj <- MetacellsByGroups(
    seurat_obj = obj,
    min_cells = min_cells,
    group.by = c("Sample", "cell_type"),
    reduction = "harmony",
    k = k,
    max_shared = max_shared,
    ident.group = "cell_type"
  )
  obj <- NormalizeMetacells(obj)
  obj <- SetDatExpr(obj, group_name = NULL, group.by = NULL, assay = "RNA", layer = "data")
  obj <- TestSoftPowers(obj, networkType = "signed")

  pdf(file.path(out_dir, paste0(prefix, "_soft_power.pdf")), width = 8, height = 6)
  print(wrap_plots(PlotSoftPowers(obj), ncol = 2))
  dev.off()

  obj <- ConstructNetwork(obj, tom_name = "all")
  pdf(file.path(out_dir, paste0(prefix, "_dendrogram.pdf")), width = 10, height = 5)
  PlotDendrogram(obj, main = paste(prefix, "hdWGCNA dendrogram"))
  dev.off()

  obj <- ModuleEigengenes(obj, group.by.vars = "Sample")
  obj <- ModuleConnectivity(obj, group.by = NULL, group_name = NULL)
  obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell")

  modules <- GetModules(obj) %>% subset(module != "grey")
  hub_df <- GetHubGenes(obj, n_hubs = n_hubs)
  write.csv(modules, file.path(out_dir, paste0(prefix, "_modules.csv")), row.names = FALSE)
  write.csv(hub_df, file.path(out_dir, paste0(prefix, "_hub_genes.csv")), row.names = FALSE)

  pdf(file.path(out_dir, paste0(prefix, "_KMEs.pdf")), width = 9, height = 9)
  print(PlotKMEs(obj, ncol = 3))
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_module_feature_plots.pdf")), width = 10, height = 10)
  print(wrap_plots(ModuleFeaturePlot(obj, features = "hMEs", order = TRUE), ncol = 3))
  print(wrap_plots(ModuleFeaturePlot(obj, features = "scores", order = "shuffle", ucell = TRUE), ncol = 3))
  dev.off()

  MEs <- GetMEs(obj, harmonized = TRUE)
  obj@meta.data <- cbind(obj@meta.data, MEs)
  mods <- levels(GetModules(obj)$module)
  mods <- mods[mods != "grey"]

  pdf(file.path(out_dir, paste0(prefix, "_module_dotplots.pdf")), width = 7, height = 8)
  print(DotPlot(obj, features = mods, group.by = "celltype") + RotatedAxis() +
          scale_color_gradient2(high = "red", mid = "grey95", low = "blue"))
  if ("group" %in% colnames(obj@meta.data)) print(DotPlot(obj, features = mods, group.by = "group") + RotatedAxis())
  if ("name" %in% colnames(obj@meta.data)) print(DotPlot(obj, features = mods, group.by = "name") + RotatedAxis())
  if ("geneotype" %in% colnames(obj@meta.data)) print(DotPlot(obj, features = mods, group.by = "geneotype") + RotatedAxis())
  dev.off()

  dme_group_1 <- rownames(obj@meta.data)[obj$group == "JS"]
  dme_group_2 <- rownames(obj@meta.data)[obj$group != "JS"]
  if (length(dme_group_1) > 0 && length(dme_group_2) > 0) {
    DMEs <- FindDMEs(obj, barcodes1 = dme_group_1, barcodes2 = dme_group_2,
                     test.use = "wilcox", wgcna_name = "tutorial")
    write.csv(DMEs, file.path(out_dir, paste0(prefix, "_DMEs_JS_vs_ctrl.csv")), row.names = FALSE)
    pdf(file.path(out_dir, paste0(prefix, "_DMEs_lollipop.pdf")), width = 6, height = 6)
    print(PlotDMEsLollipop(obj, DMEs, wgcna_name = "tutorial", pvalue = "p_val_adj"))
    dev.off()
  }

  obj
}

cerebellum_file <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
if (file.exists(cerebellum_file)) {
  cerebellum <- prepare_cerebellum_for_hdwgcna(readRDS(cerebellum_file))
  cerebellum_hdwgcna <- run_hdwgcna_all_modules(
    cerebellum,
    prefix = "cerebellum_1021_all_celltypes",
    min_cells = 50,
    k = 25,
    max_shared = 10,
    n_hubs = 30
  )
  saveRDS(cerebellum_hdwgcna,
          file.path(out_dir, "cerebellum_1021_all_celltypes_hdwgcna.rds"))
}

kidney_file <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
if (file.exists(kidney_file)) {
  kidney <- prepare_kidney_for_hdwgcna(readRDS(kidney_file))
  kidney_hdwgcna <- run_hdwgcna_all_modules(
    kidney,
    prefix = "kidney_330_or_918_all_celltypes",
    min_cells = 90,
    k = 25,
    max_shared = 10,
    n_hubs = 10
  )
  saveRDS(kidney_hdwgcna,
          file.path(out_dir, "kidney_330_or_918_all_celltypes_hdwgcna.rds"))
}
