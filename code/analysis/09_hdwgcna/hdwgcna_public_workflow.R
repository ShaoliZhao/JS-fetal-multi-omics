# hdWGCNA workflows for cerebellum and kidney co-expression modules.
#
# Source scripts:
# - code/Rpro/hdWGCNA/ceball/GC/allGC_hdwgcna.R
# - code/Rpro/hdWGCNA/ceball/VZP/allVZP_hdwgcna.R
# - code/Rpro/hdWGCNA/ceball/pkc/allpkc_hdwgcna.R
# - code/Rpr0/kid0326/hdwgcna/hdwgcna.kidney330.R

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
  library(enrichR)
})

theme_set(cowplot::theme_cowplot())
set.seed(12345)

out_dir <- file.path(paths$results, "hdwgcna")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)

prepare_cerebellum_metadata <- function(obj) {
  obj$group <- obj$orig.ident
  obj$group[obj$group %in% c("pr", "g6", "g7", "g2", "lu", "wxs")] <- "ctrl"
  obj$group[obj$group %in% c("wxj", "wxjba", "zjj", "lyh", "lyv", "g5v", "g5h", "ly1ceb")] <- "JS"

  obj$celltype <- as.character(obj$seurat_clusters)
  obj$celltype[obj$celltype %in% c("11", "6")] <- "UBCs"
  obj$celltype[obj$celltype %in% c("1", "4", "7", "13", "8")] <- "GCs"
  obj$celltype[obj$celltype %in% c("9", "12", "15", "17")] <- "INs"
  obj$celltype[obj$celltype %in% c("5")] <- "PKCs"
  obj$celltype[obj$celltype %in% c("14")] <- "Cellcycle"
  obj$celltype[obj$celltype %in% c("0", "2", "3", "10", "16")] <- "VZP"

  obj$name <- obj$orig.ident
  obj$name[obj$name %in% c("pr", "g6", "g7", "g2", "lu", "wxs")] <- "ctrl"
  obj$name[obj$name %in% c("wxj", "wxjba")] <- "CEP290"
  obj$name[obj$name %in% c("zjj")] <- "CPLANE1"
  obj$name[obj$name %in% c("g5v", "g5h")] <- "TMEM67"
  obj$name[obj$name %in% c("lyh", "lyv", "ly1ceb")] <- "KIAA0586"
  obj
}

run_hdwgcna_core <- function(obj, group_name, prefix, group.by = "cell_type",
                             metacell_k = 50, max_shared = 10,
                             network_tom_name = group_name) {
  obj <- SetupForWGCNA(obj, gene_select = "fraction", fraction = 0.05, wgcna_name = "tutorial")
  obj$cell_type <- obj$celltype
  obj$Sample <- obj$orig.ident
  obj <- RunHarmony(obj, group.by.vars = "Sample")
  obj <- MetacellsByGroups(
    seurat_obj = obj,
    group.by = c("cell_type", "Sample"),
    reduction = "harmony",
    k = metacell_k,
    max_shared = max_shared,
    ident.group = "cell_type"
  )
  obj <- NormalizeMetacells(obj)
  obj <- SetDatExpr(obj, group_name = group_name, group.by = group.by, assay = "RNA", layer = "data")
  obj <- TestSoftPowers(obj, networkType = "signed")

  pdf(file.path(out_dir, paste0(prefix, "_soft_power.pdf")), width = 8, height = 6)
  print(wrap_plots(PlotSoftPowers(obj), ncol = 2))
  dev.off()

  obj <- ConstructNetwork(obj, tom_name = network_tom_name)
  pdf(file.path(out_dir, paste0(prefix, "_dendrogram.pdf")), width = 10, height = 5)
  PlotDendrogram(obj, main = paste(prefix, "hdWGCNA dendrogram"))
  dev.off()

  obj <- ModuleEigengenes(obj, group.by.vars = "Sample")
  obj <- ModuleConnectivity(obj, group.by = group.by, group_name = group_name)
  obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell")

  modules <- GetModules(obj) %>% subset(module != "grey")
  hub_df <- GetHubGenes(obj, n_hubs = 100)
  write.csv(modules, file.path(out_dir, paste0(prefix, "_modules.csv")), row.names = FALSE)
  write.csv(hub_df, file.path(out_dir, paste0(prefix, "_hub_genes.csv")), row.names = FALSE)

  pdf(file.path(out_dir, paste0(prefix, "_KMEs.pdf")), width = 12, height = 9)
  print(PlotKMEs(obj, ncol = 3))
  dev.off()

  pdf(file.path(out_dir, paste0(prefix, "_module_feature_plots.pdf")), width = 10, height = 8)
  print(wrap_plots(ModuleFeaturePlot(obj, features = "hMEs", order = TRUE), ncol = 5))
  print(wrap_plots(ModuleFeaturePlot(obj, features = "scores", order = "shuffle", ucell = TRUE), ncol = 5))
  dev.off()

  obj
}

cerebellum_file <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
if (file.exists(cerebellum_file)) {
  cerebellum <- prepare_cerebellum_metadata(readRDS(cerebellum_file))
  for (ct in c("GCs", "VZP", "PKCs")) {
    ct_obj <- run_hdwgcna_core(cerebellum, group_name = ct, prefix = paste0("cerebellum_", ct))
    saveRDS(ct_obj, file.path(out_dir, paste0("cerebellum_", ct, "_hdWGCNA.rds")))
  }
}

kidney_file <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
if (file.exists(kidney_file)) {
  kidney <- readRDS(kidney_file)
  kidney_cells <- c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT", "UB_CD", "Stromal", "Endo", "Cycling")
  kidney <- subset(kidney, celltype %in% kidney_cells)
  kidney$celltype[kidney$celltype == "Cycling"] <- "NPC"
  kidney$celltype <- droplevels(kidney$celltype)

  kidney <- SetupForWGCNA(kidney, gene_select = "fraction", fraction = 0.05, wgcna_name = "tutorial")
  kidney$cell_type <- kidney$celltype
  kidney$Sample <- kidney$orig.ident
  kidney <- RunHarmony(kidney, group.by.vars = "Sample")
  kidney <- MetacellsByGroups(
    seurat_obj = kidney,
    group.by = c("Sample", "cell_type"),
    reduction = "harmony",
    ident.group = "cell_type"
  )
  kidney <- NormalizeMetacells(kidney)
  kidney <- SetDatExpr(kidney, group_name = NULL, group.by = NULL, assay = "RNA", layer = "data")
  kidney <- TestSoftPowers(kidney, networkType = "signed")
  kidney <- ConstructNetwork(kidney, tom_name = "all")
  kidney <- ModuleEigengenes(kidney, group.by.vars = "Sample")
  kidney <- ModuleConnectivity(kidney, group.by = NULL, group_name = NULL)
  kidney <- ModuleExprScore(kidney, n_genes = 25, method = "UCell")

  modules <- GetModules(kidney) %>% subset(module != "grey")
  write.csv(modules, file.path(out_dir, "kidney_modules.csv"), row.names = FALSE)
  write.csv(GetHubGenes(kidney, n_hubs = 100), file.path(out_dir, "kidney_hub_genes.csv"), row.names = FALSE)

  MEs <- GetMEs(kidney, harmonized = TRUE)
  kidney@meta.data <- cbind(kidney@meta.data, MEs)
  mods <- levels(GetModules(kidney)$module)
  mods <- mods[mods != "grey"]

  pdf(file.path(out_dir, "kidney_module_dotplots.pdf"), width = 7, height = 8)
  print(DotPlot(kidney, features = mods, group.by = "celltype") + RotatedAxis() +
          scale_color_gradient2(high = "#E25857", mid = "lightgrey", low = "#6EB65A"))
  if ("group" %in% colnames(kidney@meta.data)) print(DotPlot(kidney, features = mods, group.by = "group") + RotatedAxis())
  if ("geneotype" %in% colnames(kidney@meta.data)) print(DotPlot(kidney, features = mods, group.by = "geneotype") + RotatedAxis())
  dev.off()

  saveRDS(kidney, file.path(out_dir, "kidney_hdwgcna.rds"))
}
