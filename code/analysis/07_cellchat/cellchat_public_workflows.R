# CellChat JS/control comparison workflows.

source("analysis/00_setup/packages.R")
source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
})

out_root <- file.path(paths$results, "07_cellchat")
fig_root <- file.path(paths$figures, "07_cellchat")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  cerebellum = list(
    object = file.path(paths$data_processed, "seurat", "cerebellum_final.rds"),
    celltypes = c("VZP", "IN", "PKC", "GC", "UBC", "Cycling"),
    group_col = "group",
    celltype_col = "celltype",
    trim = 0.01
  ),
  kidney = list(
    object = file.path(paths$data_processed, "seurat", "kidney_final.rds"),
    celltypes = c("NPC", "Podocyte", "PT", "LOH", "PEC", "DCT", "UB_CD", "Stromal", "Endo", "Cycling"),
    group_col = "group",
    celltype_col = "celltype",
    trim = 0.01
  )
)

CellChatDB <- CellChatDB.human

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
  obj[[cfg$celltype_col]][, 1] <- droplevels(factor(obj[[cfg$celltype_col]][, 1], levels = cfg$celltypes))
  DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)

  JS <- subset(obj, subset = group == "JS")
  ctrl <- subset(obj, subset = group == "ctrl")
  object.list <- list()

  for (group_name in c("ctrl", "JS")) {
    sub <- if (group_name == "ctrl") ctrl else JS
    cellchat <- createCellChat(sub, group.by = cfg$celltype_col)
    cellchat@DB <- CellChatDB
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = cfg$trim)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    df.net <- subsetCommunication(cellchat)
    write.csv(df.net, file.path(out_dir, paste0(group_name, "_communication_table.csv")), row.names = FALSE)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, file.path(out_dir, paste0(group_name, "_cellchat.rds")))
    object.list[[group_name]] <- cellchat
  }

  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  saveRDS(cellchat, file.path(out_dir, "ctrl_JS_merged_cellchat.rds"))

  pdf(file.path(fig_dir, "compare_interaction_number_and_weight.pdf"), width = 7, height = 4)
  p1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
  p2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
  print(p1 + p2)
  dev.off()

  pdf(file.path(fig_dir, "differential_network_circle.pdf"), width = 8, height = 4)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, edge.width.max = 2)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", edge.width.max = 2)
  dev.off()

  pdf(file.path(fig_dir, "merged_heatmap_count_weight.pdf"), width = 8, height = 4)
  h1 <- netVisual_heatmap(cellchat)
  h2 <- netVisual_heatmap(cellchat, measure = "weight")
  print(h1 + h2)
  dev.off()

  pdf(file.path(fig_dir, "rankNet_and_signaling_role_heatmaps.pdf"), width = 12, height = 10)
  print(rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE))
  print(rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE))

  pathway.union <- union(object.list$ctrl@netP$pathways, object.list$JS@netP$pathways)
  ht1 <- netAnalysis_signalingRole_heatmap(object.list$ctrl, pattern = "outgoing", signaling = pathway.union, title = "ctrl")
  ht2 <- netAnalysis_signalingRole_heatmap(object.list$JS, pattern = "outgoing", signaling = pathway.union, title = "JS")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

  ht1 <- netAnalysis_signalingRole_heatmap(object.list$ctrl, pattern = "incoming", signaling = pathway.union, title = "ctrl")
  ht2 <- netAnalysis_signalingRole_heatmap(object.list$JS, pattern = "incoming", signaling = pathway.union, title = "JS")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()

  pdf(file.path(fig_dir, "bubble_all_and_dysregulated_LR_pairs.pdf"), width = 14, height = 10)
  print(netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE))
  print(netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, angle.x = 45, remove.isolate = TRUE, title.name = "Higher in JS"))
  print(netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, angle.x = 45, remove.isolate = TRUE, title.name = "Higher in ctrl"))

  pos.dataset <- "JS"
  features.name <- paste0(pos.dataset, ".merged")
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.05,
    thresh.p = 0.05,
    group.DE.combined = FALSE
  )
  net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(cellchat, net = net, datasets = "JS", ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(cellchat, net = net, datasets = "ctrl", ligand.logFC = 0.05, receptor.logFC = NULL)
  write.csv(net.up, file.path(out_dir, "JS_higher_ligand_receptor_pairs.csv"), row.names = FALSE)
  write.csv(net.down, file.path(out_dir, "ctrl_higher_ligand_receptor_pairs.csv"), row.names = FALSE)

  if (nrow(net.up) > 0) {
    print(netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = FALSE], comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE))
  }
  if (nrow(net.down) > 0) {
    print(netVisual_bubble(cellchat, pairLR.use = net.down[, "interaction_name", drop = FALSE], comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE))
  }
  dev.off()

  saveRDS(cellchat, file.path(out_dir, "ctrl_JS_merged_cellchat_with_DEG_mapping.rds"))
}
