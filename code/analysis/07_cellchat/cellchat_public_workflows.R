# CellChat workflows curated from code/Rpr0.
#
# Main source scripts:
# - code/Rpr0/cebJS/cellchat/0423/cc0423.R
# - code/Rpr0/kid0326/cellchat/cckid0330.R
# - code/Rpr0/kidney/cellchat/1208final/cc1208.R
#
# Earlier exploratory scripts under code/Rpr0/cebJS/cellchat/GC and
# code/Rpr0/kidney/cellchat are summarized here rather than copied verbatim.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ComplexHeatmap)
  library(wordcloud)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

out_dir <- file.path(paths$results, "cellchat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_cellchat <- function(obj, group.by = "celltype", trim = 0.01,
                         db_search = NULL, prefix) {
  chat <- createCellChat(obj, group.by = group.by)
  db <- CellChatDB.human
  if (!is.null(db_search)) db <- subsetDB(db, search = db_search)
  chat@DB <- db

  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(chat, type = "truncatedMean", trim = trim)
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  write.csv(subsetCommunication(chat),
            file.path(out_dir, paste0(prefix, "_LR_edges.csv")),
            row.names = FALSE)
  chat <- aggregateNet(chat)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP")
  saveRDS(chat, file.path(out_dir, paste0(prefix, "_CellChat.rds")))
  chat
}

plot_merged_cellchat <- function(cellchat, object.list, prefix, heatmap_height = 24) {
  pdf(file.path(out_dir, paste0(prefix, "_CellChat_comparison.pdf")), width = 12, height = 10)
  print(compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2)) +
          compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight"))
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, edge.width.max = 2)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", edge.width.max = 2)
  print(netVisual_heatmap(cellchat) + netVisual_heatmap(cellchat, measure = "weight"))
  print(rankNet(cellchat, mode = "comparison", measure = "weight",
                stacked = TRUE, do.stat = TRUE))
  print(netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45))

  pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing",
                                           signaling = pathway.union,
                                           title = names(object.list)[1],
                                           width = 5, height = heatmap_height)
  ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing",
                                           signaling = pathway.union,
                                           title = names(object.list)[2],
                                           width = 5, height = heatmap_height)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming",
                                           signaling = pathway.union,
                                           title = names(object.list)[1],
                                           width = 5, height = heatmap_height)
  ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming",
                                           signaling = pathway.union,
                                           title = names(object.list)[2],
                                           width = 5, height = heatmap_height)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
}

map_cerebellum_celltypes <- function(obj) {
  if (!"celltype" %in% colnames(obj@meta.data)) obj$celltype <- as.character(obj$seurat_clusters)
  obj$celltype <- as.character(obj$celltype)
  obj$celltype[obj$seurat_clusters %in% c("9", "12")] <- "UBC"
  obj$celltype[obj$seurat_clusters %in% c("0", "3", "6", "13", "8")] <- "GC"
  obj$celltype[obj$seurat_clusters %in% c("7", "15", "16")] <- "IN"
  obj$celltype[obj$seurat_clusters %in% c("5")] <- "PKC"
  obj$celltype[obj$seurat_clusters %in% c("14")] <- "Cycling"
  obj$celltype[obj$seurat_clusters %in% c("1", "4", "2", "10", "11")] <- "VZP"
  obj$celltype[obj$seurat_clusters %in% c("17")] <- "DN"
  obj$celltype[obj$seurat_clusters %in% c("18")] <- "OPC"
  obj$celltype[obj$seurat_clusters %in% c("21")] <- "ODC"
  obj$celltype[obj$seurat_clusters %in% c("19")] <- "MG"
  obj$celltype[obj$seurat_clusters %in% c("20")] <- "Endo"
  obj$celltype <- factor(obj$celltype,
                         levels = c("VZP", "IN", "PKC", "GC", "UBC", "DN",
                                    "OPC", "ODC", "MG", "Endo", "Cycling"))
  obj
}

# 1. Cerebellum JS/control CellChat, from Rpr0/cebJS/cellchat/0423/cc0423.R.
cerebellum_file <- file.path(paths$data, "controlled/seurat/cerebellum_final_reference.rds")
if (file.exists(cerebellum_file)) {
  cc14 <- map_cerebellum_celltypes(readRDS(cerebellum_file))
  cc14_sub <- subset(cc14, celltype %in% c("VZP", "IN", "PKC", "GC", "UBC", "Cycling"))
  cc14_sub$celltype <- droplevels(cc14_sub$celltype)

  ceb_js <- subset(cc14_sub, group == "JS")
  ceb_ctrl <- subset(cc14_sub, group == "ctrl")
  ceb_js_chat <- run_cellchat(ceb_js, trim = 0.01, prefix = "cerebellum_JS")
  ceb_ctrl_chat <- run_cellchat(ceb_ctrl, trim = 0.01, prefix = "cerebellum_ctrl")
  ceb_merged <- mergeCellChat(list(ctrl = ceb_ctrl_chat, JS = ceb_js_chat),
                              add.names = c("ctrl", "JS"))
  saveRDS(ceb_merged, file.path(out_dir, "cerebellum_ctrl_JS_CellChat_merged.rds"))
  plot_merged_cellchat(ceb_merged, list(ctrl = ceb_ctrl_chat, JS = ceb_js_chat),
                       prefix = "cerebellum_ctrl_JS", heatmap_height = 48)

  pos.dataset <- "JS"
  features.name <- paste0(pos.dataset, ".merged")
  ceb_merged <- identifyOverExpressedGenes(
    ceb_merged,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.05,
    thresh.p = 0.05,
    group.DE.combined = FALSE
  )
  net <- netMappingDEG(ceb_merged, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(ceb_merged, net = net, datasets = "JS",
                                ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(ceb_merged, net = net, datasets = "ctrl",
                                  ligand.logFC = 0.05, receptor.logFC = NULL)
  write.csv(net.up, file.path(out_dir, "cerebellum_CellChat_JS_up_LR.csv"), row.names = FALSE)
  write.csv(net.down, file.path(out_dir, "cerebellum_CellChat_ctrl_up_LR.csv"), row.names = FALSE)
}

# 2. Kidney updated CellChat, from Rpr0/kid0326/cellchat/cckid0330.R.
kidney_file <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
if (file.exists(kidney_file)) {
  kd6rm <- readRDS(kidney_file)
  kidney_cells <- c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT",
                    "UB_CD", "Stromal", "Endo", "Cycling")
  kidney_sub <- subset(kd6rm, celltype %in% kidney_cells)
  kidney_sub$celltype <- droplevels(kidney_sub$celltype)

  kidney_js <- subset(kidney_sub, group == "JS")
  kidney_ctrl <- subset(kidney_sub, group == "ctrl")
  kidney_js_chat <- run_cellchat(kidney_js, trim = 0.01, prefix = "kidney_JS")
  kidney_ctrl_chat <- run_cellchat(kidney_ctrl, trim = 0.01, prefix = "kidney_ctrl")
  kidney_merged <- mergeCellChat(list(ctrl = kidney_ctrl_chat, JS = kidney_js_chat),
                                 add.names = c("ctrl", "JS"))
  saveRDS(kidney_merged, file.path(out_dir, "kidney_ctrl_JS_CellChat_merged.rds"))
  plot_merged_cellchat(kidney_merged, list(ctrl = kidney_ctrl_chat, JS = kidney_js_chat),
                       prefix = "kidney_ctrl_JS", heatmap_height = 48)

  pos.dataset <- "JS"
  features.name <- paste0(pos.dataset, ".merged")
  kidney_merged <- identifyOverExpressedGenes(
    kidney_merged,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.05,
    thresh.p = 0.05,
    group.DE.combined = FALSE
  )
  net <- netMappingDEG(kidney_merged, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(kidney_merged, net = net, datasets = "JS",
                                ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(kidney_merged, net = net, datasets = "ctrl",
                                  ligand.logFC = 0.05, receptor.logFC = NULL)
  write.csv(net.up, file.path(out_dir, "kidney_CellChat_JS_up_LR.csv"), row.names = FALSE)
  write.csv(net.down, file.path(out_dir, "kidney_CellChat_ctrl_up_LR.csv"), row.names = FALSE)

  # Figure-focused SPP1 ligand-receptor bubble used in the kidney CellChat script.
  spp1_pairs <- kidney_merged@DB$interaction[kidney_merged@DB$interaction$ligand == "SPP1",
                                             "interaction_name", drop = FALSE]
  if (nrow(spp1_pairs) > 0) {
    pdf(file.path(out_dir, "kidney_SPP1_UB_CD_to_microenvironment_bubble.pdf"),
        width = 6, height = 4)
    print(netVisual_bubble(
      kidney_merged,
      sources.use = "UB_CD",
      targets.use = c("Podocyte", "NPC", "Stromal", "UB_CD", "DCT",
                      "MAC", "LOH", "Endo", "LOH_DTL", "PT"),
      pairLR.use = spp1_pairs,
      comparison = c(1, 2),
      angle.x = 45,
      remove.isolate = TRUE,
      title.name = "SPP1 signaling from UB/CD"
    ))
    dev.off()
  }
}

# 3. ECM-focused kidney CellChat from Rpr0/kidney/cellchat/1208final/cc1208.R.
# This is kept as an optional focused branch because the later kd0326 workflow
# uses the full CellChatDB. It can be useful for collagen/ECM figure panels.
early_kidney_file <- file.path(paths$data, "controlled/seurat/kidney_early_kd917_or_kd918.rds")
if (file.exists(early_kidney_file)) {
  kd7 <- readRDS(early_kidney_file)
  kd7$celltype[kd7$celltype == "Juxtaglomerular"] <- "IC"
  kd7 <- subset(kd7, celltype %in% c("CD-PC", "Cycle", "DCT", "EC",
                                     "LOHP", "PODO", "PT", "SC", "TAL", "MAC",
                                     "CD", "Prolif", "Endo", "IC", "M", "NPC"))
  kd7$celltype <- droplevels(kd7$celltype)
  early_js <- subset(kd7, group == "JS")
  early_ctrl <- subset(kd7, group == "ctrl")
  early_js_chat <- run_cellchat(early_js, trim = 0.1,
                                db_search = "ECM-Receptor",
                                prefix = "kidney_ECM_JS")
  early_ctrl_chat <- run_cellchat(early_ctrl, trim = 0.1,
                                  db_search = "ECM-Receptor",
                                  prefix = "kidney_ECM_ctrl")
  early_merged <- mergeCellChat(list(JS = early_js_chat, ctrl = early_ctrl_chat),
                                add.names = c("JS", "ctrl"))
  saveRDS(early_merged, file.path(out_dir, "kidney_ECM_JS_ctrl_CellChat_merged.rds"))
}
