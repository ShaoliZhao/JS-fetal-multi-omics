# CellChat workflows kept as separate analysis modules.
#
# Source scripts:
# - code/Rpr0/kid0326/cellchat/cckid0330.R
# - code/ST_code_backup/Rpro/ana/cellbin/0127/cellchat/cebcc.st.R

source("analysis/00_setup/project_config.R")
source("analysis/04_spatial_transcriptomics/STimport_public.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(FNN)
})

out_dir <- file.path(paths$results, "cellchat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)
set.seed(12345)

run_cellchat_singlecell <- function(obj, group.by = "celltype", trim = 0.01) {
  chat <- createCellChat(obj, group.by = group.by)
  chat@DB <- CellChatDB.human
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(chat, type = "truncatedMean", trim = trim)
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP")
  chat
}

# Kidney single-nucleus CellChat, from cckid0330.R.
kidney_obj_file <- file.path(paths$data, "controlled/seurat/kidney_final.rds")
if (file.exists(kidney_obj_file)) {
  kd6rm <- readRDS(kidney_obj_file)
  kidney_cells <- c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT",
                    "UB_CD", "Stromal", "Endo", "Cycling")
  kidney_sub <- subset(kd6rm, celltype %in% kidney_cells)
  kidney_sub$celltype <- droplevels(kidney_sub$celltype)
  kidney_sub$celltype[kidney_sub$celltype == "Cycling"] <- "NPC"

  kidney_js <- subset(kidney_sub, group == "JS")
  kidney_ctrl <- subset(kidney_sub, group == "ctrl")

  kidney_js_chat <- run_cellchat_singlecell(kidney_js, trim = 0.01)
  kidney_ctrl_chat <- run_cellchat_singlecell(kidney_ctrl, trim = 0.01)
  kidney_merged <- mergeCellChat(list(ctrl = kidney_ctrl_chat, JS = kidney_js_chat),
                                 add.names = c("ctrl", "JS"))

  write.csv(subsetCommunication(kidney_js_chat), file.path(out_dir, "kidney_JS_CellChat_edges.csv"))
  write.csv(subsetCommunication(kidney_ctrl_chat), file.path(out_dir, "kidney_ctrl_CellChat_edges.csv"))
  saveRDS(kidney_js_chat, file.path(out_dir, "kidney_JS_CellChat.rds"))
  saveRDS(kidney_ctrl_chat, file.path(out_dir, "kidney_ctrl_CellChat.rds"))
  saveRDS(kidney_merged, file.path(out_dir, "kidney_ctrl_JS_CellChat_merged.rds"))

  pdf(file.path(out_dir, "kidney_CellChat_comparison.pdf"), width = 8, height = 8)
  print(compareInteractions(kidney_merged, show.legend = FALSE, group = c(1, 2)) +
          compareInteractions(kidney_merged, show.legend = FALSE, group = c(1, 2), measure = "weight"))
  netVisual_diffInteraction(kidney_merged, weight.scale = TRUE, edge.width.max = 2)
  netVisual_diffInteraction(kidney_merged, weight.scale = TRUE, measure = "weight", edge.width.max = 2)
  print(rankNet(kidney_merged, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE))
  print(netVisual_bubble(kidney_merged, comparison = c(1, 2), angle.x = 45))
  dev.off()
}

run_cellchat_spatial <- function(obj, sample_name, trim) {
  data.input <- GetAssayData(obj, assay = "SCT", slot = "data")
  meta <- obj@meta.data
  coordinates <- data.frame(x = meta$Spatial_1, y = meta$Spatial_2)
  rownames(coordinates) <- colnames(obj)

  nn <- FNN::get.knn(coordinates, k = 2)
  median_nn_um <- median(nn$nn.dist[, 2]) * 0.5
  message(sample_name, " median nearest-neighbor distance: ", round(median_nn_um, 3), " um")

  keep <- meta$celltype != "low"
  chat <- createCellChat(
    object = data.input[, keep],
    meta = meta[keep, ],
    group.by = "celltype",
    datatype = "spatial",
    coordinates = coordinates[keep, ],
    spatial.factors = data.frame(ratio = 0.5, tol = 10, row.names = sample_name)
  )
  chat@DB <- CellChatDB.human
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(
    chat,
    type = "truncatedMean",
    trim = trim,
    distance.use = TRUE,
    interaction.range = 20,
    scale.distance = 20,
    contact.dependent = TRUE,
    contact.range = 10
  )
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  chat <- netAnalysis_computeCentrality(chat, slot.name = "netP")
  chat
}

# Spatial cell-bin CellChat, from cebcc.st.R.
spatial_obj_file <- file.path(paths$results, "spatial_cellbin/cellbin_brain_merge_annotated.rds")
if (file.exists(spatial_obj_file)) {
  brain_merge <- readRDS(spatial_obj_file)
  ofd1 <- subset(brain_merge, orig.ident == "OFD1")
  ctrl <- subset(brain_merge, orig.ident == "ctrl")

  ofd1_chat <- run_cellchat_spatial(ofd1, "OFD1", trim = 0.1)
  ctrl_chat <- run_cellchat_spatial(ctrl, "ctrl", trim = 0.5)
  spatial_merged <- mergeCellChat(list(JS = ofd1_chat, ctrl = ctrl_chat),
                                  add.names = c("JS", "ctrl"),
                                  cell.prefix = TRUE)

  saveRDS(ofd1_chat, file.path(out_dir, "spatial_OFD1_CellChat.rds"))
  saveRDS(ctrl_chat, file.path(out_dir, "spatial_ctrl_CellChat.rds"))
  saveRDS(spatial_merged, file.path(out_dir, "spatial_JS_ctrl_CellChat_merged.rds"))

  pos.dataset <- "JS"
  features.name <- paste0(pos.dataset, ".merged")
  spatial_merged <- identifyOverExpressedGenes(
    spatial_merged,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.05,
    thresh.p = 0.05,
    group.DE.combined = FALSE
  )
  net <- netMappingDEG(spatial_merged, features.name = features.name, variable.all = TRUE)
  net_up <- subsetCommunication(spatial_merged, net = net, datasets = "JS", ligand.logFC = 0.05)
  net_down <- subsetCommunication(spatial_merged, net = net, datasets = "ctrl", ligand.logFC = -0.05)
  write.csv(net_up, file.path(out_dir, "spatial_CellChat_JS_up_LR.csv"), row.names = FALSE)
  write.csv(net_down, file.path(out_dir, "spatial_CellChat_ctrl_up_LR.csv"), row.names = FALSE)
}
