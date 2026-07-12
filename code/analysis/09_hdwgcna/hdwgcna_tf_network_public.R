# hdWGCNA TF-network / regulon analysis curated from code/Rpr0.
#
# Main source scripts:
# - code/Rpr0/cebJS/hdWGCNA/1021/TF/TFceb1021.R
# - code/Rpr0/kidney/hdwgcna/TF/TFhdWGCNAkidney.R

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(JASPAR2020)
  library(motifmatchr)
  library(TFBSTools)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(xgboost)
  library(enrichR)
})

theme_set(theme_cowplot())
set.seed(12345)

out_dir <- file.path(paths$results, "hdwgcna_tf_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_hdwgcna_tf_network <- function(hdwgcna_file, prefix) {
  seurat_obj <- readRDS(hdwgcna_file)
  if ("name" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$name[seurat_obj$name %in% c("KIAA0586", "ly1ceb")] <- "OFD1"
  }

  pfm_core <- TFBSTools::getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
  )
  seurat_obj <- MotifScan(
    seurat_obj,
    species_genome = "hg38",
    pfm = pfm_core,
    EnsDb = EnsDb.Hsapiens.v86
  )
  motif_df <- GetMotifs(seurat_obj)
  write.csv(motif_df, file.path(out_dir, paste0(prefix, "_motif_df.csv")), row.names = FALSE)

  tf_genes <- unique(motif_df$gene_name)
  modules <- GetModules(seurat_obj)
  nongrey_genes <- subset(modules, module != "grey")$gene_name
  seurat_obj <- SetWGCNAGenes(seurat_obj, unique(c(tf_genes, nongrey_genes)))
  seurat_obj <- SetDatExpr(seurat_obj, group.by = NULL, group_name = NULL)

  model_params <- list(
    objective = "reg:squarederror",
    max_depth = 1,
    eta = 0.1,
    nthread = 16,
    alpha = 0.5
  )
  seurat_obj <- ConstructTFNetwork(seurat_obj, model_params)
  tf_network <- GetTFNetwork(seurat_obj)
  write.csv(tf_network, file.path(out_dir, paste0(prefix, "_TF_network.csv")), row.names = FALSE)

  seurat_obj <- AssignTFRegulons(seurat_obj, strategy = "C", reg_thresh = 0.1)
  seurat_obj <- RegulonScores(seurat_obj, target_type = "positive", ncores = 8)
  seurat_obj <- RegulonScores(seurat_obj, target_type = "negative", cor_thresh = -0.05, ncores = 8)

  tf_regulons <- GetTFRegulons(seurat_obj)
  hub_tf <- GetHubGenes(seurat_obj, n_hubs = 100) %>% subset(gene_name %in% tf_regulons$tf)
  write.csv(hub_tf, file.path(out_dir, paste0(prefix, "_hub_TFs.csv")), row.names = FALSE)

  Idents(seurat_obj) <- seurat_obj$cell_type
  marker_tfs <- FindAllMarkers(seurat_obj, features = unique(tf_regulons$tf))
  write.csv(marker_tfs, file.path(out_dir, paste0(prefix, "_marker_TFs.csv")), row.names = FALSE)

  group1 <- rownames(seurat_obj@meta.data)[seurat_obj$group == "JS"]
  group2 <- rownames(seurat_obj@meta.data)[seurat_obj$group == "ctrl"]
  dregs <- FindDifferentialRegulons(seurat_obj, barcodes1 = group1, barcodes2 = group2)
  write.csv(dregs, file.path(out_dir, paste0(prefix, "_differential_regulons.csv")), row.names = FALSE)

  pdf(file.path(out_dir, paste0(prefix, "_differential_regulons.pdf")), width = 7, height = 7)
  print(PlotDifferentialRegulons(seurat_obj, dregs))
  dev.off()

  seurat_obj <- RunEnrichrRegulons(seurat_obj, wait_time = 1)
  pdf(file.path(out_dir, paste0(prefix, "_module_regulatory_heatmap.pdf")), width = 12, height = 7)
  print(ModuleRegulatoryHeatmap(seurat_obj, feature = "delta", dendrogram = FALSE) +
          ggtitle("TFs only"))
  print(ModuleRegulatoryHeatmap(seurat_obj, feature = "delta", TFs_only = FALSE,
                                max_val = 5, dendrogram = FALSE) +
          ggtitle("All target genes"))
  dev.off()

  saveRDS(seurat_obj, file.path(out_dir, paste0(prefix, "_hdwgcna_TF_network.rds")))
}

cerebellum_hdwgcna <- file.path(paths$results, "hdwgcna/cerebellum_1021_all_celltypes_hdwgcna.rds")
if (file.exists(cerebellum_hdwgcna)) run_hdwgcna_tf_network(cerebellum_hdwgcna, "cerebellum")

kidney_hdwgcna <- file.path(paths$results, "hdwgcna/kidney_330_or_918_all_celltypes_hdwgcna.rds")
if (file.exists(kidney_hdwgcna)) run_hdwgcna_tf_network(kidney_hdwgcna, "kidney")
