# Export annotated cell-bin Seurat data for the browser-based spatial atlas.
# The web page stores one coordinate table plus one expression vector per gene.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(jsonlite)
})

cellbin_rds <- file.path(paths$results, "spatial_label_transfer_deg/cellbin_label_transfer.rds")
out_dir <- file.path(paths$results, "spatial_point_atlas_export")
dir.create(file.path(out_dir, "genes"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(cellbin_rds)
DefaultAssay(obj) <- "Spatial"

meta <- obj@meta.data
meta$cell_id <- rownames(meta)

if ("Spatial_1" %in% colnames(meta)) {
  meta$x_raw <- meta$Spatial_1
  meta$y_raw <- meta$Spatial_2
} else if (all(c("x", "y") %in% colnames(meta))) {
  meta$x_raw <- meta$x
  meta$y_raw <- meta$y
} else {
  stop("Spatial coordinates were not found in metadata.")
}

meta$sample <- as.character(meta$orig.ident)
meta$sample[meta$sample %in% c("control", "Control", "ctrl", "w18", "W18")] <- "ctrl"
meta$sample[meta$sample != "ctrl"] <- "JS"

if ("predicted_class" %in% colnames(meta)) {
  meta$celltype <- as.character(meta$predicted_class)
} else if ("celltype" %in% colnames(meta)) {
  meta$celltype <- as.character(meta$celltype)
} else {
  meta$celltype <- as.character(meta$seurat_clusters)
}

meta$x <- ave(meta$x_raw, meta$sample, FUN = function(v) {
  if (max(v) == min(v)) return(rep(0.5, length(v)))
  (v - min(v)) / (max(v) - min(v))
})
meta$y <- ave(meta$y_raw, meta$sample, FUN = function(v) {
  if (max(v) == min(v)) return(rep(0.5, length(v)))
  (v - min(v)) / (max(v) - min(v))
})

meta <- meta[colnames(obj), ]
points <- meta[, c("cell_id", "sample", "x", "y", "celltype")]
colnames(points)[1] <- "id"
write_json(list(points = points), file.path(out_dir, "points.json"), dataframe = "rows", auto_unbox = TRUE)

counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
counts <- counts[, points$id]

genes_to_export <- intersect(
  c("ATOH1", "PCP4", "WNT7B", "SOX2", "MKI67", "FOXJ1", "OFD1", "CEP290", "TMEM67"),
  rownames(counts)
)

# To export every gene for a production atlas, replace the line above with:
# genes_to_export <- rownames(counts)

manifest <- list(
  demo = FALSE,
  samples = list(
    list(id = "ctrl", label = "Control"),
    list(id = "JS", label = "JS / OFD1")
  ),
  genes = sort(genes_to_export)
)
write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

for (gene in genes_to_export) {
  values <- as.numeric(log1p(counts[gene, ]))
  names(values) <- NULL
  file_gene <- gsub("[^A-Za-z0-9_.-]+", "_", gene)
  write_json(
    list(gene = gene, values = values),
    file.path(out_dir, "genes", paste0(file_gene, ".json")),
    auto_unbox = TRUE
  )
}

message("Spatial atlas data written to: ", out_dir)
message("Copy manifest.json, points.json and the genes folder into github_release_xulabb/spatial_atlas/data/")
