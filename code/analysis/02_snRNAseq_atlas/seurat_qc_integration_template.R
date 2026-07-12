# Public snRNA-seq atlas workflow for the fetal Joubert syndrome study.
#
# This script is designed for GitHub presentation and controlled re-analysis.
# It mirrors the core Methods without publishing private sample paths or full
# human fetal Seurat objects. It is runnable after the user fills in
# metadata/snrna_sample_sheet_public_template.csv with local 10x output paths.
#
# Important reproducibility note:
# The public Data Availability for this manuscript points to Cell Ranger BAMs.
# BAM files document aligned reads, but this workflow starts from Cell Ranger
# filtered feature-barcode matrices or SoupX-corrected counts. Reviewers cannot
# rebuild the final Seurat objects from BAM files alone unless the matching
# Cell Ranger output matrices or controlled-access processed objects are also
# available.

source("analysis/00_setup/project_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

optional_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this branch. Install it or set the sample sheet to skip this branch.", call. = FALSE)
  }
}

sample_sheet_file <- file.path(paths$metadata, "snrna_sample_sheet_public_template.csv")

expected_columns <- c(
  "tissue", "sample_id", "public_label", "group", "genotype",
  "gestational_weeks", "n_cells_final", "median_nCount_RNA",
  "median_nFeature_RNA", "mean_nCount_RNA", "mean_nFeature_RNA",
  "include_in_public_workflow", "cellranger_out_dir", "matrix_dir",
  "ambient_method", "manual_contamination_fraction", "doublet_method",
  "min_features", "notes"
)

read_public_sample_sheet <- function(path = sample_sheet_file) {
  samples <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(expected_columns, colnames(samples))
  if (length(missing) > 0) {
    stop("Sample sheet is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  samples <- samples %>%
    mutate(
      include_in_public_workflow = tolower(include_in_public_workflow) %in% c("true", "t", "yes", "y", "1"),
      tissue_key = tolower(tissue),
      manual_contamination_fraction = suppressWarnings(as.numeric(manual_contamination_fraction)),
      min_features = suppressWarnings(as.integer(min_features)),
      gestational_weeks = suppressWarnings(as.numeric(gestational_weeks))
    )
  samples
}

resolve_project_path <- function(path) {
  if (is.na(path) || path == "") {
    return(path)
  }
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(project_root, path)
}

add_qc_metrics <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

read_filtered_matrix <- function(matrix_dir) {
  matrix_dir <- resolve_project_path(matrix_dir)
  if (!dir.exists(matrix_dir)) {
    stop("Matrix directory does not exist: ", matrix_dir, call. = FALSE)
  }
  Read10X(data.dir = matrix_dir)
}

load_cellranger_soupx_counts <- function(cellranger_out_dir, manual_contamination_fraction = NA_real_) {
  optional_package("SoupX")
  cellranger_out_dir <- resolve_project_path(cellranger_out_dir)
  if (!dir.exists(cellranger_out_dir)) {
    stop("Cell Ranger output directory does not exist: ", cellranger_out_dir, call. = FALSE)
  }

  soup_channel <- SoupX::load10X(cellranger_out_dir)
  if (is.na(manual_contamination_fraction)) {
    soup_channel <- SoupX::autoEstCont(soup_channel)
  } else {
    soup_channel <- SoupX::setContaminationFraction(soup_channel, manual_contamination_fraction)
  }
  SoupX::adjustCounts(soup_channel)
}

make_object_from_sample <- function(sample_row) {
  min_features <- ifelse(is.na(sample_row$min_features), 600L, sample_row$min_features)

  counts <- if (tolower(sample_row$ambient_method) == "soupx") {
    load_cellranger_soupx_counts(
      cellranger_out_dir = sample_row$cellranger_out_dir,
      manual_contamination_fraction = sample_row$manual_contamination_fraction
    )
  } else {
    read_filtered_matrix(sample_row$matrix_dir)
  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_row$sample_id,
    min.features = min_features
  )
  obj <- add_qc_metrics(obj)
  obj$sample_id <- sample_row$sample_id
  obj$public_label <- sample_row$public_label
  obj$tissue <- sample_row$tissue
  obj$group <- sample_row$group
  obj$genotype <- sample_row$genotype
  obj$gestational_weeks <- sample_row$gestational_weeks
  obj
}

basic_qc_filter <- function(obj, min_features = 600L, max_features = Inf, max_percent_mt = Inf) {
  subset(
    obj,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      percent.mt <= max_percent_mt
  )
}

preprocess_for_review <- function(obj, dims = 1:30, resolution = 0.6) {
  obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = dims, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE) %>%
    RunUMAP(dims = dims, verbose = FALSE)
}

run_scdblfinder <- function(obj, expected_rate_per_cell = 8e-6) {
  optional_package("scDblFinder")
  optional_package("SingleCellExperiment")
  optional_package("scater")

  sce <- Seurat::as.SingleCellExperiment(obj)
  expected_rate <- ncol(sce) * expected_rate_per_cell
  sce <- scDblFinder::scDblFinder(sce, dbr = expected_rate)
  singlet_cells <- colnames(sce)[sce$scDblFinder.class == "singlet"]
  subset(obj, cells = singlet_cells)
}

run_doubletfinder_branch <- function(obj, pcs = 1:10, expected_doublet_rate = 0.075, pN = 0.25, pK = 0.09) {
  optional_package("DoubletFinder")

  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)

  homotypic_prop <- DoubletFinder::modelHomotypic(obj$seurat_clusters)
  n_exp <- round(expected_doublet_rate * ncol(obj))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))

  obj <- DoubletFinder::doubletFinder(
    obj,
    PCs = pcs,
    pN = pN,
    pK = pK,
    nExp = n_exp_adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )

  df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(df_col) == 0) {
    stop("DoubletFinder did not create a DF.classifications column.", call. = FALSE)
  }
  final_df_col <- df_col[length(df_col)]
  singlet_cells <- colnames(obj)[obj@meta.data[[final_df_col]] == "Singlet"]

  list(
    object = obj,
    singlets = subset(obj, cells = singlet_cells),
    expected_doublets = n_exp,
    homotypic_adjusted_doublets = n_exp_adj,
    classification_column = final_df_col
  )
}

integrate_seurat_objects <- function(objects, dims = 1:30, resolution = 0.5, add_cell_ids = names(objects)) {
  stopifnot(length(objects) >= 2)

  atlas <- merge(objects[[1]], y = objects[-1], add.cell.ids = add_cell_ids)
  atlas <- atlas %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = max(dims), verbose = FALSE)

  atlas <- IntegrateLayers(
    object = atlas,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )

  atlas <- FindNeighbors(atlas, reduction = "integrated.cca", dims = dims, verbose = FALSE)
  atlas <- FindClusters(atlas, resolution = resolution, verbose = FALSE)
  atlas <- RunUMAP(atlas, reduction = "integrated.cca", dims = dims, verbose = FALSE)
  JoinLayers(atlas)
}

plot_atlas_qc <- function(atlas, prefix) {
  p_cluster <- DimPlot(atlas, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    theme_js_atlas() +
    ggtitle(paste(prefix, "clusters"))
  p_group <- DimPlot(atlas, group.by = "group") +
    theme_js_atlas() +
    ggtitle(paste(prefix, "group"))
  p_sample <- DimPlot(atlas, group.by = "sample_id") +
    theme_js_atlas() +
    ggtitle(paste(prefix, "sample"))
  p_feature <- VlnPlot(atlas, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) +
    theme_js_atlas() +
    NoLegend()

  save_panel(p_cluster, paste0(prefix, "_umap_clusters.pdf"), width = 6, height = 5)
  save_panel(p_group, paste0(prefix, "_umap_group.pdf"), width = 5, height = 4)
  save_panel(p_sample, paste0(prefix, "_umap_sample.pdf"), width = 6, height = 5)
  save_panel(p_feature, paste0(prefix, "_nFeature_by_sample.pdf"), width = 7, height = 3)
}

annotate_kidney_public_celltypes <- function(kidney_obj) {
  # Mapping taken from the final kd6rm object used for the manuscript figures.
  kidney_obj$ct <- as.character(kidney_obj$seurat_clusters)
  fine_map <- c(
    "0" = "Str_Cor_Out", "1" = "LOH", "2" = "UB", "3" = "Endo_Cap",
    "4" = "Str_Cor_Inn", "5" = "Podo_early", "6" = "DCT", "7" = "Str_var",
    "8" = "CD_PC", "9" = "Str_Med", "10" = "LOH_DTL", "11" = "PT",
    "12" = "NPC", "13" = "Cycling", "14" = "Podo_mature", "15" = "UB_Tip",
    "16" = "Fibro", "17" = "MAC", "18" = "Endo_Art", "19" = "SMC_JGC",
    "20" = "Endo_Lymph", "21" = "Schwann", "22" = "CD_IC"
  )
  kidney_obj$ct <- unname(fine_map[kidney_obj$ct])

  kidney_obj$celltype <- as.character(kidney_obj$seurat_clusters)
  kidney_obj$celltype[kidney_obj$celltype %in% c("0", "4", "7", "9", "16")] <- "Stromal"
  kidney_obj$celltype[kidney_obj$celltype %in% c("5", "14")] <- "Podocyte"
  kidney_obj$celltype[kidney_obj$celltype == "1"] <- "LOH"
  kidney_obj$celltype[kidney_obj$celltype == "10"] <- "LOH_DTL"
  kidney_obj$celltype[kidney_obj$celltype == "6"] <- "DCT"
  kidney_obj$celltype[kidney_obj$celltype == "11"] <- "PT"
  kidney_obj$celltype[kidney_obj$celltype == "12"] <- "NPC"
  kidney_obj$celltype[kidney_obj$celltype %in% c("2", "8", "15", "22")] <- "UB_CD"
  kidney_obj$celltype[kidney_obj$celltype %in% c("3", "18", "20")] <- "Endo"
  kidney_obj$celltype[kidney_obj$celltype == "13"] <- "Cycling"
  kidney_obj$celltype[kidney_obj$celltype == "17"] <- "MAC"
  kidney_obj$celltype[kidney_obj$celltype == "19"] <- "SMC_JGC"
  kidney_obj$celltype[kidney_obj$celltype == "21"] <- "Schwann"

  kidney_obj$celltype <- factor(
    kidney_obj$celltype,
    levels = c("NPC", "Podocyte", "PT", "LOH", "LOH_DTL", "DCT", "UB_CD",
               "Stromal", "SMC_JGC", "Endo", "MAC", "Schwann", "Cycling")
  )
  kidney_obj
}

write_object_summary <- function(obj, prefix) {
  summary_file <- file.path(paths$results, paste0(prefix, "_object_summary.txt"))
  capture.output({
    cat("Object:", prefix, "\n")
    cat("Genes x cells:", paste(dim(obj), collapse = " x "), "\n")
    cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
    cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n\n")
    cat("Samples:\n")
    print(table(obj$sample_id, useNA = "ifany"))
    cat("\nGroups:\n")
    print(table(obj$group, useNA = "ifany"))
    if ("celltype" %in% colnames(obj@meta.data)) {
      cat("\nCell types:\n")
      print(table(obj$celltype, useNA = "ifany"))
    }
    cat("\nQC summary:\n")
    print(summary(obj@meta.data[, intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt"), colnames(obj@meta.data)), drop = FALSE]))
  }, file = summary_file)
  summary_file
}

build_cerebellum_atlas <- function(samples) {
  ceb_samples <- samples %>%
    filter(include_in_public_workflow, tissue_key == "cerebellum")

  ceb_objects <- purrr::pmap(
    ceb_samples,
    function(...) {
      sample_row <- as.data.frame(list(...), stringsAsFactors = FALSE)
      obj <- make_object_from_sample(sample_row)
      obj <- basic_qc_filter(obj, min_features = ifelse(is.na(sample_row$min_features), 800L, sample_row$min_features))

      if (tolower(sample_row$doublet_method) == "doubletfinder") {
        df_result <- run_doubletfinder_branch(obj, pcs = 1:10, expected_doublet_rate = 0.075, pN = 0.25, pK = 0.09)
        obj <- df_result$singlets
        obj <- CreateSeuratObject(GetAssayData(obj, assay = "RNA", layer = "counts"), project = sample_row$sample_id)
        obj <- add_qc_metrics(obj)
        obj$sample_id <- sample_row$sample_id
        obj$public_label <- sample_row$public_label
        obj$tissue <- sample_row$tissue
        obj$group <- sample_row$group
        obj$genotype <- sample_row$genotype
        obj$gestational_weeks <- sample_row$gestational_weeks
      }
      obj
    }
  )
  names(ceb_objects) <- ceb_samples$sample_id

  ceb_atlas <- integrate_seurat_objects(
    objects = ceb_objects,
    dims = 1:30,
    resolution = 0.4,
    add_cell_ids = ceb_samples$sample_id
  )

  plot_atlas_qc(ceb_atlas, "cerebellum")
  write_object_summary(ceb_atlas, "cerebellum")
  saveRDS(ceb_atlas, file.path(paths$results, "cerebellum_atlas_public_workflow.rds"))
  ceb_atlas
}

build_kidney_atlas <- function(samples) {
  kidney_samples <- samples %>%
    filter(include_in_public_workflow, tissue_key == "kidney")

  kidney_objects <- purrr::pmap(
    kidney_samples,
    function(...) {
      sample_row <- as.data.frame(list(...), stringsAsFactors = FALSE)
      obj <- make_object_from_sample(sample_row)

      # The renal branch used SoupX correction, low-quality-cluster review,
      # and scDblFinder singlet selection before final CCA integration.
      obj <- preprocess_for_review(obj, dims = 1:30, resolution = 0.6)
      if (tolower(sample_row$doublet_method) == "scdblfinder") {
        obj <- run_scdblfinder(obj, expected_rate_per_cell = 8e-6)
      }

      obj <- CreateSeuratObject(round(GetAssayData(obj, assay = "RNA", layer = "counts")), project = sample_row$sample_id)
      obj <- add_qc_metrics(obj)
      obj$sample_id <- sample_row$sample_id
      obj$public_label <- sample_row$public_label
      obj$tissue <- sample_row$tissue
      obj$group <- sample_row$group
      obj$genotype <- sample_row$genotype
      obj$gestational_weeks <- sample_row$gestational_weeks
      obj
    }
  )
  names(kidney_objects) <- kidney_samples$sample_id

  kidney_atlas <- integrate_seurat_objects(
    objects = kidney_objects,
    dims = 1:30,
    resolution = 0.6,
    add_cell_ids = kidney_samples$sample_id
  )

  # In the private analysis, clusters with low complexity, erythroid/urothelial
  # contamination or likely doublet profiles were reviewed and excluded, then
  # the remaining six-library object was reclustered. The final public mapping
  # below documents the annotation used by kd6rm.
  kidney_atlas <- annotate_kidney_public_celltypes(kidney_atlas)

  plot_atlas_qc(kidney_atlas, "kidney")
  write_object_summary(kidney_atlas, "kidney")
  saveRDS(kidney_atlas, file.path(paths$results, "kidney_atlas_public_workflow.rds"))
  kidney_atlas
}

main <- function() {
  samples <- read_public_sample_sheet()

  message("Building cerebellum atlas from filtered feature-barcode matrices.")
  cerebellum_atlas <- build_cerebellum_atlas(samples)

  message("Building kidney atlas with SoupX/scDblFinder branch where requested.")
  kidney_atlas <- build_kidney_atlas(samples)

  invisible(list(cerebellum = cerebellum_atlas, kidney = kidney_atlas))
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
