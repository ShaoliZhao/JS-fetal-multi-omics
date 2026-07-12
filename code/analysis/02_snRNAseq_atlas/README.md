# snRNA-seq Atlas

This folder contains the public-facing Seurat workflow for the cerebellum and kidney snRNA-seq atlases.

The script is intentionally more detailed than a minimal tutorial. It documents the real computational decisions behind the manuscript while replacing private server paths with a metadata-driven sample sheet.

## Files

- `seurat_qc_integration_template.R`: main public workflow for matrix loading, QC, doublet branches, CCA integration, clustering and atlas export.
- `final_object_summary.md`: inspected structure of the two final local Seurat objects used for figure-level analyses.
- `../../metadata/snrna_sample_sheet_public_template.csv`: de-identified sample sheet template.

## How to run

Run from the `curated_code/` folder after replacing the placeholder paths in the sample sheet:

```bash
Rscript analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R
```

For GitHub, keep the sample sheet as a template or point it to a tiny toy dataset. Do not commit real matrix folders, `.RDS` objects, BAM files or FASTQ files.

## What the sample sheet is

`snrna_sample_sheet_public_template.csv` is derived from `figure/suptable/JS_supplementary_sample_information.xlsx`, sheet `Sample_summary`.

It is not raw expression data. It is the table that lets the code loop over samples without hard-coded server paths.

Each row describes one sample:

- sample label and tissue;
- ctrl/JS group;
- genotype;
- gestational week;
- final nuclei count and QC summaries from the supplementary table;
- local path to a Cell Ranger `filtered_feature_bc_matrix` folder;
- local path to a Cell Ranger `outs/` folder when SoupX is needed;
- ambient RNA method;
- doublet method;
- minimum feature threshold.

For the public repository, the paths should remain placeholders or point to a tiny toy dataset. The real human fetal matrices and final Seurat objects should not be committed to GitHub.

## Cerebellum branch

The cerebellar atlas was built from Cell Ranger filtered feature-barcode matrices.

Public workflow highlights:

- loads samples by looping over the sample sheet;
- applies per-sample `nFeature_RNA` filtering;
- records mitochondrial percentage;
- documents the DoubletFinder branch with `pN = 0.25`, `pK = 0.09`, expected doublet rate 7.5 percent and PCs 1-10;
- integrates samples with Seurat v5 `IntegrateLayers(..., method = CCAIntegration)`;
- uses integrated CCA dimensions 1-30;
- clusters at resolution 0.4 for the main cerebellar atlas.

The final local cerebellar object inspected for this repository contains 194,911 nuclei and retains `pca`, `integrated.cca`, `umap` and `tsne` reductions.

## Kidney branch

The renal atlas includes an ambient RNA correction branch.

Public workflow highlights:

- loads Cell Ranger output folders with SoupX when `ambient_method = soupx`;
- uses `autoEstCont` or a documented manual contamination fraction when automatic estimates are unstable;
- applies scDblFinder singlet selection using the project branch rate `ncol(sce) * 8e-6`;
- integrates six final libraries with Seurat v5 CCA;
- uses integrated CCA dimensions 1-30;
- records clustering at resolutions 0.5 and 0.6 during iterative refinement;
- documents the broad and fine cell-type annotation fields used in `kd6rm`.

The final local kidney object inspected for this repository contains 66,906 nuclei and includes `celltype`, `ct`, `group` and genotype fields.

## Reproducibility boundary

The public Data Availability points to Cell Ranger BAM files. BAMs alone cannot exactly recreate the final Seurat objects, because this workflow starts from filtered matrices and processed count objects after QC, SoupX, singlet filtering and cluster review.

The code is therefore public-safe and methods-aligned, but full rerun requires controlled-access processed matrices or local Seurat objects.
