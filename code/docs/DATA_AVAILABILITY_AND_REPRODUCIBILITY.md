# Data Availability and Reproducibility Notes

The public code is organized to show the computational logic of the manuscript while protecting human fetal sample-level data and private server paths.

## What the public BAM files can support

Cell Ranger BAM files document aligned reads and can support read-level inspection, locus-level checks and some downstream analyses that explicitly consume BAM/loom-style inputs.

However, the Seurat atlas workflow in this repository starts from:

- Cell Ranger `filtered_feature_bc_matrix` folders for the cerebellar atlas.
- Cell Ranger output folders plus SoupX-corrected count matrices for the kidney atlas.
- Singlet-filtered and iteratively reviewed Seurat objects for figure-level analyses.

Therefore, BAM files alone are not sufficient to reproduce every final Seurat object exactly. The missing pieces include Cell Ranger matrix outputs, SoupX correction decisions, singlet calls, cluster-review exclusions and manual annotation decisions.

## Recommended wording for GitHub

Use language like this in the repository README:

```text
This repository provides public-safe analysis templates and figure-generation workflows corresponding to the manuscript. Human fetal expression matrices and final Seurat objects are not distributed in this repository. The snRNA-seq atlas scripts are runnable with controlled-access filtered feature-barcode matrices or processed count objects matching the sample sheet template.
```

## Why not publish every historical script

The original analysis folder contains many one-off commands, absolute server paths, repeated Cell Ranger jobs and exploratory branches. Those are useful provenance records, but they make a public repository harder to evaluate. The curated code keeps the computational decisions that matter:

- matrix-level input conventions;
- per-sample QC thresholds;
- kidney ambient RNA correction by SoupX;
- doublet removal branches;
- CCA integration dimensions;
- clustering resolutions;
- iterative cluster review;
- final annotation fields used in the figures.

## Minimal data needed to rerun the public workflow

To run `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`, provide:

- one row per sample in `metadata/snrna_sample_sheet_public_template.csv`, derived from the supplementary sample information table and extended with local input paths;
- local paths to `filtered_feature_bc_matrix` folders for cerebellum;
- local paths to Cell Ranger `outs/` folders for kidney SoupX correction;
- R packages listed in `analysis/00_setup/packages.R` plus optional branch-specific packages such as SoupX, scDblFinder and DoubletFinder.

The repository can also include a tiny synthetic 10x matrix for CI/testing, but not the real fetal expression matrices.

Downstream scripts are also controlled-input workflows. They are runnable when supplied with local processed Seurat objects, h5ad/loom files or public-safe summary tables, but they should not be expected to run from BAM files alone.
