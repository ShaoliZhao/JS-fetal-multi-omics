# Public Workflow Runbook

This runbook describes how to execute the public-safe workflows after controlled-access inputs are placed on a local workstation or analysis server.

## 1. Configure Project Root

Run commands from the `curated_code/` folder, or set:

```bash
export JS_MULTIOMICS_ROOT=/path/to/curated_code
```

The scripts write outputs to:

```text
results/
figures/
```

These folders are ignored by `.gitignore`.

## 2. Update Sample Sheet

Edit:

```text
metadata/snrna_sample_sheet_public_template.csv
```

The sample and QC columns are derived from `JS_supplementary_sample_information.xlsx`, sheet `Sample_summary`.

Replace the placeholder paths:

```text
data/controlled/cellranger/<sample_id>/outs
data/controlled/cellranger/<sample_id>/outs/filtered_feature_bc_matrix
```

with local Cell Ranger output paths.

## 3. Build snRNA-seq Atlases

```bash
Rscript analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R
```

Expected controlled inputs:

- Cell Ranger filtered matrices for cerebellum.
- Cell Ranger `outs/` folders for kidney SoupX correction.

Main outputs:

- `results/cerebellum_atlas_public_workflow.rds`
- `results/kidney_atlas_public_workflow.rds`
- QC UMAP and violin plots in `figures/`

## 4. Run Annotation, Markers, DE and GO

```bash
Rscript analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R \
  --object results/cerebellum_atlas_public_workflow.rds \
  --prefix cerebellum \
  --celltype-col celltype \
  --group-col group \
  --ident1 JS \
  --ident2 ctrl
```

For kidney:

```bash
Rscript analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R \
  --object results/kidney_atlas_public_workflow.rds \
  --prefix kidney \
  --celltype-col celltype \
  --group-col group \
  --ident1 JS \
  --ident2 ctrl
```

## 5. Run Spatial Label Transfer

```bash
Rscript analysis/04_spatial_transcriptomics/stereo_label_transfer_spatial_plots.R \
  --reference results/cerebellum_atlas_public_workflow.rds \
  --spatial results/stereo_cellbin_object.rds \
  --prefix ofd1_cellbin \
  --celltype-col celltype \
  --x-col spatial_x \
  --y-col spatial_y \
  --genes NRN1,BARHL1,S100B,CALB1
```

Expected controlled inputs:

- annotated cerebellar snRNA-seq Seurat object;
- Stereo-seq/cell-bin Seurat object with expression matrix and coordinate metadata.

## 6. Run Bulk Multi-omics Summaries

```bash
Rscript analysis/05_bulk_multiomics/bulk_multiomics_summary_plots.R \
  --sample-meta metadata/bulk_sample_metadata_public_template.csv \
  --matrix data/processed/bulk_rna_normalized_matrix.csv \
  --de-dir data/processed/bulk_de_tables \
  --prefix bulk_multiomics
```

This module works on public-safe processed summary tables, not raw FASTQ/BAM files.

## 7. Run Lineage and Velocity

Monocle2:

```bash
Rscript analysis/06_lineage_velocity/monocle_lineage_template.R \
  --object results/kidney_atlas_public_workflow.rds \
  --prefix kidney_npc_nephron \
  --subset-col celltype \
  --subset-values NPC,Podocyte,PT,LOH,LOH_DTL,DCT \
  --group-col group \
  --celltype-col celltype
```

scVelo:

```bash
python analysis/06_lineage_velocity/scvelo_velocity_template.py \
  --h5ad results/kidney_velocity_input.h5ad \
  --loom data/controlled/velocity/sample1.loom data/controlled/velocity/sample2.loom \
  --celltype-col celltype \
  --basis umap \
  --output-prefix results/kidney_velocity
```

## 8. Run Network and Signaling Modules

CellChat:

```bash
Rscript analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R \
  --mode cellchat \
  --object results/kidney_atlas_public_workflow.rds \
  --prefix kidney \
  --celltype-col celltype \
  --group-col group
```

pySCENIC input export:

```bash
Rscript analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R \
  --mode pyscenic_export \
  --object results/kidney_atlas_public_workflow.rds \
  --prefix kidney
```

hdWGCNA setup:

```bash
Rscript analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R \
  --mode hdwgcna \
  --object results/kidney_atlas_public_workflow.rds \
  --prefix kidney \
  --celltype-col celltype \
  --group-col group
```

## Reproducibility Boundary

Public BAM files support read-level provenance, but the final atlas workflows require matrix-level outputs and processed objects. This repository is therefore designed as a controlled-input workflow rather than a one-click public reanalysis package.
