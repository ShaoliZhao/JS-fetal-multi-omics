# Methods-to-code Map

This map turns the manuscript Methods into a concise public code structure.

## snRNA-seq processing

Methods:

- Cell Ranger v9/v10
- GRCh38-compatible reference
- intronic reads retained
- filtered feature-barcode matrices used for atlas construction
- kidney count matrices additionally corrected with SoupX before final atlas construction

Curated code:

- `analysis/01_raw_processing/cellranger_count_template.sh`
- `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`
- `metadata/snrna_sample_sheet_public_template.csv`
- `docs/PUBLIC_WORKFLOW_RUNBOOK.md`

## QC, doublet detection, and integration

Methods:

- Seurat v5
- per-sample QC
- DoubletFinder/scDblFinder branches
- CCA integration
- iterative removal of low-quality clusters
- cerebellum integration dimensions 1-30 and main resolution 0.4
- kidney integration dimensions 1-30 and iterative resolutions 0.5-0.8

Curated code:

- `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`
- `analysis/02_snRNAseq_atlas/final_object_summary.md`

The detailed historical scripts remain private because they contain sample-specific and server-specific paths. The curated code now records the Methods-level parameters and final-object metadata needed to evaluate the analysis logic.

## Cell-type annotation, markers, DE, and enrichment

Methods:

- COSG and Seurat marker identification
- FindMarkers for JS-control contrasts
- clusterProfiler/org.Hs.eg.db for GO biological process enrichment

Curated code:

- `analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R`

The annotation script now runs marker discovery, cell-type-specific JS-control DE, GO enrichment and composition summaries from a supplied Seurat object.

## Spatial transcriptomics

Methods:

- Stereo-seq FF V1.3
- cell-bin gene-expression matrices
- Seurat label transfer from snRNA-seq to spatial bins
- marker and module-score spatial maps

Curated code:

- `analysis/04_spatial_transcriptomics/stereo_label_transfer_spatial_plots.R`

The spatial script now performs Seurat label transfer, exports predicted cell-bin metadata, plots selected genes and plots cerebellar program scores.

## Bulk multi-omics

Methods:

- bulk RNA-seq, DIA proteomics, untargeted metabolomics
- PCA, differential feature analysis, GO/pathway enrichment
- cross-omic interpretation

Curated code:

- `analysis/05_bulk_multiomics/bulk_multiomics_summary_plots.R`

The bulk script now accepts processed feature matrices and DE tables for PCA, volcano plotting, GO enrichment and cross-table set-size summaries.

## Trajectory and velocity

Methods:

- Monocle 2/3 for GC, VZP, NPC, nephron, and UB/CD trajectories
- scVelo for RNA velocity from loom and h5ad objects

Curated code:

- `analysis/06_lineage_velocity/monocle_lineage_template.R`
- `analysis/06_lineage_velocity/scvelo_velocity_template.py`

These scripts now accept command-line inputs for Seurat objects, h5ad files and loom files.

## Networks and signaling

Methods:

- CellChat for ligand-receptor signaling
- hdWGCNA for co-expression modules
- pySCENIC for regulon inference

Curated code:

- `analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R`

This script now supports `cellchat`, `pyscenic_export` and `hdwgcna` modes.

## Figure visualization

Methods:

- ggplot2, patchwork, pheatmap, ComplexHeatmap, CellChat, hdWGCNA plotting outputs

Curated code:

- `analysis/08_visualization/publication_plot_helpers.R`
