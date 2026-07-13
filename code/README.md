# Curated Code for the Fetal Joubert Syndrome Multi-omics Atlas

This folder is the public-facing, controlled-input code portfolio for the fetal Joubert syndrome multi-omics project.

It intentionally does **not** include every historical script. The original workspace contains many exploratory files, sample-specific Cell Ranger commands, absolute server paths, and intermediate analysis branches. Those are useful as a lab notebook, but they are not ideal for GitHub, manuscript review, or internship presentation.

## What is kept

- Core analysis workflows that reflect the manuscript Methods, supplementary sample table and final object structure.
- Runnable R/Python/Shell workflows for the major computational steps once controlled-access matrices or processed objects are supplied.
- Figure-oriented code organization.
- Public-safe relative input paths instead of server paths.
- Notes explaining which code was summarized rather than retained.

## What is not kept

- Per-sample Cell Ranger command files.
- Hard-coded `/data/StoredData/...` or `/mnt/...` server paths.
- Raw human fetal data paths and real filtered matrices.
- Old duplicated trial scripts.
- Large intermediate objects such as `.RDS`, `.h5ad`, `.loom`, BAM, FASTQ, or raw matrices.

## Folder layout

```text
curated_code/
├── README.md
├── analysis/
│   ├── 00_setup/
│   ├── 01_raw_processing/
│   ├── 02_snRNAseq_atlas/
│   ├── 03_annotation_de_enrichment/
│   ├── 04_spatial_transcriptomics/
│   ├── 05_bulk_multiomics/
│   ├── 06_lineage_velocity/
│   ├── 07_cellchat/
│   ├── 08_pyscenic/
│   └── 09_hdwgcna/
├── metadata/
│   └── snrna_sample_sheet_public_template.csv
└── docs/
    ├── METHODS_TO_CODE_MAP.md
    ├── CURATION_SUMMARY.md
    ├── DATA_AVAILABILITY_AND_REPRODUCIBILITY.md
    ├── PUBLIC_WORKFLOW_RUNBOOK.md
    ├── WHAT_WAS_SUMMARIZED_NOT_INCLUDED.md
    └── INTERACTIVE_SPATIAL_WEB_RECOMMENDATION.md
```

## Figure-level map

- Figure 1: cohort design and omics overview.
- Figure 2: snRNA-seq atlas, spatial mapping, and cross-omic overview.
- Figure 3: bulk RNA-seq, proteomics, metabolomics, enrichment, and module summaries.
- Figure 4: cerebellar granule lineage, proliferation, trajectory, regulon, and spatial validation.
- Figure 5: VZP/glial scaffold programs, hdWGCNA, CellChat, and pseudo-space analysis.
- Figure 6: kidney stromal remodeling, hdWGCNA, regulons, and collagen signaling.
- Figure 7: nephron progenitor-derived output failure, gene-set scoring, TF networks, and IHC-linked interpretation.
- Figure 8: UB/CD branch trajectory, fibrosis/injury programs, WNT and SPP1 signaling.

## How to use this as a portfolio

Use this folder for GitHub presentation. Keep the full `organized_code/` or original `code/` as private provenance. In the public repository, pair these curated scripts with a data availability statement and small non-identifying example inputs.

The snRNA-seq workflow is not advertised as runnable from public BAM files alone. It is runnable from local filtered feature-barcode matrices, SoupX-compatible Cell Ranger output folders, or controlled-access processed Seurat objects that match the sample sheet.

## Run Model

The public repository should contain code, metadata templates and approved summary tables. It should not contain real human fetal matrices or final Seurat objects. To run the workflows locally or on a server, set:

```bash
export JS_MULTIOMICS_ROOT=/path/to/curated_code
```

Then replace the relative paths in `metadata/snrna_sample_sheet_public_template.csv` with the local controlled-access Cell Ranger output folders.

## Main Entry Points

- snRNA-seq atlas construction: `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`
- annotation, marker, DE and GO analysis: `analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R`
- Stereo-seq/cell-bin construction, annotation, label transfer and spatial DEG:
  `analysis/04_spatial_transcriptomics/`
- bulk RNA-seq DESeq2 and GO summaries:
  `analysis/05_bulk_multiomics/bulk_rnaseq_deseq2_go.R`
- proteomics, metabolomics and RNA-protein overlap figures:
  `analysis/05_bulk_multiomics/proteomics_metabolomics_figures.R`
- Monocle2 trajectories and scVelo velocity: `analysis/06_lineage_velocity/`
- CellChat workflows: `analysis/07_cellchat/cellchat_public_workflows.R`
- pySCENIC/SCENIC export, command template and RSS downstream analysis:
  `analysis/08_pyscenic/`
- hdWGCNA co-expression modules and hdWGCNA TF networks:
  `analysis/09_hdwgcna/`

This is the version I would show to reviewers, collaborators, or pharmaceutical internship interviewers.
