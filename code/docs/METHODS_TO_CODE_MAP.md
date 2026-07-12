# Methods-to-Code Map

This map links the manuscript computational methods to the curated public code.

## Raw processing and sample metadata

- Sample metadata:
  `metadata/snrna_sample_sheet_public_template.csv`
- Raw processing notes:
  `analysis/01_raw_processing/`

Public BAM files provide read-level provenance. The analysis workflows begin
from matrix-level Cell Ranger outputs or processed controlled-access objects.

## snRNA-seq atlas construction

- Main script:
  `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`

Methods represented:

- looping over Cell Ranger filtered matrices
- kidney SoupX ambient RNA removal
- doublet removal
- Seurat QC, integration, PCA/UMAP and clustering
- final object summary

## Cell-type annotation, DEG and GO

- Main script:
  `analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R`

Methods represented:

- marker gene discovery
- cell-type-specific JS/control contrasts
- GO Biological Process enrichment
- cell-type composition summaries

## Spatial transcriptomics

- Spatial helper functions:
  `analysis/04_spatial_transcriptomics/STimport_public.R`
- Cell-bin object construction and annotation:
  `analysis/04_spatial_transcriptomics/cellbin_build_and_annotation.R`
- Label transfer and spatial cell-type DEG:
  `analysis/04_spatial_transcriptomics/st_label_transfer_and_celltype_deg.R`

Methods represented:

- Stereo-seq cell-bin count matrix loading
- ENSEMBL-to-symbol conversion
- spatial coordinate embedding
- SCT/PCA/UMAP/clustering
- manual cerebellar spatial domain annotation
- snRNA-seq-to-cell-bin label transfer
- cell-type-specific spatial DEG and GO

## Bulk RNA-seq

- Main script:
  `analysis/05_bulk_multiomics/bulk_deseq2_go_figures.R`

Methods represented:

- featureCounts matrix assembly
- DESeq2 JS/control contrasts
- GO enrichment and DEG set export

## Lineage and velocity

- Monocle2:
  `analysis/06_lineage_velocity/monocle2_ddrtree_lineage.R`
- Monocle3:
  `analysis/06_lineage_velocity/monocle3_graph_lineage.R`
- scVelo:
  `analysis/06_lineage_velocity/scvelo_velocity_template.py`

## CellChat

- Main script:
  `analysis/07_cellchat/cellchat_public_workflows.R`

Methods represented:

- cerebellum single-nucleus CellChat
- kidney single-nucleus CellChat
- optional kidney ECM-Receptor CellChat branch
- merged JS/control signaling comparison
- differential ligand-receptor mapping

## pySCENIC

- Downstream regulon analysis:
  `analysis/08_pyscenic/pyscenic_downstream_public.R`
- Seurat counts export:
  `analysis/08_pyscenic/export_seurat_counts_for_pyscenic.R`
- pySCENIC shell template:
  `analysis/08_pyscenic/run_pyscenic_cli_template.sh`
- Spatial cell-bin export:
  `analysis/08_pyscenic/spatial_cellbin_export_for_pyscenic.R`
- Loom creation:
  `analysis/08_pyscenic/make_loom_from_matrix.py`

Methods represented:

- regulon AUC extraction from loom files
- RSS analysis by cell type/group/genotype
- TF-target GO enrichment
- TF-target network visualization

## hdWGCNA

- Main script:
  `analysis/09_hdwgcna/hdwgcna_public_workflow.R`
- TF-network script:
  `analysis/09_hdwgcna/hdwgcna_tf_network_public.R`

Methods represented:

- metacell construction
- signed co-expression network inference
- module eigengenes and kME
- hub gene export
- cerebellum 1021 all-major-celltype module analysis
- kidney all-major-celltype module analysis
- JASPAR motif scan and hdWGCNA differential regulons
