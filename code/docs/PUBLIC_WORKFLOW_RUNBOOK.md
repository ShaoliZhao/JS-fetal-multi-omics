# Public Workflow Runbook

This repository is a controlled-input workflow. It is not meant to rebuild the
full atlas from BAM files alone. To rerun the analysis, provide local
Cell Ranger filtered matrices, SoupX-compatible kidney Cell Ranger output
folders, processed Seurat objects, featureCounts tables and pySCENIC/CellChat
intermediate files where noted.

Set the project root before running:

```bash
export JS_MULTIOMICS_ROOT=/path/to/curated_code
```

## 1. snRNA-seq atlas

Input:

- `metadata/snrna_sample_sheet_public_template.csv`
- local `filtered_feature_bc_matrix` folders
- SoupX-compatible Cell Ranger output folders for kidney samples

Run:

```bash
Rscript analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R
```

Key settings:

- cerebellum minimum features: 800
- kidney minimum features: 600
- kidney ambient RNA correction: SoupX
- doublet handling: DoubletFinder for cerebellum, scDblFinder for kidney
- Seurat integration/clustering: PCA/CCA, dims 1:30, resolution recorded in the script

## 2. Annotation, markers and cell-type DEG

Input:

- processed cerebellum or kidney Seurat object

Run:

```bash
Rscript analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R
```

This script writes marker tables, cell-type-specific JS/control contrasts, GO
tables and cell-type composition summaries.

## 3. Cell-bin spatial transcriptomics

Input:

- OFD1 and control Stereo-seq cell-bin Seurat objects with `Spatial` counts and
  `x/y` coordinates
- final cerebellar snRNA-seq reference object for label transfer

Run:

```bash
Rscript analysis/04_spatial_transcriptomics/cellbin_build_and_annotation.R
Rscript analysis/04_spatial_transcriptomics/st_label_transfer_and_celltype_deg.R
```

Key settings:

- `SCTransform(..., assay = "Spatial", return.only.var.genes = FALSE)`
- PCA/neighbor/UMAP dims 1:30
- merged cell-bin marker threshold: `logfc.threshold = 0.5`
- spatial cell-type DEG threshold: `p_val < 0.05` and `abs(avg_log2FC) > 0.5`
- PKC inner/outer comparison: `logfc.threshold = 0.3`, `min.diff.pct = 0.1`

## 4. Bulk RNA-seq

Input:

- featureCounts `.txt` files for kidney, cerebellum, frontal lobe and occipital lobe bulk RNA-seq
- `metadata/bulk_rnaseq_sample_sheet_template.csv` with final sample inclusion and JS/ctrl group labels
- proteomics differential tables and metabolomics heatmap tables under `data/processed/`

Run:

```bash
Rscript analysis/05_bulk_multiomics/bulk_rnaseq_deseq2_go.R
Rscript analysis/05_bulk_multiomics/proteomics_metabolomics_figures.R
```

Key settings:

- rows retained with `rowSums(counts) > 10`
- DESeq2 `fitType = "mean"`
- JS/control DEG threshold: adjusted P value < 0.05 and abs(log2FoldChange) > 0.5
- protein threshold: `Significant == "yes"` and abs(log2FC) > 0.5

## 5. Trajectory and RNA velocity

Run the relevant script after providing processed Seurat/h5ad/loom inputs:

```bash
Rscript analysis/06_lineage_velocity/monocle2_ddrtree_lineage.R
Rscript analysis/06_lineage_velocity/monocle3_graph_lineage.R
python analysis/06_lineage_velocity/scvelo_velocity_template.py
```

## 6. CellChat

Run:

```bash
Rscript analysis/07_cellchat/cellchat_public_workflows.R
```

Key settings:

- cerebellum CellChat: VZP/IN/PKC/GC/UBC/Cycling subset, full CellChatDB,
  `type = "truncatedMean"`, `trim = 0.01`, `min.cells = 10`
- kidney updated CellChat: NPC/Podocyte/PT/LOH/PEC/DCT/UB_CD/Stromal/Endo/Cycling,
  full CellChatDB, `trim = 0.01`, `min.cells = 10`
- kidney ECM branch: optional older ECM-Receptor workflow, `trim = 0.1`
- differential LR mapping: `thresh.pc = 0.1`, `thresh.fc = 0.05`, `thresh.p = 0.05`

## 7. pySCENIC

Run after pySCENIC loom and adjacency outputs are available:

```bash
Rscript analysis/08_pyscenic/export_seurat_counts_for_pyscenic.R
python analysis/08_pyscenic/make_loom_from_matrix.py
bash analysis/08_pyscenic/run_pyscenic_cli_template.sh
Rscript analysis/08_pyscenic/pyscenic_downstream_public.R
Rscript analysis/08_pyscenic/spatial_cellbin_export_for_pyscenic.R
```

Key outputs:

- regulon AUC added to Seurat metadata
- RSS plots by cell type/group/genotype
- top TF-target GO enrichment
- top TF-target network plots
- cerebellum RSS thresholds: celltype 2, group 1, name 1.2
- kidney RSS thresholds: celltype 2.5, group 1, geneotype 1.2

## 8. hdWGCNA

Run:

```bash
Rscript analysis/09_hdwgcna/hdwgcna_public_workflow.R
Rscript analysis/09_hdwgcna/hdwgcna_tf_network_public.R
```

Key settings:

- `gene_select = "fraction"`, `fraction = 0.05`
- Harmony by sample before metacell construction
- cerebellum 1021 branch uses Cellcycle/GCs/INs/MG/OPC&ODC/PKCs/UBCs/VZP
- cerebellum metacells: `min_cells = 50`, `k = 25`, `max_shared = 10`
- kidney branch uses nephron/stromal/endothelial major cell states
- kidney metacells: `min_cells = 90`, `k = 25`, `max_shared = 10`
- hdWGCNA TF branch uses JASPAR motif scan, `AssignTFRegulons(strategy = "C", reg_thresh = 0.1)`,
  positive regulons and negative regulons with `cor_thresh = -0.05`
