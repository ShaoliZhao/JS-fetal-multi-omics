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

Key settings retained from the original scripts:

- `SCTransform(..., assay = "Spatial", return.only.var.genes = FALSE)`
- PCA/neighbor/UMAP dims 1:30
- merged cell-bin marker threshold: `logfc.threshold = 0.5`
- spatial cell-type DEG threshold: `p_val < 0.05` and `abs(avg_log2FC) > 0.5`
- PKC inner/outer comparison: `logfc.threshold = 0.3`, `min.diff.pct = 0.1`

## 4. Bulk RNA-seq

Input:

- featureCounts `.txt` files for cerebellum and kidney bulk RNA-seq

Run:

```bash
Rscript analysis/05_bulk_multiomics/bulk_deseq2_go_figures.R
```

Key settings:

- rows retained with `rowSums(counts) > 10`
- DESeq2 `fitType = "mean"`
- cerebellum DEG: `padj < 0.05` and `abs(log2FoldChange) > 0.5`
- kidney genotype contrasts: `pvalue < 0.05` and `abs(log2FoldChange) > 0.5`

## 5. Trajectory and RNA velocity

Run the relevant script after providing processed Seurat/h5ad/loom inputs:

```bash
Rscript analysis/06_lineage_velocity/monocle_lineage_template.R
python analysis/06_lineage_velocity/scvelo_velocity_template.py
```

## 6. CellChat

Run:

```bash
Rscript analysis/07_cellchat/cellchat_public_workflows.R
```

Key settings:

- kidney CellChat: `type = "truncatedMean"`, `trim = 0.01`, `min.cells = 10`
- spatial CellChat: cell-bin coordinates, `ratio = 0.5`, `tol = 10`,
  `interaction.range = 20`, `contact.range = 10`
- spatial OFD1 trim: 0.1
- spatial control trim: 0.5

## 7. pySCENIC

Run after pySCENIC loom and adjacency outputs are available:

```bash
Rscript analysis/08_pyscenic/pyscenic_downstream_public.R
Rscript analysis/08_pyscenic/spatial_cellbin_export_for_pyscenic.R
python analysis/08_pyscenic/make_loom_from_matrix.py
```

Key outputs:

- regulon AUC added to Seurat metadata
- RSS plots by cell type/group/genotype
- top TF-target GO enrichment
- top TF-target network plots

## 8. hdWGCNA

Run:

```bash
Rscript analysis/09_hdwgcna/hdwgcna_public_workflow.R
```

Key settings:

- `gene_select = "fraction"`, `fraction = 0.05`
- Harmony by sample before metacell construction
- cerebellum modules are run separately for GCs, VZP and PKCs
- kidney module analysis uses nephron/stromal/endothelial major cell states
