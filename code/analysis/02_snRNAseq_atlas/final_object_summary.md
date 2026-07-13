# Final snRNA-seq Objects Used by the Manuscript

This page documents the processed Seurat objects used locally for figure-level analyses. The objects themselves are not included in the public GitHub repository because they contain human fetal single-nucleus expression data.

## Cerebellum object


Summary:

- Seurat object: 38,606 genes x 194,911 nuclei.
- Assay: `RNA`.
- Reductions: `pca`, `integrated.cca`, `umap`, `tsne`.
- Cluster field: `RNA_snn_res.0.4` / `seurat_clusters`.
- Group field: `group`.
- Final group counts: 93,677 control nuclei and 101,234 JS nuclei.
- Integration branch: Seurat v5 CCA integration using dimensions 1-30.
- Initial clustering resolution recorded in the object: 0.4.
- QC after filtering: median 2,213 genes per nucleus; mean 2,465 genes per nucleus.

Sample-level counts retained in the final local object:

```text
g2      9,341
g5h    13,179
g5v    13,963
g6     25,551
g7     23,390
lu     12,067
ly1ceb 16,697
lyh    14,672
lyv    11,357
pr     11,873
wxj    11,766
wxjba   9,134
wxs    11,455
zjj    10,466
```

Public repository note: these internal sample IDs should be replaced by de-identified labels in public files.

## Kidney object


Summary:

- Seurat object: 38,377 genes x 66,906 nuclei.
- Assay: `RNA`.
- Reductions: `pca`, `integrated.cca`, `umap`.
- Cluster fields: `RNA_snn_res.0.5`, `RNA_snn_res.0.6`, `seurat_clusters`.
- Main cell-type field: `celltype`.
- Fine annotation field: `ct`.
- Genotype field: `geneotype` in the local object. Public code should use the corrected spelling `genotype` for new outputs.
- Final group counts: 39,356 control nuclei and 27,550 JS nuclei.
- Final genotype counts: 39,356 control, 6,250 OFD1, 9,658 CEP290 and 11,642 TMEM67 nuclei.
- Integration branch: Seurat v5 CCA integration using dimensions 1-30.
- Final clustering resolution recorded in object: 0.5 and 0.6 metadata are both present after iterative reclustering.
- QC after filtering: median 2,107 genes per nucleus; mean 2,159 genes per nucleus.

Final broad kidney cell-type counts:

```text
NPC        2,842
Podocyte   7,059
PT         2,885
LOH        5,364
PEC    2,983
DCT        3,783
UB_CD     10,177
Stromal   21,009
SMC_JGC      786
Endo       5,998
MAC        1,083
Schwann      184
Cycling    2,753
```

Public repository note: the kidney workflow includes SoupX ambient RNA correction, scDblFinder singlet selection and iterative removal of low-quality or contaminating clusters before the final six-library object is annotated.
