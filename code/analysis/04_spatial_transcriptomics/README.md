# Spatial Transcriptomics

This folder contains the public cell-bin Stereo-seq workflow.

- `STimport_public.R`: colors and spatial plotting helpers.
- `export_cellbin_gef_to_h5ad.py`: Stereo-seq `.gef` cell-bin export to
  Seurat-style `.h5ad`.
- `cellbin_build_and_annotation.R`: OFD1/control cell-bin object construction,
  SCT/PCA/UMAP/clustering, manual cell-type annotation and PKC comparison.
- `st_label_transfer_and_celltype_deg.R`: snRNA-seq reference label transfer
  and cell-type-specific OFD1 vs control spatial DEG/GO.
- `export_spatial_point_atlas_data.R`: exports cell-bin coordinates and
  per-gene expression vectors for the browser-based spatial atlas.

The spatial CellChat workflow is separated into `analysis/07_cellchat/`.
