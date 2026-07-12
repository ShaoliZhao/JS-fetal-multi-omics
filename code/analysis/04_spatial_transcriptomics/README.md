# Spatial Transcriptomics

This folder contains the public cell-bin Stereo-seq workflow curated from
`code/ST_code_backup`.

`ST_code_backup` is not a final public analysis folder. It is a historical
backup of spatial scripts, including cell-bin object construction, spatial
label transfer, spatial CellChat, pySCENIC input conversion and older trials.
The useful manuscript-facing pieces were kept here as shorter, path-safe
scripts:

- `STimport_public.R`: colors and spatial plotting helpers.
- `cellbin_build_and_annotation.R`: OFD1/control cell-bin object construction,
  SCT/PCA/UMAP/clustering, manual cell-type annotation and PKC comparison.
- `st_label_transfer_and_celltype_deg.R`: snRNA-seq reference label transfer
  and cell-type-specific OFD1 vs control spatial DEG/GO.

The spatial CellChat workflow is separated into `analysis/07_cellchat/`.
