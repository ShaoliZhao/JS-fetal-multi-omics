# Curation Summary

The original project code contains exploratory scripts, server paths,
sample-specific command fragments and repeated figure trials. This public
folder keeps the core manuscript-facing analysis in a cleaner form.

## Kept

- snRNA-seq loading, QC, doublet handling, integration and annotation logic
- kidney SoupX branch
- cell-bin Stereo-seq object construction and label transfer
- bulk RNA-seq DESeq2/GO workflow
- CellChat, pySCENIC/SCENIC and hdWGCNA as separate modules
- figure-oriented outputs and public-safe metadata templates

## Simplified

- repeated per-sample Cell Ranger commands are summarized in
  `analysis/01_raw_processing/`
- old trial scripts are represented by the final workflow steps
- absolute server paths are replaced with project-relative input variables
- spatial helper functions are reduced to the colors and plotting utilities
  actually needed by the cell-bin workflow

## Removed from the public code

- historical backup folders
- unused visualization helper folder
- duplicated GO/plot trials that were not tied to manuscript outputs
- large intermediate objects and private human fetal data matrices

## Important boundary

The code is runnable after controlled-access inputs are supplied locally, but
it is not intended to reconstruct the full atlas from BAM files alone.
