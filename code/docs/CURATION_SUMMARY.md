# Curation Summary

## Starting point

The private workspace contained 152 detected code files:

- 79 R scripts
- 13 Python scripts
- 51 shell scripts
- 9 JavaScript/Illustrator helper scripts

Those files are valuable as provenance, but they are not ideal as a public GitHub repository because many are duplicated, exploratory, sample-specific, or tied to private server paths.

## Public-facing result

The curated version keeps representative workflow scripts, a metadata template and documentation files:

- `analysis/00_setup/project_config.R`
- `analysis/00_setup/packages.R`
- `analysis/01_raw_processing/cellranger_count_template.sh`
- `analysis/02_snRNAseq_atlas/seurat_qc_integration_template.R`
- `analysis/02_snRNAseq_atlas/final_object_summary.md`
- `analysis/03_annotation_de_enrichment/annotation_markers_de_enrichment.R`
- `analysis/04_spatial_transcriptomics/stereo_label_transfer_spatial_plots.R`
- `analysis/05_bulk_multiomics/bulk_multiomics_summary_plots.R`
- `analysis/06_lineage_velocity/monocle_lineage_template.R`
- `analysis/06_lineage_velocity/scvelo_velocity_template.py`
- `analysis/07_networks_signaling/hdwgcna_pyscenic_cellchat_templates.R`
- `analysis/08_visualization/publication_plot_helpers.R`
- `metadata/snrna_sample_sheet_public_template.csv`
- `docs/PUBLIC_WORKFLOW_RUNBOOK.md`

## What this curation is meant to show

This code presents the project as a coherent computational workflow:

1. Process nuclear RNA-seq data.
2. Build integrated cerebellum and kidney atlases from matrix-level inputs.
3. Annotate cell types and identify disease-associated cell states.
4. Integrate spatial transcriptomics with snRNA-seq labels.
5. Summarize bulk transcriptomic, proteomic, and metabolomic signals.
6. Reconstruct lineage dynamics and RNA velocity.
7. Infer regulons, co-expression modules, and ligand-receptor signaling.
8. Generate publication-style visualizations.

## Key refinement made after inspecting the final local objects

The snRNA-seq template was expanded from a short generic script into a workflow that records the real analysis decisions:

- cerebellum: 14-sample matrix loading, `nFeature_RNA > 800`, CCA dimensions 1-30, clustering resolution 0.4 and DoubletFinder settings;
- kidney: SoupX ambient RNA correction, scDblFinder singlet filtering, six-library final atlas construction, CCA dimensions 1-30 and clustering resolutions 0.5-0.6;
- final object summaries for `ceb14s` and `kd6rm` are documented without publishing the objects.

The downstream modules were also expanded from short examples into controlled-input scripts with command-line arguments, file checks and defined outputs:

- annotation, marker, cell-type DE, GO and composition summaries;
- Stereo-seq/cell-bin label transfer, spatial gene plots and program-score plots;
- bulk matrix PCA, volcano plots, GO and cross-table DE set summaries;
- Monocle2 lineage reconstruction and scVelo velocity;
- CellChat comparison, pySCENIC input export and hdWGCNA setup.

## What remains private

- Full historical code provenance: `organized_code/`
- Original lab scripts: `code/` and code-like files under `figure/`
- Raw human fetal sequencing data
- Full Seurat, h5ad, loom, BAM, and FASTQ objects
- Sample-specific command files
- Files containing private server paths or sensitive sample identifiers
- Full processed matrices needed for exact atlas re-analysis

## Recommended GitHub strategy

Use `curated_code/` for the public repository, and keep `organized_code/` private. The curated version is cleaner, easier to evaluate, and better aligned with the manuscript story.
