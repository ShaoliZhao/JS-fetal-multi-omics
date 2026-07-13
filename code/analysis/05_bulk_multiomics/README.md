# Bulk RNA-seq, Proteomics and Metabolomics

This folder is split into two workflows.

- `bulk_rnaseq_deseq2_go.R`: assembles featureCounts tables, runs DESeq2 JS/control contrasts, writes `diffbox` tables, DEG sets, GO results, PCA/sample-distance plots, Venn plots and cilia-gene heatmaps.
- `proteomics_metabolomics_figures.R`: analyzes protein differential tables, protein GO, RNA-protein overlap, proteomics Venn plots and metabolomics heatmaps.

The bulk RNA-seq script uses `metadata/bulk_rnaseq_sample_sheet_template.csv`.
Replace the relative `count_file` values with local controlled-access featureCounts outputs before running.
