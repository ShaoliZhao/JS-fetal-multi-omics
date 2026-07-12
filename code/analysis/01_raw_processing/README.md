# Raw Processing

The original workspace contains many sample-specific Cell Ranger commands. For a public GitHub repository, those commands are better represented by one template plus the Methods text.

Manuscript Methods summary:

- snRNA-seq reads were processed with Cell Ranger v9 or v10.
- A GRCh38-compatible human reference was used.
- Intronic reads were retained for nuclear libraries.
- Filtered feature-barcode matrices were used for the cerebellar atlas.
- Kidney libraries were further evaluated with SoupX/scDblFinder branches before final atlas construction.
- Cell Ranger BAM files were retained for Data Availability and read-level provenance, but the Seurat workflows start from matrix-level outputs.

Public repository choice:

- Keep `cellranger_count_template.sh`.
- Do not publish per-sample FASTQ paths, BAM paths, sample IDs, or server folders.
- Do not imply that the final Seurat objects can be exactly recreated from BAM files alone; exact rerun requires matching Cell Ranger matrices or controlled-access processed objects.
