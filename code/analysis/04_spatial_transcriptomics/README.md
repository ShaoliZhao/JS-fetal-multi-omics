# Spatial Transcriptomics

This module represents the Stereo-seq cerebellar vermis workflow used for Figure 2, Figure 5, Supplementary Figure 5, and Supplementary Figure 6.

Core public-facing steps:

- normalize cell-bin expression objects
- transfer cell-type labels from the snRNA-seq cerebellar atlas
- plot predicted spatial cell states
- map marker genes and module scores
- summarize local spatial organization

For an interactive web atlas, this module would become the data-preparation backend: export downsampled coordinates, predicted labels, module scores, and selected gene expression values as compressed web tables.
