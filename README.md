# Fetal Joubert Syndrome Multi-omics Atlas

This repository hosts a public-facing project website and controlled-input analysis workflows for a fetal Joubert syndrome multi-omics study.

The website focuses on:

- main manuscript figures only;
- a browser-based spatial atlas demo for cell-bin gene lookup;
- the computational workflow behind the atlas;
- public-safe code organization;
- transparent reproducibility boundaries.

## Website

Open `index.html` locally, or enable GitHub Pages from the repository root.

The interactive spatial atlas entry point is:

```text
spatial_atlas/index.html
```

## Code

The `code/` folder contains curated workflows:

- snRNA-seq atlas construction from sample-sheet-driven matrix inputs;
- annotation, marker, cell-type DE and GO enrichment;
- Stereo-seq/cell-bin object construction, label transfer, spatial DEG and atlas data export;
- bulk RNA-seq DESeq2/GO workflows separated from proteomics and metabolomics figure code;
- Monocle2 and scVelo lineage dynamics;
- CellChat, pySCENIC export and hdWGCNA setup.

Start with:

```bash
code/docs/PUBLIC_WORKFLOW_RUNBOOK.md
```

## Data Boundary

Human fetal matrices and final Seurat objects are not included. Public BAM files support read-level provenance, but exact final atlas reconstruction requires controlled-access matrix-level outputs and processed objects.
