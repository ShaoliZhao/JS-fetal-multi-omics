# Fetal Joubert Syndrome Multi-omics Atlas

This repository hosts a public-facing project website and controlled-input analysis workflows for a fetal Joubert syndrome multi-omics study.

The website focuses on:

- main manuscript figures only;
- the computational workflow behind the atlas;
- public-safe code organization;
- transparent reproducibility boundaries.

## Website

Open `index.html` locally, or enable GitHub Pages from the repository root.

## Code

The `code/` folder contains curated workflows:

- snRNA-seq atlas construction from sample-sheet-driven matrix inputs;
- annotation, marker, cell-type DE and GO enrichment;
- Stereo-seq/cell-bin label transfer and spatial plotting;
- bulk multi-omics summaries;
- Monocle2 and scVelo lineage dynamics;
- CellChat, pySCENIC export and hdWGCNA setup.

Start with:

```bash
code/docs/PUBLIC_WORKFLOW_RUNBOOK.md
```

## Data Boundary

Human fetal matrices and final Seurat objects are not included. Public BAM files support read-level provenance, but exact final atlas reconstruction requires controlled-access matrix-level outputs and processed objects.
