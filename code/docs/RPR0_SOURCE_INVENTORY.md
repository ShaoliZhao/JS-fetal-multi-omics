# Rpr0 Source Inventory

This file records how the public workflows were curated from `code/Rpr0`.

## CellChat

Kept as public workflow:

- `code/Rpr0/cebJS/cellchat/0423/cc0423.R`
  - cerebellum cluster-to-major-celltype mapping
  - VZP/IN/PKC/GC/UBC/Cycling subset
  - JS vs control CellChat
  - `computeCommunProb(type = "truncatedMean", trim = 0.01)`
  - `filterCommunication(min.cells = 10)`
  - merged comparison, bubble plots, heatmaps and DE-mapped LR pairs
- `code/Rpr0/kid0326/cellchat/cckid0330.R`
  - updated kidney CellChat
  - NPC/Podocyte/PT/LOH/LOH_DTL/DCT/UB_CD/Stromal/Endo/Cycling subset
  - JS vs control, full CellChatDB
  - SPP1-focused UB/CD bubble plot
- `code/Rpr0/kidney/cellchat/1208final/cc1208.R`
  - older ECM-Receptor focused branch
  - retained as optional ECM workflow, not the main kidney workflow

Summarized, not copied:

- `code/Rpr0/cebJS/cellchat/GC/GCcc.R`
- `code/Rpr0/cebJS/cellchat/anacellchatJSceb1013.R`
- older kidney CellChat scripts

Reason: these contain exploratory object names, repeated code and mixed tissue
branches. Their useful settings are already represented in the public workflow.

## pySCENIC / SCENIC

Kept as public workflow:

- `code/Rpr0/cebJS/pyscenic/0928/change0928.py`
  - counts-to-loom conversion for cerebellum
- `code/Rpr0/cebJS/pyscenic/0928/try0928.sh`
  - `pyscenic grn`, `ctx`, `aucell`
- `code/Rpr0/cebJS/pyscenic/ana928/scenicceb.R`
  - loom reading with SCopeLoomR
  - regulon AUC added to Seurat metadata
  - RSS by celltype/group/name
  - TF-target GO and TF network plots
- `code/Rpr0/kid0326/pyscenic/kidpyscenic.R`
  - Matrix Market export from kidney Seurat counts
- `code/Rpr0/kid0326/pyscenic/change.py`
  - Matrix Market-to-loom conversion
- `code/Rpr0/kid0326/pyscenic/try330pyscenic.kid.sh`
  - updated kidney pySCENIC commands
- `code/Rpr0/kid0326/pyscenic/rana/pyscenic.kid.downstream.R`
  - kidney SCENIC/RSS downstream analysis

Summarized, not copied:

- dated kidney pySCENIC folders `0909`, `0915`, `0917`, `0918`
- old pySCENIC folders under `code/Rpr0/kidney/old`

Reason: these are earlier attempts. The public scripts use the updated
`kid0326` branch while preserving the RSS thresholds from the older final
figure scripts where they were explicit.

## hdWGCNA

Kept as public workflow:

- `code/Rpr0/cebJS/hdWGCNA/1021/finalcebhdwgcna.R`
  - cerebellum final all-celltype hdWGCNA branch
  - ribosomal/blacklist gene removal noted in comments
  - `gene_select = "fraction"`, `fraction = 0.05`
  - Harmony by sample
  - metacells by Sample and cell_type
  - `min_cells = 50`, `k = 25`, `max_shared = 10`
  - signed network, module eigengenes, kME, hub genes, DMEs
- `code/Rpr0/kid0326/hdwgcna/hdwgcna.kidney330.R`
  - updated kidney hdWGCNA branch
  - kidney major cell states
  - all-module network and downstream enrichment/DMEs
- `code/Rpr0/kidney/hdwgcna/wgcna0918.R`
  - earlier kidney workflow with `min_cells = 90`
  - retained as parameter provenance for the kidney public workflow
- `code/Rpr0/cebJS/hdWGCNA/1021/TF/TFceb1021.R`
  - hdWGCNA TF network, JASPAR motif scan, differential regulons
- `code/Rpr0/kidney/hdwgcna/TF/TFhdWGCNAkidney.R`
  - kidney hdWGCNA TF network and differential regulons

Summarized, not copied:

- older `TF1010`, `cebjs_hdwgcna1011.R`, and `old/hdWGCNA` scripts

Reason: these repeat the same package tutorial pattern with older input objects.
