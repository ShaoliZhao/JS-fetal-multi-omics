#!/usr/bin/env bash
set -euo pipefail

# pySCENIC command template curated from:
# - code/Rpr0/cebJS/pyscenic/0928/try0928.sh
# - code/Rpr0/kid0326/pyscenic/try330pyscenic.kid.sh
#
# Required inputs:
#   1. loom file made from expression counts
#   2. human TF list
#   3. hg38 motif ranking feather file
#   4. motif annotation table

LOOM_FILE="${1:-results/pyscenic/cerebellum.loom}"
TF_LIST="${2:-resources/pyscenic/hs_hgnc_tfs.txt}"
RANKINGS="${3:-resources/pyscenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"
MOTIF_ANNOTATIONS="${4:-resources/pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl}"
PREFIX="${5:-sample}"

pyscenic grn \
  --num_workers 16 \
  --output "${PREFIX}_adj.sample.tsv" \
  --method grnboost2 \
  "${LOOM_FILE}" \
  "${TF_LIST}"

pyscenic ctx \
  "${PREFIX}_adj.sample.tsv" \
  "${RANKINGS}" \
  --annotations_fname "${MOTIF_ANNOTATIONS}" \
  --expression_mtx_fname "${LOOM_FILE}" \
  --mode "dask_multiprocessing" \
  --output "${PREFIX}_reg.csv" \
  --num_workers 8 \
  --mask_dropouts

pyscenic aucell \
  "${LOOM_FILE}" \
  "${PREFIX}_reg.csv" \
  --output "${PREFIX}_SCENIC.loom" \
  --num_workers 8
