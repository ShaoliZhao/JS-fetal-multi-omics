#!/usr/bin/env bash
set -euo pipefail

# Generic Cell Ranger template.
# The original project used one command per sample; those files are summarized here.

CELLRANGER="${CELLRANGER:-cellranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/path/to/refdata-gex-GRCh38}"
FASTQS="${FASTQS:-/path/to/fastqs}"
SAMPLE_ID="${SAMPLE_ID:-sample_id}"
OUT_ID="${OUT_ID:-run_count_sample_id}"

"${CELLRANGER}" count \
  --id="${OUT_ID}" \
  --transcriptome="${TRANSCRIPTOME}" \
  --fastqs="${FASTQS}" \
  --sample="${SAMPLE_ID}" \
  --create-bam=true \
  --localcores="${LOCAL_CORES:-16}" \
  --localmem="${LOCAL_MEM_GB:-128}"
