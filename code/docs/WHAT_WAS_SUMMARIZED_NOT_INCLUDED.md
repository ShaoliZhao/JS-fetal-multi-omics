# What Was Summarized Instead of Included

The full workspace contains 152 detected code files. For a public GitHub portfolio, keeping all of them would make the project harder to understand and would expose unnecessary server-specific details.

## Not included: per-sample Cell Ranger commands

Reason:

- They are repetitive.
- They often encode private sample paths.
- The manuscript Methods already state the version, reference, intron handling, and output used.

Public replacement:

- `analysis/01_raw_processing/cellranger_count_template.sh`

## Not included: scripts dominated by absolute server paths

Reason:

- Paths such as `/data/StoredData/...` and `/mnt/...` are not portable.
- They can accidentally expose internal server structure or sample names.
- They make the GitHub look less professional.

Public replacement:

- `analysis/00_setup/project_config.R`
- modular scripts that use `paths$data_processed`, `paths$results`, and `paths$figures`

## Not included: duplicate exploratory analysis branches

Reason:

- Several scripts represent older attempts at the same CellChat, hdWGCNA, pySCENIC, Monocle, or plotting analysis.
- A portfolio should show the final analytical idea, not every scratchpad branch.

Public replacement:

- one representative template per method family
- clear Methods-to-code mapping

## Not included: large or sensitive data objects

Reason:

- Human fetal sequencing and clinical-genetic data should not be pushed to GitHub without institutional approval.
- `.RDS`, `.h5ad`, `.loom`, BAM, FASTQ, and raw matrices are too large and often sensitive.

Public replacement:

- data availability statement
- small processed example tables if approved
- interactive web assets only after de-identification and downsampling

## Recommended private/public split

- Keep `organized_code/` privately as provenance.
- Publish `curated_code/` as the clean portfolio.
- Publish the GitHub Pages website from `github_page/`.
- Add only approved figure images and processed summary tables.
