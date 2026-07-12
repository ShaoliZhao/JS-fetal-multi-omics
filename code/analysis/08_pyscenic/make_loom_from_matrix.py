"""Create a pySCENIC loom file from Matrix Market expression files.

Inputs are produced by spatial_cellbin_export_for_pyscenic.R:
matrix.mtx, genes.tsv and barcodes.tsv.
"""

from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from scipy.io import mmread


input_dir = Path("results/pyscenic/spatial_cellbin_matrix")
output_loom = Path("results/pyscenic/spatial_cellbin.loom")
output_loom.parent.mkdir(parents=True, exist_ok=True)

matrix = mmread(input_dir / "matrix.mtx").tocsr()
genes = pd.read_csv(input_dir / "genes.tsv", header=None)[0].astype(str).to_numpy()
barcodes = pd.read_csv(input_dir / "barcodes.tsv", header=None)[0].astype(str).to_numpy()

row_attrs = {"Gene": genes}
col_attrs = {"CellID": barcodes, "nGene": np.asarray((matrix > 0).sum(axis=0)).ravel()}
loompy.create(str(output_loom), matrix, row_attrs, col_attrs)
