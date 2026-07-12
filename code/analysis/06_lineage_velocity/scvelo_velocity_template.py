#!/usr/bin/env python3
"""scVelo RNA velocity workflow for public-safe re-analysis.

Example:
python analysis/06_lineage_velocity/scvelo_velocity_template.py \
  --h5ad results/kidney_velocity_input.h5ad \
  --loom data/controlled/velocity/sample1.loom data/controlled/velocity/sample2.loom \
  --celltype-col celltype \
  --basis umap \
  --output-prefix results/kidney_velocity
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata
import scanpy as sc
import scvelo as scv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scVelo from an h5ad object and one or more loom files.")
    parser.add_argument("--h5ad", required=True, help="Annotated h5ad converted from the final Seurat object.")
    parser.add_argument("--loom", nargs="+", required=True, help="Velocyto loom files generated from Cell Ranger BAMs.")
    parser.add_argument("--celltype-col", default="celltype", help="obs column used for coloring velocity plots.")
    parser.add_argument("--group-col", default="group", help="obs column used for group summaries.")
    parser.add_argument("--basis", default="umap", help="Embedding basis, usually umap.")
    parser.add_argument("--mode", default="stochastic", choices=["stochastic", "dynamical", "deterministic"])
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--output-prefix", default="results/velocity_output")
    return parser.parse_args()


def require_file(path: str) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"Input file not found: {resolved}")
    return resolved


def read_and_merge_velocity(h5ad_path: Path, loom_paths: list[Path]) -> anndata.AnnData:
    adata = sc.read_h5ad(h5ad_path)
    loom_objects = [scv.read(str(path), cache=True) for path in loom_paths]
    if len(loom_objects) == 1:
        ldata = loom_objects[0]
    else:
        ldata = loom_objects[0].concatenate(loom_objects[1:], join="outer", batch_key="velocity_sample")
    return scv.utils.merge(adata, ldata)


def main() -> None:
    args = parse_args()
    h5ad_path = require_file(args.h5ad)
    loom_paths = [require_file(path) for path in args.loom]
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    adata = read_and_merge_velocity(h5ad_path, loom_paths)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=args.n_top_genes)
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)

    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode=args.mode)
    scv.tl.velocity_graph(adata)

    color = args.celltype_col if args.celltype_col in adata.obs.columns else None
    scv.pl.velocity_embedding_stream(
        adata,
        basis=args.basis,
        color=color,
        save=f"_{output_prefix.name}_stream.png",
        show=False,
    )
    scv.pl.velocity_embedding(
        adata,
        basis=args.basis,
        color=color,
        arrow_length=3,
        arrow_size=2,
        dpi=150,
        save=f"_{output_prefix.name}_arrows.png",
        show=False,
    )

    adata.write_h5ad(f"{output_prefix}.h5ad")
    if args.group_col in adata.obs.columns:
        adata.obs[[args.group_col, args.celltype_col]].to_csv(f"{output_prefix}_cell_metadata.csv")


if __name__ == "__main__":
    main()
