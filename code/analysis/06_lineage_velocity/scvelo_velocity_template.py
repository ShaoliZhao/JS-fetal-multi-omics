#!/usr/bin/env python3
"""scVelo RNA velocity workflow following the basic and dynamical tutorials."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv

project_root = Path(__file__).resolve().parents[2]

adata_file = project_root / "results" / "velocity" / "kidney_velocity_input.h5ad"
loom_files = [
    project_root / "data" / "external_controlled" / "velocity" / "sample1.loom",
    project_root / "data" / "external_controlled" / "velocity" / "sample2.loom",
]

analysis_name = "kidney_velocity"
celltype_col = "celltype"
group_col = "group"
basis = "umap"
velocity_mode = "stochastic"
run_dynamical_model = False

min_shared_counts = 20
n_top_genes = 3000
n_pcs = 30
n_neighbors = 30

out_dir = project_root / "results" / "06_lineage_velocity" / "scvelo" / analysis_name
fig_dir = project_root / "figures" / "06_lineage_velocity" / "scvelo" / analysis_name
out_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
scv.settings.figdir = str(fig_dir)

adata = sc.read_h5ad(adata_file)

loom_files = [x for x in loom_files if x.exists()]
if len(loom_files) > 0:
    loom_list = [scv.read(str(x), cache=True) for x in loom_files]
    if len(loom_list) == 1:
        ldata = loom_list[0]
    else:
        ldata = loom_list[0].concatenate(loom_list[1:], join="outer", batch_key="velocity_sample")
    adata = scv.utils.merge(adata, ldata)

scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=min_shared_counts,
    n_top_genes=n_top_genes,
)
scv.pp.moments(
    adata,
    n_pcs=n_pcs,
    n_neighbors=n_neighbors,
)

if run_dynamical_model:
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
else:
    scv.tl.velocity(adata, mode=velocity_mode)

scv.tl.velocity_graph(adata)
scv.tl.velocity_confidence(adata)

color = celltype_col if celltype_col in adata.obs.columns else None
scv.pl.velocity_embedding_stream(
    adata,
    basis=basis,
    color=color,
    save=f"_{analysis_name}_stream.png",
    show=False,
)
scv.pl.velocity_embedding(
    adata,
    basis=basis,
    color=color,
    arrow_length=3,
    arrow_size=2,
    dpi=150,
    save=f"_{analysis_name}_arrows.png",
    show=False,
)
scv.pl.velocity(
    adata,
    var_names=adata.var_names[:12],
    color=color,
    save=f"_{analysis_name}_selected_genes.png",
    show=False,
)

if run_dynamical_model:
    scv.tl.latent_time(adata)
    scv.pl.scatter(
        adata,
        color="latent_time",
        color_map="gnuplot",
        save=f"_{analysis_name}_latent_time.png",
        show=False,
    )

    top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(
        adata,
        var_names=top_genes,
        sortby="latent_time",
        col_color=celltype_col if celltype_col in adata.obs.columns else None,
        n_convolve=100,
        save=f"_{analysis_name}_latent_time_heatmap.png",
        show=False,
    )

    if celltype_col in adata.obs.columns:
        scv.tl.rank_dynamical_genes(adata, groupby=celltype_col)
        dynamic_genes = scv.get_df(adata, "rank_dynamical_genes/names")
        dynamic_genes.to_csv(out_dir / "rank_dynamical_genes_by_celltype.csv")

obs_cols = [x for x in [celltype_col, group_col, "velocity_length", "velocity_confidence"] if x in adata.obs.columns]
if len(obs_cols) > 0:
    adata.obs[obs_cols].to_csv(out_dir / "cell_velocity_metadata.csv")

adata.write_h5ad(out_dir / f"{analysis_name}_scvelo.h5ad", compression="gzip")
plt.close("all")
