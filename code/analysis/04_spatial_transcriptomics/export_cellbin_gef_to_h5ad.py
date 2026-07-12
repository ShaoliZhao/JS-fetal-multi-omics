"""Export Stereo-seq cell-bin GEF files to Seurat-flavored H5AD files."""

from pathlib import Path
import stereo as st

project_root = Path(__file__).resolve().parents[2]
input_dir = project_root / "data" / "external_controlled" / "spatial" / "gef"
output_dir = project_root / "data" / "external_controlled" / "spatial" / "h5ad"
output_dir.mkdir(parents=True, exist_ok=True)

samples = {
    "OFD1": input_dir / "OFD1.adjusted.cellbin.gef",
    "control_w18": input_dir / "control_w18.adjusted.cellbin.gef",
}

for sample, gef_file in samples.items():
    if not gef_file.exists():
        print(f"Skip {sample}: {gef_file} not found")
        continue

    st.io.read_gef_info(str(gef_file))
    data = st.io.read_gef(file_path=str(gef_file), bin_type="cell_bins")
    data.tl.cal_qc()
    data.tl.raw_checkpoint()
    st.io.stereo_to_anndata(
        data,
        flavor="seurat",
        output=str(output_dir / f"{sample}_cellbin.h5ad"),
    )
