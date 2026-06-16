"""
10x Visium human breast cancer — download, load, inspect, and QC spatial plots.

Dataset: V1_Breast_Cancer_Block_A_Section_1 (10x Genomics Visium demo).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

SAMPLE_ID = "V1_Breast_Cancer_Block_A_Section_1"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "figures"
DATA_ROOT = Path(__file__).resolve().parents[1] / "data"
DATA_DIR = DATA_ROOT / SAMPLE_ID


def download_visium_data(sample_id: str = SAMPLE_ID, data_root: Path = DATA_ROOT) -> Path:
    """Download 10x Visium breast cancer data and return the on-disk folder."""
    data_root.mkdir(parents=True, exist_ok=True)
    print(f"Downloading Visium dataset: {sample_id}")
    sq.datasets.visium(sample_id=sample_id, base_dir=data_root)
    print(f"Saved to: {data_root / sample_id}")
    return data_root / sample_id


def load_visium_with_scanpy(data_dir: Path, sample_id: str = SAMPLE_ID) -> sc.AnnData:
    """Load Visium count matrix and spatial metadata with Scanpy."""
    adata = sc.read_visium(path=data_dir)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_id
    return adata


def inspect_anndata(adata: sc.AnnData) -> None:
    """Print a concise summary of the AnnData object."""
    print("\n=== AnnData summary ===")
    print(adata)
    print("\n--- obs (spots) ---")
    print(adata.obs.head())
    print("\n--- var (genes) ---")
    print(adata.var.head())
    print("\n--- obsm keys ---")
    print(list(adata.obsm.keys()))
    print("\n--- uns keys ---")
    print(list(adata.uns.keys()))
    if "spatial" in adata.uns:
        library_id = list(adata.uns["spatial"].keys())[0]
        print(f"\n--- spatial library: {library_id} ---")
        print(list(adata.uns["spatial"][library_id].keys()))


def plot_qc_spatial(adata: sc.AnnData, output_dir: Path = OUTPUT_DIR) -> None:
    """Plot tissue image, total counts, and genes per spot."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(output_dir)
    sc.set_figure_params(dpi=120, facecolor="white")

    sc.pp.calculate_qc_metrics(adata, inplace=True)

    library_id = list(adata.uns["spatial"].keys())[0]

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color=None,
        ax=ax,
        show=False,
        title="H&E tissue image",
    )
    fig.savefig(output_dir / "01_tissue_image.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color="total_counts",
        ax=ax,
        show=False,
        title="Total counts per spot",
    )
    fig.savefig(output_dir / "02_total_counts_per_spot.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(
        adata,
        library_id=library_id,
        color="n_genes_by_counts",
        ax=ax,
        show=False,
        title="Number of genes per spot",
    )
    fig.savefig(output_dir / "03_genes_per_spot.png", bbox_inches="tight")
    plt.close(fig)

    print(f"\nFigures saved to: {output_dir}")


def main() -> None:
    data_dir = download_visium_data()
    adata = load_visium_with_scanpy(data_dir)
    inspect_anndata(adata)
    plot_qc_spatial(adata)


if __name__ == "__main__":
    main()
