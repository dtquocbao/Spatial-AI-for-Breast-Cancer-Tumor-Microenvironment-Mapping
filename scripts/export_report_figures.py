"""Export publication figures for reports/conference_report.qmd."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay, roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "reports" / "figures"
RAW_H5 = ROOT / "data" / "raw" / "V1_Breast_Cancer_Block_A_Section_1"
CLUSTERED_H5 = ROOT / "data" / "processed" / "breast_cancer_visium_clustered.h5ad"
SPATIAL_H5 = ROOT / "data" / "processed" / "breast_cancer_spatial_analysis.h5ad"

CLUSTER_ANNOTATION = {
    "0": "Unknown",
    "1": "Immune-rich",
    "2": "Tumor-like epithelial",
    "3": "Tumor-like epithelial",
    "4": "Mixed epithelial",
    "5": "Stromal-like",
    "6": "Inflammatory immune",
    "7": "Unknown",
    "8": "Unknown",
    "9": "Adipose-rich",
}

IMMUNE_GENES = ["PTPRC", "CD3D", "CD74", "HLA-DRA", "HLA-DPB1"]
MARKER_GENES = [
    "EPCAM", "KRT8", "KRT18", "PTPRC", "CD3D", "MS4A1",
    "COL1A1", "COL1A2", "PECAM1", "VWF",
]


def _save(fig: plt.Figure, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}")


def export_exploration_figures() -> sc.AnnData:
    adata = sc.read_visium(RAW_H5)
    adata.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    library_id = list(adata.uns["spatial"].keys())[0]

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(adata, library_id=library_id, color=None, ax=ax, show=False, title="H&E tissue")
    _save(fig, "fig01_tissue_image.png")

    for color, fname, title in [
        ("total_counts", "fig02_total_counts_per_spot.png", "Total counts per spot"),
        ("n_genes_by_counts", "fig03_genes_per_spot.png", "Genes per spot"),
    ]:
        fig, ax = plt.subplots(figsize=(6, 6))
        sc.pl.spatial(adata, library_id=library_id, color=color, ax=ax, show=False, title=title)
        _save(fig, fname)
    return adata


def export_clustering_figures(adata_clustered: sc.AnnData) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata_clustered, color="leiden_0_5", ax=ax, show=False, title="Leiden clusters (UMAP)")
    _save(fig, "fig04_umap_leiden.png")

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(
        adata_clustered, color="leiden_0_5", img_key="hires", size=1.5,
        ax=ax, show=False, title="Leiden clusters on tissue",
    )
    _save(fig, "fig05_spatial_leiden.png")


def export_marker_figure(adata_clustered: sc.AnnData) -> None:
    adata = adata_clustered.copy()
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    markers = [g for g in MARKER_GENES if g in adata.var_names]
    sc.pl.dotplot(
        adata, var_names=markers, groupby="leiden_0_5",
        standard_scale="var", show=False,
    )
    _save(plt.gcf(), "fig06_marker_dotplot.png")


def export_spatial_figures(adata_spatial: sc.AnnData) -> sc.AnnData:
    if "spatial_connectivities" not in adata_spatial.obsp:
        sq.gr.spatial_neighbors(adata_spatial)
    if "leiden_0_5_nhood_enrichment" not in adata_spatial.uns:
        sq.gr.nhood_enrichment(adata_spatial, cluster_key="leiden_0_5")

    sq.pl.nhood_enrichment(adata_spatial, cluster_key="leiden_0_5", show=False)
    _save(plt.gcf(), "fig07_nhood_enrichment.png")

    adata_spatial.obs["bio_region"] = (
        adata_spatial.obs["leiden_0_5"].astype(str).map(CLUSTER_ANNOTATION).astype("category")
    )
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.spatial(
        adata_spatial, color="bio_region", img_key="hires", size=1.5,
        ax=ax, show=False, title="Biologically annotated regions",
    )
    _save(fig, "fig08_bio_regions.png")
    return adata_spatial


def export_ml_figures(adata_spatial: sc.AnnData) -> None:
    adata = adata_spatial.copy()
    immune_genes = [g for g in IMMUNE_GENES if g in adata.var_names]
    sc.tl.score_genes(adata, immune_genes, score_name="immune_score")
    adata.obs["immune_binary"] = (adata.obs["immune_score"] > 1).astype(int)

    adata_ml = adata[:, adata.var["highly_variable"]].copy() if "highly_variable" in adata.var else adata
    X = adata_ml.X.toarray() if not isinstance(adata_ml.X, np.ndarray) else adata_ml.X
    y = adata.obs["immune_binary"].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=42, stratify=y,
    )
    rf = RandomForestClassifier(
        n_estimators=300, min_samples_leaf=5, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )
    rf.fit(X_train, y_train)
    prob = rf.predict_proba(X_test)[:, 1]

    fig, ax = plt.subplots(figsize=(6, 5))
    RocCurveDisplay.from_predictions(y_test, prob, ax=ax, name="Random Forest")
    ax.set_title("Immune-high spot prediction (ROC)")
    _save(fig, "fig09_immune_hotspot_roc.png")

    importances = pd.DataFrame({
        "gene": adata_ml.var_names,
        "importance": rf.feature_importances_,
    }).sort_values("importance", ascending=False)
    top = importances.head(20)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["gene"][::-1], top["importance"][::-1])
    ax.set_xlabel("Feature importance")
    ax.set_title("Top genes predicting immune-high spots")
    _save(fig, "fig10_feature_importance.png")

    # Spatial context comparison (notebook 07)
    epithelial = [g for g in ["EPCAM", "KRT8", "KRT18"] if g in adata.var_names]
    stromal = [g for g in ["COL1A1", "COL1A2", "MGP"] if g in adata.var_names]
    sc.tl.score_genes(adata, epithelial, score_name="epithelial_score")
    sc.tl.score_genes(adata, stromal, score_name="stromal_score")

    W = adata.obsp["spatial_connectivities"]
    def neighbor_avg(values):
        counts = np.asarray(W.sum(axis=1)).flatten()
        return (W @ values) / np.maximum(counts, 1)

    context = np.column_stack([
        neighbor_avg(adata.obs["immune_score"].values),
        neighbor_avg(adata.obs["epithelial_score"].values),
        neighbor_avg(adata.obs["stromal_score"].values),
        adata.obsm["spatial"],
    ])
    X_context = np.hstack([X, context])

    indices = np.arange(adata.n_obs)
    train_idx, test_idx, y_train, y_test = train_test_split(
        indices, y, test_size=0.25, random_state=42, stratify=y,
    )
    rf_gene = RandomForestClassifier(
        n_estimators=300, min_samples_leaf=5, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )
    rf_ctx = RandomForestClassifier(
        n_estimators=300, min_samples_leaf=5, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )
    rf_gene.fit(X[train_idx], y_train)
    rf_ctx.fit(X_context[train_idx], y_train)
    prob_gene = rf_gene.predict_proba(X[test_idx])[:, 1]
    prob_ctx = rf_ctx.predict_proba(X_context[test_idx])[:, 1]
    auc_gene = roc_auc_score(y_test, prob_gene)
    auc_ctx = roc_auc_score(y_test, prob_ctx)

    comparison = pd.DataFrame({
        "Model": ["Gene expression only", "Gene + spatial context"],
        "ROC AUC": [auc_gene, auc_ctx],
    })
    comparison.to_csv(FIG_DIR / "model_comparison_auc.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(comparison["Model"], comparison["ROC AUC"])
    ax.set_ylabel("ROC AUC")
    ax.set_ylim(0, 1)
    ax.set_title("Spatial context vs. gene expression")
    plt.xticks(rotation=15, ha="right")
    _save(fig, "fig11_model_comparison_auc.png")

    fig, ax = plt.subplots(figsize=(6, 5))
    RocCurveDisplay.from_predictions(y_test, prob_gene, ax=ax, name="Gene-only")
    RocCurveDisplay.from_predictions(y_test, prob_ctx, ax=ax, name="Gene + spatial context")
    ax.set_title("ROC comparison")
    _save(fig, "fig12_roc_gene_vs_context.png")


def main() -> None:
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    export_exploration_figures()
    adata_clustered = sc.read_h5ad(CLUSTERED_H5)
    export_clustering_figures(adata_clustered)
    export_marker_figure(adata_clustered)
    adata_spatial = sc.read_h5ad(SPATIAL_H5)
    export_spatial_figures(adata_spatial)
    export_ml_figures(adata_spatial)
    print(f"Figures written to {FIG_DIR}")


if __name__ == "__main__":
    main()
