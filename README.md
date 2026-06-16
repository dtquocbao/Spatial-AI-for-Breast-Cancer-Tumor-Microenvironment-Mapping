# Spatial AI for Breast Cancer Tumor Microenvironment Mapping

Exploratory analysis of **10x Visium** spatial transcriptomics data from human breast cancer tissue.

## Dataset

- **Sample:** `V1_Breast_Cancer_Block_A_Section_1`
- **Source:** [10x Genomics — Human Breast Cancer (Block A Section 1)](https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0)
- **Tissue:** Fresh frozen invasive ductal carcinoma (~3,798 spots under tissue)

## Progress

| Step | Notebook | Status | Output |
|------|----------|--------|--------|
| 1. Data exploration | `01_visium_breast_cancer_exploration.ipynb` | Done | Raw QC spatial plots |
| 2. QC & preprocessing | `02_qc_preprocessing.ipynb` | Done | `data/processed/breast_cancer_visium_qc.h5ad` |
| 3. Clustering | `03_clustering.ipynb` | Done | `data/processed/breast_cancer_visium_clustered.h5ad` |
| 4. Marker genes | `04_marker_genes.ipynb` | Done | DE tables & marker plots |
| 5. Spatial neighborhoods | `05_spatial_neighborhoods.ipynb` | Done | `data/processed/breast_cancer_spatial_analysis.h5ad` |
| 6. Baseline ML model | `06_baseline_ml_model.ipynb` | Planned | — |

## Project layout

```
data/
  raw/          # Downloaded 10x Visium files
  processed/    # AnnData objects between pipeline steps
notebooks/      # Analysis notebooks (run in order)
outputs/figures/
scripts/        # Batch scripts
```

## Setup (Anaconda)

Anaconda is installed on this machine, but `conda` may not be on your PATH in every terminal (e.g. Cursor, plain PowerShell). Use one of the options below.

### Option A — Anaconda Prompt (easiest)

Open **Anaconda Prompt** from the Start menu, then:

```bash
cd path\to\Spatial-AI-for-Breast-Cancer-Tumor-Microenvironment-Mapping
conda env create -f environment.yml
conda activate spatial-breast-cancer
```

### Option B — PowerShell with full path

```powershell
& "$env:USERPROFILE\anaconda3\Scripts\conda.exe" env create -f environment.yml
& "$env:USERPROFILE\anaconda3\Scripts\activate" spatial-breast-cancer
```

### Option C — Add conda to PATH permanently

Run once in **Anaconda Prompt**:

```bash
conda init powershell
```

Restart the terminal, then `conda activate spatial-breast-cancer` will work everywhere.

---

Register the Jupyter kernel (optional):

```bash
python -m ipykernel install --user --name spatial-breast-cancer --display-name "spatial-breast-cancer"
```

To update the environment after dependency changes:

```bash
conda env update -f environment.yml --prune
```

### Conda solver error (`libmamba` / `QueryFormat`)

If `conda install` fails with a libmamba solver error, use the classic solver:

```bash
conda install -c conda-forge <package> --solver=classic
```

Clustering requires **`leidenalg`** and **`python-igraph`** (both listed in `environment.yml`).

## Run

Run notebooks **in order** with the `spatial-breast-cancer` kernel:

```bash
jupyter notebook notebooks/
```

| Notebook | What it does |
|----------|--------------|
| `01_visium_breast_cancer_exploration.ipynb` | Download data, load with Scanpy, inspect AnnData, plot tissue image and per-spot QC |
| `02_qc_preprocessing.ipynb` | Filter spots/genes, normalize, find HVGs, PCA → save QC object |
| `03_clustering.ipynb` | Scale, PCA, neighbors, UMAP, Leiden clustering, resolution comparison, spatial cluster plots |
| `04_marker_genes.ipynb` | Wilcoxon DE per cluster, known TME marker spatial/UMAP/dot plots |
| `05_spatial_neighborhoods.ipynb` | Spatial graph, neighborhood enrichment, co-occurrence, centrality → save spatial object |
| `06_baseline_ml_model.ipynb` | Baseline machine learning on spot-level features *(planned)* |

**Script (batch, step 1 only):**

```bash
python scripts/01_visium_breast_cancer_qc.py
```

## Pipeline workflow

### Step 1 — Exploration (`01`)

1. Download Visium breast cancer data (Squidpy → 10x CDN) into `data/raw/`
2. Load with `scanpy.read_visium`
3. Inspect the `AnnData` object (`obs`, `var`, `obsm`, `uns['spatial']`)
4. Plot H&E tissue image, total UMI counts per spot, and genes detected per spot

### Step 2 — QC & preprocessing (`02`)

1. Load raw Visium data and compute QC metrics
2. Violin/scatter plots for `total_counts` and `n_genes_by_counts`
3. Filter spots (`min_counts=1000`) and genes (`min_cells=3`)
4. Normalize (`target_sum=1e4`) and log-transform
5. Select 3,000 highly variable genes
6. Run PCA and plot variance ratio
7. Save → `data/processed/breast_cancer_visium_qc.h5ad`

### Step 3 — Clustering (`03`)

1. Load preprocessed `h5ad` from step 2
2. Scale expression (`max_value=10`)
3. PCA (30 PCs used for graph)
4. Build kNN graph (`n_neighbors=10`)
5. UMAP embedding
6. Leiden clustering (`resolution=0.5`, key `leiden_0_5`); compare resolutions 0.2–1.0
7. Plot clusters on UMAP and H&E tissue image
8. Save → `data/processed/breast_cancer_visium_clustered.h5ad`

### Step 4 — Marker genes (`04`)

1. Load clustered `h5ad` from step 3
2. Run `sc.tl.rank_genes_groups` (Wilcoxon) per `leiden_0_5`
3. Visualize top markers per cluster
4. Score known biology markers (epithelial/tumor, immune, stromal, endothelial)
5. Dot plot of marker expression by cluster

### Step 5 — Spatial neighborhoods (`05`)

1. Load clustered `h5ad` from step 3
2. Build spatial neighbor graph (`sq.gr.spatial_neighbors`)
3. Visualize spatial connectivity on tissue
4. Neighborhood enrichment between Leiden clusters
5. Co-occurrence and centrality analysis
6. Save → `data/processed/breast_cancer_spatial_analysis.h5ad`

### Next step (planned)

- **06** — Baseline machine learning on spot-level features (cluster labels / TME phenotypes)
