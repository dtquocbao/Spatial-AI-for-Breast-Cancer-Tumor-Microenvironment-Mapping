# Spatial AI for Breast Cancer Tumor Microenvironment Mapping

Exploratory analysis of **10x Visium** spatial transcriptomics data from human breast cancer tissue.

## Dataset

- **Sample:** `V1_Breast_Cancer_Block_A_Section_1`
- **Source:** [10x Genomics — Human Breast Cancer (Block A Section 1)](https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0)
- **Tissue:** Fresh frozen invasive ductal carcinoma (~3,798 spots under tissue)

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

## Run

**Notebook (interactive):**

```bash
jupyter notebook notebooks/01_visium_breast_cancer_exploration.ipynb
```

**Script (batch):**

```bash
python scripts/01_visium_breast_cancer_qc.py
```

Downloads are cached under `data/`. Figures are written to `outputs/figures/`.

## Workflow

1. Download Visium breast cancer data (Squidpy helper → 10x CDN)
2. Load with `scanpy.read_visium`
3. Inspect the `AnnData` object (`obs`, `var`, `obsm`, `uns['spatial']`)
4. Plot H&E tissue image
5. Plot total UMI counts per spot
6. Plot number of detected genes per spot
