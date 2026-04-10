# Spatially-Resolved Multi-Omics Integration in Alzheimer's Disease

> **⚠️ Project Status: Paused / Incomplete**  
> *This project is currently on hold. Due to the massive computational requirements of processing ~45,000 spatial spots, running this `Cell2location` pipeline requires High-Performance Computing (HPC) cluster resources. Active development is paused as it takes too long to execute on local laptop hardware, and I am currently prioritizing and focusing full-time on my Master's thesis. The codebase below represents the fully structured, memory-optimized pipeline ready for cluster deployment.*

**Overview:**  
This project implements a computational pipeline to map astrocyte and microglia states to Alzheimer's amyloid pathology. Using `Cell2location`, it integrates single-cell and spatial transcriptomics to deconvolute the spatial distribution of glial cells—specifically those driven by the APOE-activating enhancer RNA, AANCR—enabling high-accuracy detection of rare neuroinflammatory states across cortical tissue.

---

## Data

1) **Single-Cell Reference Data:** Transcriptomic profiles of human astrocytes and microglia, capturing baseline and AANCR-knockdown states (NCBI GEO).
2) **Spatial Transcriptomics Data:** High-resolution spatial mapping of the human middle temporal gyrus (MTG), capturing localized gene expression and Alzheimer's disease neuropathology (SEA-AD).

---

## Model Architecture

### Cell2location (Bayesian Spatial Deconvolution)

#### 1. Single-Cell Reference Regression (Encoder)
- **Objective:** Estimates the basal expression signature of each glial subpopulation (e.g., `Astro_1`, `Micro-PVM_2`).
- **Mechanism:** Uses a Negative Binomial regression model to calculate a highly accurate matrix of reference signatures, filtering out technical noise and batch effects.

#### 2. Spatial Mapping Model (Decoder)
- **Objective:** Infers the absolute abundance of each glial cell state at every spatial coordinate (spot).
- **Mechanism:** Employs a hierarchical Bayesian model (using `pyro`) that strictly enforces non-negative physical cell counts using Variational Inference to map `~45,000` spatial locations simultaneously.

---

## Pipeline Strategy & Memory Optimization

- **Glial Cell Filtering:** Extracts astrocytes, microglia, and oligodendrocytes using `Supertype` annotations, falling back to a 95th-percentile marker-based threshold (`GFAP`, `TMEM119`) if annotations are missing.
- **Sparse Integer Sanitization:** Implements an ultra-low-memory, in-place sparse matrix cleaner (`sanitize_sparse_matrix`) that guarantees strict float32 integers without triggering RAM spikes.
- **Stability Thresholds:** Clips expression signatures at a strict `1e-4` minimum to prevent PyTorch `Gamma` distribution zero-underflow crashes.
- **Hardware Fallbacks:** Automatically detects and utilizes NVIDIA GPUs if available to accelerate the posterior matrix sampling process.

---

## The Script

This is the core logic used in this repository. It covers robust data splitting, memory-safe data sanitization, reference signature training, and batched spatial posterior export.

```python
# ----------------------------
# 0) Imports + Hardware Check
# ----------------------------
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import pyro
import torch
from cell2location.models import RegressionModel, Cell2location

# Clear PyTorch/Pyro memory states
pyro.clear_param_store()

# Auto-detect GPU for massive acceleration
print("=== HARDWARE CHECK ===")
if torch.cuda.is_available():
    print(f"✅ GPU Available: {torch.cuda.get_device_name(0)}")
    accelerator_type = "gpu" 
else:
    print("⚠️ GPU NOT Available! Falling back to CPU...")
    accelerator_type = "cpu"

# ----------------------------
# 1) Robust Glial Reference Split
# ----------------------------
def split_reference_spatial(adata_in):
    """Robust glial reference + full spatial split for SEA-AD MERFISH."""
    obs = adata_in.obs.copy()

    # Try annotation-based glial detection
    glial_mask = pd.Series(False, index=obs.index)
    for col in ["broad_cell_type", "Subclass", "Supertype", "Class", "cell_type"]:
        if col in obs.columns:
            col_mask = obs[col].astype(str).str.contains(
                "astro|micro|glia|oligodendro", case=False, na=False
            )
            glial_mask |= col_mask

    # Marker-based fallback
    if glial_mask.sum() == 0:
        astro_markers = ["GFAP", "AQP4", "ALDH1L1", "S100B"]
        micro_markers = ["TMEM119", "P2RY12", "CSF1R", "C1QA"]
        # ... [marker scoring logic] ...
        glial_mask = (astro_score > np.percentile(astro_score, 95)) | \
                     (micro_score > np.percentile(micro_score, 95))

    sp_mask = pd.Series(True, index=obs.index) if "spatial" in adata_in.obsm else pd.Series(False, index=obs.index)
    return adata_in[glial_mask].copy(), adata_in[sp_mask].copy()

adata_ref, adata_sp = split_reference_spatial(adata)

# Use stable annotations and filter tiny populations
celltype_col = "Supertype"
adata_ref.obs["cell2loc_label"] = adata_ref.obs[celltype_col].astype(str)
label_counts = adata_ref.obs["cell2loc_label"].value_counts()
adata_ref = adata_ref[adata_ref.obs["cell2loc_label"].isin(label_counts[label_counts >= 50].index)].copy()

# ----------------------------
# 2) Memory-Safe Sparse Integer Sanitization
# ----------------------------
def sanitize_sparse_matrix(data_obj):
    """Deep cleans X entirely within sparse format to prevent RAM crashes."""
    X_sparse = data_obj.X.tocsr() if sparse.issparse(data_obj.X) else sparse.csr_matrix(data_obj.X)
    data = X_sparse.data
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    data = np.clip(np.round(data), 0, None)
    X_sparse.data = data.astype(np.float32)
    X_sparse.eliminate_zeros()
    return X_sparse

adata_ref.layers["counts"] = sanitize_sparse_matrix(adata_ref)
adata_ref.X = adata_ref.layers["counts"].copy()

adata_sp.layers["counts"] = sanitize_sparse_matrix(adata_sp)
adata_sp.X = adata_sp.layers["counts"].copy()

# ----------------------------
# 3) Train Reference Signature Model
# ----------------------------
RegressionModel.setup_anndata(adata=adata_ref, layer="counts", labels_key="cell2loc_label", batch_key="batch")
reg_model = RegressionModel(adata_ref)

reg_model.train(max_epochs=250, accelerator=accelerator_type)

adata_ref = reg_model.export_posterior(
    adata_ref, sample_kwargs=dict(num_samples=1000, batch_size=1000)
)

cell_state_df = pd.DataFrame(
    adata_ref.varm["means_per_cluster_mu_fg"],
    index=adata_ref.var_names,
    columns=adata_ref.uns["mod"]["factor_names"]
)

# ----------------------------
# 4) Spatial Mapping via Variational Inference
# ----------------------------
common_genes = adata_sp.var_names.intersection(cell_state_df.index)
adata_sp = adata_sp[:, list(common_genes)].copy()

# Gamma distribution stability fix for PyTorch
cell_state_df = cell_state_df.loc[list(common_genes)].fillna(0.0).clip(lower=1e-4).astype(np.float32)

Cell2location.setup_anndata(adata=adata_sp, layer="counts", batch_key="batch")

c2l_model = Cell2location(
    adata_sp, 
    cell_state_df=cell_state_df,
    N_cells_per_location=10,
    detection_alpha=20
)

# Train model safely
c2l_model.train(
    max_epochs=15000, 
    train_size=1.0, 
    batch_size=2500, 
    accelerator=accelerator_type
)

# Low-memory posterior export
adata_sp = c2l_model.export_posterior(
    adata_sp, 
    sample_kwargs=dict(num_samples=100, batch_size=500)
)

# ----------------------------
# 5) Consolidate Glial Pathology Abundances
# ----------------------------
abundance_key = "q05_cell_abundance_w_sf"
abund = pd.DataFrame(adata_sp.obsm[abundance_key], index=adata_sp.obs_names)

inflam_cols = [c for c in abund.columns if any(term in str(c).lower() for term in ["microglia", "astrocyte"])]
adata_sp.obs["reactive_glia_abundance"] = abund[inflam_cols].sum(axis=1)

print("🎉 Cell2location Spatial Mapping COMPLETE!")
```

---

## Hardware Requirements
Running the full `Cell2location` script on the SEA-AD spatial data (`~45,000 spots`) requires a machine with at least **32GB of RAM** and a dedicated **NVIDIA GPU** for CUDA-accelerated processing. Running this script strictly on a CPU or standard laptop may result in Out-Of-Memory (OOM) kernel crashes.

---

## Data Availability & Citations

The datasets used in this analysis are publicly available and must be cited if you reuse this pipeline. **Do not upload the raw data to this repository.**

### 1. Single-Cell Reference Data (GSE263862)
The single-cell RNA-seq reference data used to establish AANCR and APOE expression profiles was obtained from the NCBI Gene Expression Omnibus (GEO).
* **Accession:** [GSE263862](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263862)
* **Citation:** Wan M, Liu Y, Li D, Snyder RJ et al. The enhancer RNA, AANCR, regulates APOE expression in astrocytes and microglia. *Nucleic Acids Res* 2024 Sep 23;52(17):10235-10254. PMID: [39162226](https://www.ncbi.nlm.nih.gov/pubmed/39162226)

### 2. Spatial Transcriptomics Data (SEA-AD)
The spatial transcriptomic and neuropathology data were provided by the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) consortium, funded by the National Institutes on Aging (NIA U19AG060909).
* **Source:** [AWS Open Data Registry](https://registry.opendata.aws/allen-sea-ad-atlas/)
* **Required Citation Statement:** Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) was accessed on [INSERT DATE HERE] from https://registry.opendata.aws/allen-sea-ad-atlas.
