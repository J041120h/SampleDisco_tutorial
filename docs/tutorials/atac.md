# ATAC

This page runs the scATAC-seq branch on the COVID demo, one function at a time: load the data, preprocess it, call cell types, build the sample embedding and plot it. The steps match the [RNA page](rna.md); only preprocessing differs. The arguments are the demo-config values, so you get the ATAC part of the [quick start](quickstart.md) run. It runs in about 3 minutes on 8 CPU cores and needs about 6 GB of RAM.

## 1. Load the data

Download the demo file (260 MB) and look at its samples; `preprocess` reads the file itself, so this step is only for a first look.

```bash
mkdir -p data
curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
```

```python
import anndata as ad

adata = ad.read_h5ad("data/test_ATAC.h5ad", backed="r")  # reads .obs; the counts stay on disk
print(adata.shape)
print(adata.obs.groupby("sample", observed=True)[["sev.level"]].first())
```

**Output**

- `data/test_ATAC.h5ad`: raw counts, 29,021 cells × 230,356 peaks from 8 samples (SRA runs), 4 healthy (`sev.level` `1`) and 4 severe COVID-19 (`4`). There is **no batch column**, so no batch key is used below.

## 2. Preprocess

`preprocess` filters cells and peaks, runs TF-IDF and LSI (dropping the first, depth-driven component), and then Harmony on `sample` to give `Z_comp` for cell typing; with no batch column, `Z_rmd` is the uncorrected LSI, which keeps the sample differences the embedding measures.

```python
from sampledisco.preparation.atac_preprocess_cpu import preprocess

OUT = "sampledisco_demo_output/atac"

adata = preprocess(
    h5ad_path="data/test_ATAC.h5ad",
    sample_meta_path=None,            # sample and sev.level are already in .obs
    output_dir=OUT,
    sample_column="sample",
    cell_level_batch_key=None,        # the ATAC demo has no batch column
    sample_level_batch_key=None,
    min_cells=10,
    min_features=200,
    max_features=50000,
    min_cells_per_sample=50,
    exclude_features=None,
    doublet_detection=False,
    num_cell_hvfs=10000,
    cell_embedding_num_PCs=30,
    num_harmony_iterations=20,
    tfidf_scale_factor=10000.0,
    log_transform=True,
    drop_first_lsi=True,
    verbose=True,
)
print(adata.shape)
print({k: adata.obsm[k].shape for k in sorted(adata.obsm.keys())})
```

**Output**

- `adata`, also written to `sampledisco_demo_output/atac/preprocess/adata_preprocessed.h5ad`: 29,021 cells × 228,367 peaks after QC, with `obsm['Z_comp']` and `obsm['Z_rmd']` (29 dimensions each: 30 LSI components minus the first).
- Scanpy's warning that `flavor='seurat_v3'` "expects raw count data" is expected: peaks are ranked on the TF-IDF matrix.

## 3. Cell types

`cell_types`, the same function as for RNA, runs Leiden on `Z_comp` (cosine distance, all 29 dimensions), brings the clusters to 8 cell types and draws a UMAP.

```python
from sampledisco.preparation.cell_type_cpu import cell_types

adata = cell_types(
    anndata_cell=adata,
    cell_type_column="cell_type",
    existing_cell_types=False,
    n_target_clusters=8,
    umap=True,
    save=True,
    output_dir=OUT,
    defined_output_path=f"{OUT}/preprocess/adata_preprocessed.h5ad",  # the default, written out
    leiden_cluster_resolution=0.8,
    cell_embedding_column=None,       # None = Z_comp
    cell_embedding_num_PCs=30,        # passed by the wrapper; not used for ATAC
    verbose=True,
)
print(adata.obs["cell_type"].value_counts().sort_index())
```

**Output**

- `adata.obs['cell_type']`: labels `"1"` to `"8"`. When Leiden finds more clusters, the closest ones are merged down to 8. The labels are also saved in `adata_preprocessed.h5ad` with the UMAP.
- `sampledisco_demo_output/atac/preprocess/cell_type.csv` (`cell_id`, `cell_type`) and `umap_cell_type.png` (below).

![UMAP of the demo ATAC cells colored by the 8 cell types](../resource/atac/umap_cell_type.png)
<div class="figure-caption">UMAP of the 29,021 cells, colored by the 8 cell types from this run. Types 7 and 8 hold only 46 and 16 cells, so they barely show. The cells form a continuum, so the borders between types move from run to run.</div>

## 4. Sample embedding

`compute_sample_embedding` turns each sample's cells into one row, exactly as for RNA; with `batch_col=None` there is no sample-level batch correction. Arguments not shown keep their defaults, which are also the demo-config values; [Sample embedding in depth](sample_embedding.md) explains the blocks and every option.

```python
from sampledisco.sample_embedding import compute_sample_embedding

adata = compute_sample_embedding(
    adata,
    OUT,
    sample_col="sample",
    celltype_col="cell_type",
    batch_col=None,                   # no batch column in the ATAC demo
    medium_K=120,
    fine_K=300,
    rmd_weight="equal",
    pca_components=10,
    seed=42,
)
emb = adata.uns["X_DR_sample"]
print(emb.shape)
```

**Output**

- `adata.uns['X_DR_sample']`: a DataFrame, 8 samples × 7 PCs (`pca_components` is capped at samples − 1), also saved in `adata_preprocessed.h5ad`.
- `sampledisco_demo_output/atac/sample_embedding/sample_embedding.csv` (the same table) and `sample_embedding_blocks.npz` (the blocks before weighting).

!!! note "Your numbers will differ a little"
    Reruns on the same kind of machine give identical results. On a different CPU type the ATAC cell-type borders can move noticeably, while the sample embedding stays close. See [Run-to-run variation](own_data.md#run-to-run-variation).

## 5. Plot the sample embedding

`build_sample_adata` turns `uns['X_DR_sample']` into a one-row-per-sample AnnData that keeps only the sample-level `obs` columns, such as `sev.level`, so the samples can be drawn on PC1 and PC2.

```python
import matplotlib.pyplot as plt
from sampledisco.sample_embedding.sample_embedding import build_sample_adata

sample_adata = build_sample_adata(adata, sample_col="sample")  # one row per sample
pcs, sev = sample_adata.obsm["X_DR_sample"], sample_adata.obs["sev.level"]
fig, ax = plt.subplots(figsize=(5, 4))
for level, color in zip(sorted(sev.unique()), ["#2a78d6", "#eb6834"]):
    ax.scatter(*pcs[sev == level, :2].T, s=60, color=color, edgecolor="white", label=level)
ax.set(xlabel="PC1", ylabel="PC2", title="Sample embedding")
ax.legend(title="sev.level")
fig.savefig(f"{OUT}/sample_embedding/sample_embedding_pc12.png", dpi=150, bbox_inches="tight")
```

**Output**

- `sample_adata`: 8 samples × 7 PCs in `obsm['X_DR_sample']`, with the sample-level `obs` columns; `build_sample_adata` prints the 23 per-cell columns it drops.
- `sampledisco_demo_output/atac/sample_embedding/sample_embedding_pc12.png` (below): the 8 samples on PC1 and PC2, colored by `sev.level`.

![The 8 demo ATAC samples on PC1 and PC2 of the sample embedding, colored by severity](../resource/atac/sample_embedding_pc12.png)
<div class="figure-caption">The 8 samples on PC1 and PC2 of this run's sample embedding, blue for healthy (1) and orange for severe (4). The four severe samples sit close together; the healthy samples spread out.</div>

For the package's built-in plots (a cell-type dendrogram and the PC1–PC2 scatter for each grouping column), see [Visualization](downstream/index.md#visualization) on the Downstream page.

Each argument above is an `atac_*` key of the demo config ([Configuration guide](configuration.md)); GPU versions of these calls are in [GPU acceleration](gpu.md).

**Next:** [Downstream analysis](downstream/index.md) runs the same analyses on this embedding; load `sampledisco_demo_output/atac/preprocess/adata_preprocessed.h5ad` instead of the RNA file.
