# `cell_types_multiomics`

Joint cell typing of an RNA + ATAC union in the sample-removed embedding (`Z_comp`):

1. **RNA clustering.** Leiden (resolution `cluster_resolution`, seed 42) on the RNA cells, using the first `num_PCs` columns of `use_rep`. Labels are `"1"`, `"2"`, … .
2. **Label transfer to ATAC.** k-nearest-neighbor graphs (`k_neighbors`, metric `transfer_metric`) within and across modalities are combined into a Jaccard-weighted shared-nearest-neighbor score; each ATAC cell takes the RNA label with the highest score. ATAC cells with no overlap at all (score 0, more common at small `k_neighbors`) are labelled **`"unassigned"`**. The score is stored in `obs['label_transfer_confidence']`.
3. **Modality-balance filter.** A cluster is **removed, with all its cells**, when either modality's share of it falls below 5% of that modality's overall share (for example, an RNA-only cluster). The `"unassigned"` group has no RNA cells, so it is always removed here. The function therefore returns a **new, smaller AnnData**.
4. Optional UMAP on all remaining cells, a CSV of labels, plots, and saving.

**Import:** `from sampledisco.preparation.multi_omics_cell_type_cpu import cell_types_multiomics`

**Source:** `preparation/multi_omics_cell_type_cpu.py:20`. GPU variant: `cell_types_multiomics_gpu` in `preparation/multi_omics_cell_type_gpu.py:19`, same signature.

## Signature

```python
def cell_types_multiomics(
    adata,
    modality_column="modality",
    rna_modality_value="RNA",
    atac_modality_value="ATAC",
    cell_type_column="cell_type",
    cluster_resolution=0.8,
    use_rep=None,
    num_PCs=50,
    k_neighbors=15,
    transfer_metric="cosine",
    compute_umap=True,
    save=False,
    output_dir=None,
    defined_output_path=None,
    verbose=True,
    generate_plots=True,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Union of RNA and ATAC cells with `obs[modality_column]` and a joint embedding in `.obsm`. |
| `modality_column` | str | `"modality"` | Modality column. |
| `rna_modality_value` / `atac_modality_value` | str | `"RNA"` / `"ATAC"` | Values marking the two modalities. Both must be present. |
| `cell_type_column` | str | `"cell_type"` | Output column for the labels. |
| `cluster_resolution` | float | `0.8` | Leiden resolution for the RNA cells. Wrapper key: `multiomics_cluster_resolution`. |
| `use_rep` | str or None | `None` | `.obsm` key to cluster and transfer in. `None` resolves `Z_comp`, then the legacy composition names (`Z_clust`, …), then `X_glue` as a last resort. Wrapper key: `multiomics_use_rep_celltype`. |
| `num_PCs` | int or None | `50` | Use only the first `num_PCs` columns of the embedding. |
| `k_neighbors` | int | `15` | Neighbors per cell in the transfer graphs. The wrapper passes `3`. |
| `transfer_metric` | str | `"cosine"` | Distance for the transfer graphs. Wrapper key: `multiomics_metric`. |
| `compute_umap` | bool | `True` | Neighbors + UMAP (`min_dist=0.5`) on all remaining cells. Wrapper key: `multiomics_generate_umap_celltype`. |
| `save` | bool | `False` | Write the AnnData to `{output_dir}/preprocess/adata_preprocessed.h5ad` (needs `output_dir`; skipped when `defined_output_path` is set). |
| `output_dir` | str or None | `None` | Root output directory for the CSV, plots and `save=True`. |
| `defined_output_path` | str or None | `None` | Write the AnnData to exactly this path. This happens **whenever it is set**, even with `save=False`. |
| `verbose` | bool | `True` | Print progress. |
| `generate_plots` | bool | `True` | Write the plots below (needs `output_dir`). |

## Returns

A **new** `AnnData`: the input minus the cells of modality-imbalanced clusters. It carries `obs[cell_type_column]` (categorical), `obs['label_transfer_confidence']` (ATAC cells; `NaN` for RNA), and `obsm['X_umap']` when `compute_umap=True`. Keep the return value: `adata = cell_types_multiomics(adata, ...)`.

## Output files

With `output_dir` set:

- `{output_dir}/preprocess/cell_type.csv`: columns `cell_id`, `modality`, `cell_type`.
- `{output_dir}/preprocess/visualization/` (with `generate_plots=True`): `label_transfer_confidence.png`, `celltype_modality_heatmap.png` and `celltype_modality_distribution.csv` (cells per cell type × modality), plus `umap_cell_type.png`, `umap_modality.png` and `umap_split_by_modality.png` when `obsm['X_umap']` exists.
- The AnnData at `defined_output_path`, or with `save=True` at `{output_dir}/preprocess/adata_preprocessed.h5ad`.

!!! note "Path change in 0.6.0"
    In 0.6.0, `save=True` without `defined_output_path` writes `{output_dir}/preprocess/adata_preprocessed.h5ad`. In 0.5.0 it wrote `{output_dir}/preprocess/preprocess/adata_preprocessed.h5ad` (a doubled folder). Pass `defined_output_path`, as below and as the wrapper does, to get the same path on both versions.

## Usage

Multi-omics chain, after [`harmonize_xglue`](harmonize_xglue.md). The union file is updated in place, as in the wrapper:

```python
import anndata as ad
import pandas as pd
from sampledisco.preparation.multi_omics_cell_type_cpu import cell_types_multiomics

union_path = "api_demo/multiomics/preprocess/adata_sample.h5ad"
mo = ad.read_h5ad(union_path)
print(mo.n_obs, "cells before")
mo = cell_types_multiomics(
    mo,
    modality_column="modality",
    cluster_resolution=0.8,
    k_neighbors=3,                 # what the wrapper uses
    output_dir="api_demo/multiomics",
    defined_output_path=union_path,
)
print(mo.n_obs, "cells after")
print(pd.crosstab(mo.obs["cell_type"], mo.obs["modality"]))
```

**Writes** → `api_demo/multiomics/preprocess/cell_type.csv`; under `api_demo/multiomics/preprocess/visualization/`: `umap_cell_type.png`, `umap_modality.png`, `umap_split_by_modality.png`, `label_transfer_confidence.png`, `celltype_modality_heatmap.png`, `celltype_modality_distribution.csv`; and the updated union file. The log reports 16 RNA clusters and a mean transfer confidence of 0.940; no cluster is imbalanced, so all 59,010 cells are kept. The crosstab (first and last rows):

```text
59010 cells after
modality   ATAC   RNA
cell_type
1          1436  1063
10         2848  1990
11         1854  1186
12         3093  5860
13         4453  5929
14         2729  3670
15           19   100
16            2    15
...
```

About 2.5 minutes (mostly the UMAP of 59,010 cells; 75 s on a faster node) and 2.8 GB.

!!! warning "Modality values other than RNA/ATAC on 0.5.0"
    In 0.5.0 the balance filter in step 3 always compared the values `"RNA"` and `"ATAC"`, whatever `rna_modality_value` / `atac_modality_value` said. With any other second modality, every cluster looked imbalanced and all cells were removed (6,000 of 6,000 on a demo subset relabelled `RNA` / `FLOW`). 0.6.0 uses the function's own modality values: the same `RNA` / `FLOW` subset with `atac_modality_value="FLOW"` now keeps all 6,000 cells. On 0.5.0, rename the second modality to `"ATAC"` before calling the function.

!!! note "Wrapper behaviour"
    The wrapper resolves `use_rep` (`multiomics_use_rep_celltype`, default `Z_comp`), passes `k_neighbors=3`, `save=True`, `output_dir=<multiomics output dir>` and `defined_output_path=<union file>`, and then copies the labels onto the per-modality files `preprocess/adata_rna_preprocessed.h5ad` / `adata_atac_preprocessed.h5ad` for DGE. It uses `cell_types_multiomics_gpu` only when `use_gpu` and `multiomics_use_gpu` are both true.
