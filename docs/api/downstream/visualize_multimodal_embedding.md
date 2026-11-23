# `visualize_multimodal_embedding`

Scatter plots (first two PCs) of the multi-omics sample embedding, one point per `<sample>_<modality>` unit. It has two modes:

- **Default mode**, used when `show_default=True` (the default) or when `modality_col`, `color_col` and `target_modality` are all `None`. Without `visualization_grouping_column`, all units are drawn in one color. With `visualization_grouping_column`, the plain plots are **replaced** by one set of plots per column, colored by that column (categorical or numeric, detected automatically).
- **Modality mode**, used only with `show_default=False` and all three of `modality_col`, `color_col`, `target_modality`: units of `target_modality` are colored by `color_col`, the other modality is drawn in grey. `visualization_grouping_column` is not used here.

`expression_key` and `proportion_key` are leftovers of the old two-embedding layout. Both default to `X_DR_sample`, so the "Expression" and "Proportion" panels show the **same** embedding and differ only in title. The key is looked up in `.obsm` first, then `.uns`.

**Import:** `from sampledisco.visualization.multi_omics_visualization import visualize_multimodal_embedding`

**Source:** `visualization/multi_omics_visualization.py:512`

## Signature

```python
def visualize_multimodal_embedding(
    adata,
    modality_col=None,
    color_col=None,
    target_modality=None,
    expression_key="X_DR_sample",
    proportion_key="X_DR_sample",
    figsize=(20, 8),
    point_size=60,
    alpha=0.8,
    colormap="viridis",
    categorical_cmap="tab10",
    output_dir=None,
    show_sample_names=False,
    force_data_type=None,
    show_default=True,
    verbose=True,
    visualization_grouping_column=None,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Sample-level AnnData of the units, e.g. [`build_sample_adata(mo, modality_col="modality")`](../shared/build_sample_adata.md). |
| `modality_col` | str or None | `None` | Modality column (modality mode). |
| `color_col` | str or None | `None` | Column to color `target_modality` by (modality mode). |
| `target_modality` | str or None | `None` | Modality to highlight (modality mode). |
| `expression_key`, `proportion_key` | str | `"X_DR_sample"` | Embedding keys of the two panels; leave both at the default. |
| `figsize` | tuple | `(20, 8)` | Size of side-by-side figures. |
| `point_size`, `alpha` | int, float | `60`, `0.8` | Marker size and opacity. |
| `colormap` | str | `"viridis"` | Colormap for numeric colors. |
| `categorical_cmap` | str | `"tab10"` | Colormap for categorical colors (default mode). |
| `output_dir` | str or None | `None` | A **folder** (an existing directory or a path without extension) becomes `{output_dir}/visualization/`; a path with an extension is used as a file name base. `None`: nothing is saved and figures are returned. |
| `show_sample_names` | bool | `False` | Write unit ids next to the points. |
| `force_data_type` | str or None | `None` | `"numerical"` or `"categorical"` to override the detection. |
| `show_default` | bool | `True` | Default mode; set `False` for modality mode. |
| `verbose` | bool | `True` | Print progress. |
| `visualization_grouping_column` | str, list or None | `None` | Default mode: one set of colored plots per column. Missing columns raise `ValueError`. |

## Returns

With a folder `output_dir`: `(None, None)`. Without `output_dir`: `(fig, axes)`; in default mode with several grouping columns, a dict `{column: (fig, axes)}`.

## Output files

With `output_dir` a folder, under `{output_dir}/visualization/`:

| Mode | Files |
| --- | --- |
| default, no grouping columns | `all_samples_expression.png`, `all_samples_proportion.png` |
| default, with grouping columns | per column: `all_samples_expression_by_{col}.png`, `all_samples_proportion_by_{col}.png`, `all_samples_combined_by_{col}.png` |
| modality mode | `{target_modality}_{color_col}_expression.png`, `{target_modality}_{color_col}_proportion.png` |

## Usage

Last step of the multi-omics chain (after [`save_embedding_to_h5ad`](../shared/save_embedding_to_h5ad.md) stored the embedding in the union file):

```python
import anndata as ad
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.visualization.multi_omics_visualization import visualize_multimodal_embedding

mo = ad.read_h5ad("api_demo/multiomics/preprocess/adata_sample.h5ad")
pseudo_mo = build_sample_adata(mo, sample_col="sample", modality_col="modality")

# default mode, colored by severity and by modality
visualize_multimodal_embedding(
    pseudo_mo,
    visualization_grouping_column=["sev.level", "modality"],
    output_dir="api_demo/multiomics",
)
# modality mode: ATAC units colored by severity, RNA units in grey
visualize_multimodal_embedding(
    pseudo_mo,
    modality_col="modality",
    color_col="sev.level",
    target_modality="ATAC",
    show_default=False,
    output_dir="api_demo/multiomics",
)
```

**Writes** → under `api_demo/multiomics/visualization/` (about 3 s in total):

- default mode: `all_samples_expression_by_sev.level.png`, `all_samples_proportion_by_sev.level.png`, `all_samples_combined_by_sev.level.png`, and the same three `*_by_modality.png`;
- modality mode: `ATAC_sev.level_expression.png`, `ATAC_sev.level_proportion.png`.

`sev.level` is detected as numeric (2 values) and `modality` as categorical. Both calls return `(None, None)` because `output_dir` is set.

!!! note "Wrapper behaviour"
    With `multiomics_visualize_embedding: true`, `downstream_analysis` calls this function twice, both times with `output_dir=<multiomics output dir>`:
    
    - in default mode, with `visualization_grouping_column=multiomics_visualization_grouping_column`;
    - when `multiomics_target_modality` is set (default `"ATAC"`), again in modality mode. The colour column is `multiomics_color_col`, or, if that is unset, the first of `multiomics_visualization_grouping_column` / `multiomics_grouping_columns` present in `.obs`.
