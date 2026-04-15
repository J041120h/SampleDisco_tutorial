# `visualize_multimodal_embedding`

Flexible scatter-plot visualization of the multi-omics sample embeddings. When called with `show_default=True` (the default), renders side-by-side plots of `X_DR_expression` and `X_DR_proportion` with every sample visible. When `modality_col`, `color_col`, and `target_modality` are all provided, highlights one modality and colors by an arbitrary `.obs` column. `visualization_grouping_column` fans out additional colored panels — one per column — on top of whatever base layout you requested.

**Source:** `visualization/multi_omics_visualization.py:628`

## Signature

```python
def visualize_multimodal_embedding(
    adata,
    modality_col=None,
    color_col=None,
    target_modality=None,
    expression_key="X_DR_expression",
    proportion_key="X_DR_proportion",
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
| `adata` | AnnData | — | Sample-level AnnData with the two DR embeddings in `.obsm`. |
| `modality_col` | str, optional | `None` | `.obs` column marking modality (required for modality-specific plots). |
| `color_col` | str, optional | `None` | `.obs` column for point coloring (required when not using default plot). |
| `target_modality` | str, optional | `None` | Which modality to highlight in modality-specific plots. |
| `expression_key` | str | `"X_DR_expression"` | Key in `.obsm` for the expression embedding. |
| `proportion_key` | str | `"X_DR_proportion"` | Key in `.obsm` for the proportion embedding. |
| `figsize` | tuple | `(20, 8)` | Matplotlib figure size. |
| `point_size` | int | `60` | Scatter marker size. |
| `alpha` | float | `0.8` | Marker opacity. |
| `colormap` | str | `"viridis"` | Colormap for numerical `color_col` values. |
| `categorical_cmap` | str | `"tab10"` | Colormap for categorical `color_col` values. |
| `output_dir` | str, optional | `None` | File or directory. Each generated panel is saved under a descriptive filename. |
| `show_sample_names` | bool | `False` | Overlay sample IDs — best left off unless you have few samples. |
| `force_data_type` | `"numerical"` / `"categorical"` / `None` | `None` | Override the auto-detection for `color_col`. |
| `show_default` | bool | `True` | Always produce the default side-by-side view without requiring modality args. |
| `verbose` | bool | `True` | Print progress. |
| `visualization_grouping_column` | list, optional | `None` | Extra `.obs` columns; generates one colored panel per entry. |

## Returns

`Tuple[matplotlib.figure.Figure, numpy.ndarray]` — the figure and axes array when the plots are rendered, `None` when saved separately.

## Output files

Under `{output_dir}/`:

- `all_samples_expression.png`, `all_samples_proportion.png` — default panel pair.
- `all_samples_expression_by_{col}.png`, `all_samples_proportion_by_{col}.png` — per grouping column.
- `all_samples_combined_by_{col}.png` — combined side-by-side with shared colorbar.
- `{target_modality}_modality_split_by_{col}.png` — modality-specific view when all of `modality_col`, `color_col`, `target_modality` are supplied.

## Usage

```python
from genodistance.visualization import visualize_multimodal_embedding

visualize_multimodal_embedding(
    adata=pseudo_adata,
    modality_col="modality",
    target_modality="ATAC",
    color_col=None,
    expression_key="X_DR_expression",
    proportion_key="X_DR_proportion",
    visualization_grouping_column=["sev.level"],
    figsize=(20, 8),
    point_size=60,
    alpha=0.8,
    colormap="viridis",
    output_dir="/results/multiomics/visualization",
    show_default=True,
)
```
