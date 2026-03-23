# `visualize_multimodal_embedding`

**Module**: `code/visualization/multi_omics_visualization.py`

```python
visualize_multimodal_embedding(
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

Visualizes multi-omics embeddings with modality-aware coloring.

## Parameters

| Parameter | Default | Controls |
| --- | --- | --- |
| `adata` | required | AnnData containing embeddings/metadata. |
| `modality_col` | `None` | Modality column in `obs`. |
| `color_col` | `None` | Metadata column used for coloring points. |
| `target_modality` | `None` | Modality to highlight for modality-specific views. |
| `expression_key` | `"X_DR_expression"` | Expression embedding key. |
| `proportion_key` | `"X_DR_proportion"` | Proportion embedding key. |
| `figsize` | `(20, 8)` | Figure size. |
| `point_size` | `60` | Scatter point size. |
| `alpha` | `0.8` | Point transparency. |
| `colormap` | `"viridis"` | Colormap for continuous data. |
| `categorical_cmap` | `"tab10"` | Colormap for categorical labels. |
| `output_dir` | `None` | Directory/path for saved visualizations. |
| `show_sample_names` | `False` | Draw sample labels on plot points. |
| `force_data_type` | `None` | Force coloring mode (`numerical`/`categorical`). |
| `show_default` | `True` | Create default all-sample panel. |
| `verbose` | `True` | Print progress logs. |
| `visualization_grouping_column` | `None` | Additional columns for separate coloring panels. |

