# API Reference: Preparation

This page groups the preparation and orchestration entry points that run before sample-level downstream analysis.

## `wrapper(...)`

```python
wrapper(
    output_dir: str,
    run_rna_pipeline: bool = True,
    run_atac_pipeline: bool = False,
    run_multiomics_pipeline: bool = False,
    ...
) -> dict
```

The top-level orchestrator in `code/wrapper/wrapper.py`. It accepts a large YAML-compatible signature that fans out into RNA, ATAC, and multi-omics sub-pipelines before dispatching to `downstream_analysis(...)`.

**Key parameters**

- `output_dir`: common output root.
- `run_rna_pipeline`, `run_atac_pipeline`, `run_multiomics_pipeline`: modality selection flags.
- `use_gpu`, `initialization`, `verbose`, `save_intermediate`: global execution controls.

**Returns**

- A dictionary of pipeline results and status flags per modality.

**Example**

```python
import yaml
from wrapper.wrapper import wrapper

with open("config.yaml", "r", encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

results = wrapper(**config)
```

**See also**

- `rna_wrapper(...)`
- `atac_wrapper(...)`
- `multiomics_wrapper(...)`

## `rna_wrapper(...)`

```python
rna_wrapper(
    rna_count_data_path: str = None,
    rna_output_dir: str = None,
    preprocessing: bool = True,
    cell_type_cluster: bool = True,
    derive_sample_embedding: bool = True,
    cca_based_cell_resolution_selection: bool = False,
    ...
) -> dict
```

Runs RNA preprocessing, cell-type clustering, sample embedding, and optional CCA-based resolution selection.

**Key parameters**

- `rna_count_data_path`, `rna_output_dir`: required RNA input and output locations.
- `preprocessing`, `cell_type_cluster`, `derive_sample_embedding`: major step gates.
- `min_cells`, `min_genes`, `pct_mito_cutoff`: RNA QC controls.
- `leiden_cluster_resolution`, `existing_cell_types`: cell annotation controls.
- `sample_hvg_number`, `sample_embedding_dimension`: sample representation settings.

**Returns**

- `adata_cell`, `adata_sample`, `pseudo_adata`, and status flags.

**See also**

- `preprocess(...)` in `preparation/rna_preprocess_cpu.py`
- `cell_types(...)` in `preparation/cell_type_cpu.py`
- `calculate_sample_embedding(...)`

## `atac_wrapper(...)`

```python
atac_wrapper(
    atac_count_data_path: str = None,
    atac_output_dir: str = None,
    preprocessing: bool = True,
    cell_type_cluster: bool = True,
    derive_sample_embedding: bool = True,
    cca_based_cell_resolution_selection: bool = False,
    ...
) -> dict
```

Runs ATAC preprocessing, cell-type clustering, sample embedding, and optional resolution optimization.

**Key parameters**

- `min_features`, `max_features`, `min_cells_per_sample`: ATAC QC controls.
- `doublet_detection`, `tfidf_scale_factor`, `drop_first_lsi`: ATAC latent-space tuning.
- `cell_embedding_column`: reuse an existing embedding instead of auto-detecting.

**Returns**

- `adata_cell`, `adata_sample`, `pseudo_adata`, and status flags.

**See also**

- `preprocess(...)` in `preparation/atac_preprocess_cpu.py`
- `calculate_sample_embedding(...)`
- `find_optimal_cell_resolution(...)`

## `multiomics_wrapper(...)`

```python
multiomics_wrapper(
    rna_file=None,
    atac_file=None,
    multiomics_output_dir=None,
    integration=True,
    integration_preprocessing=True,
    dimensionality_reduction=True,
    find_optimal_resolution=False,
    ...
) -> Dict[str, Any]
```

Runs GLUE integration, integrated preprocessing, cell typing, optional resolution search, and multi-omics sample embedding.

**Key parameters**

- `rna_file`, `atac_file`, `multiomics_output_dir`: required inputs.
- `run_glue_preprocessing`, `run_glue_training`, `run_glue_gene_activity`: GLUE stage toggles.
- `cell_type_cluster`, `cluster_resolution`, `use_rep_celltype`: integrated cell typing controls.
- `optimization_target`, `coarse_start`, `coarse_end`, `fine_step`: resolution search controls.

**Returns**

- A dictionary containing integrated AnnData objects, pseudobulk outputs, and status flags.

**See also**

- `multiomics_preparation(...)`
- `integrate_preprocess(...)`
- `calculate_multiomics_sample_embedding(...)`

## Pipeline entry point

`code/SampleDisc.py` exposes the current command-line interface:

```python
parse_args()
load_config(config_path)
validate_config(config, func)
main()
```

Use `main()` through the CLI for the most stable user-facing execution path.
