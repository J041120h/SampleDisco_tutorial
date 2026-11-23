# `wrapper`

Runs the whole pipeline for RNA, ATAC and/or multi-omics from keyword arguments. The YAML config used by the CLI (`sampledisco -m complex --config config.yaml`) is passed to this function as `wrapper(**config)`, so every config key is a keyword argument here. For each enabled modality it runs two phases:

1. **Preparation**: preprocess → [`cell_types`](../rna/cell_types_linux.md) → [`compute_sample_embedding`](../shared/calculate_sample_embedding.md) (or [`run_autotune`](../shared/run_autotune.md) when `*_autotune_enable` is true). Multi-omics: [`multiomics_preparation`](../multiomics/multiomics_preparation.md) → [`harmonize_xglue`](../multiomics/harmonize_xglue.md) → [`cell_types_multiomics`](../multiomics/cell_types_multiomics_linux.md) → sample embedding → [`save_embedding_to_h5ad`](../shared/save_embedding_to_h5ad.md).
2. **Downstream**: [`build_sample_adata`](../shared/build_sample_adata.md), then [`downstream_analysis`](downstream_analysis.md) with the modality's step flags.

**Import:** `from sampledisco.wrapper.wrapper import wrapper`

**Source:** `wrapper/wrapper.py:538`

!!! warning "Import it from `sampledisco.wrapper.wrapper`"
    On 0.5.0, `from sampledisco import wrapper` gives the `sampledisco.wrapper` subpackage, which is not callable (`TypeError: 'module' object is not callable`). In 0.6.0 that subpackage is callable, so `sampledisco.wrapper(**cfg)` runs the pipeline, and `sampledisco.wrapper.wrapper.wrapper` is still this function. `from sampledisco.wrapper.wrapper import wrapper` works on both.

## Signature (abridged)

The function has 304 keyword arguments. All are optional except `output_dir`; the per-modality ones are prefixed `rna_`, `atac_` or `multiomics_`. Every argument, its default and the function that receives it are listed in the [config key reference](../../tutorials/config_reference.md).

```python
def wrapper(
    output_dir: str,
    # pipeline selection
    run_rna_pipeline: bool = True,
    run_atac_pipeline: bool = False,
    run_multiomics_pipeline: bool = False,
    # general
    use_gpu: bool = False,
    initialization: bool = True,
    verbose: bool = True,
    save_intermediate: bool = True,
    large_data_need_extra_memory: bool = False,
    continue_on_error: bool = False,
    random_state: int = 42,
    # rna_* (78 keys), atac_* (83 keys), multiomics_* (132 keys): see the config key reference
    ...
) -> Dict[str, Any]
```

!!! note "CLI vs Python"
    In Python every argument is optional except `output_dir`. The CLI is stricter: `sampledisco -m complex --config config.yaml` rejects unknown keys and requires every argument to appear in the YAML, except the three `*_save_cell_adata_after_embedding` keys. `sampledisco --init-config [PATH]` writes a complete template (the demo config) to start from.

## General arguments

| Name | Default | Description |
| --- | --- | --- |
| `output_dir` | — | Root output directory. Modalities write to `rna_output_dir` / `atac_output_dir` / `multiomics_output_dir`, which default to `{output_dir}/rna`, `{output_dir}/atac`, `{output_dir}/multiomics`. |
| `run_rna_pipeline` / `run_atac_pipeline` / `run_multiomics_pipeline` | `True` / `False` / `False` | Which modalities to run, in that order. |
| `use_gpu` | `False` | Use the RAPIDS code paths. On Linux, the wrapper first checks that `rmm`, `cupy` and `rapids_singlecell` import and sets up the RMM allocator; if that fails it logs a warning and runs on CPU. RAPIDS being installed does nothing unless this is `true`. The multi-omics branch additionally needs `multiomics_use_gpu` (default `true`). |
| `initialization` | `True` | `True`: start a fresh status file. `False`: load the step flags from an existing `{output_dir}/sys_log/main_process_status.json` (for bookkeeping; which steps run is still decided by the step flags you pass). |
| `verbose` | `True` | Print progress. |
| `save_intermediate` | `True` | Multi-omics only: re-save the union file after the `harmonize_xglue` step and store `uns['X_DR_sample']` in it after the sample embedding. |
| `large_data_need_extra_memory` | `False` | GPU only: RMM managed (unified) memory instead of a pool allocator. |
| `continue_on_error` | `False` | `False`: an exception in a modality is re-raised. `True`: it is recorded as `results['<modality>_error']` and the next modality runs. |
| `random_state` | `42` | Seed for the RNA/ATAC sample embedding and for RAISIN (all modalities). The multi-omics `harmonize_xglue` step uses `multiomics_random_state` instead. |

Required per modality: `rna_count_data_path` for RNA and `atac_count_data_path` for ATAC (both even when resuming with `*_preprocessing: false`); `multiomics_rna_file` and `multiomics_atac_file` for multi-omics (also when `multiomics_integration: false`).

## Returns

A `dict`:

| Key | Value |
| --- | --- |
| `status_flags` | `{"rna": {...}, "atac": {...}, "multiomics": {...}, "system_info": {...}}`: which steps completed. |
| `system_info` | `platform`, `python_version`, `use_gpu`, `gpu_available` (whether the GPU path was actually enabled). |
| `rna_results` / `atac_results` / `multiomics_results` | For each modality that ran: `adata` (the cell-level AnnData, or the multi-omics union, with `uns['X_DR_sample']`), `pseudo_adata` (the sample-level AnnData from `build_sample_adata`, or `None` if no embedding was derived) and `status_flags`. In 0.6.0, `rna_results` / `atac_results` also hold `celltype_col`: the cell-type column the embedding and downstream steps used (`"cell_type"` after cell typing). |
| `rna_error` / `atac_error` / `multiomics_error` | Only with `continue_on_error=True` and a failure: the error message. |

## Output files

Everything goes under the modality output directories; see the step pages, [`downstream_analysis`](downstream_analysis.md) and [Understanding the outputs](../../tutorials/outputs.md). The wrapper itself writes `{output_dir}/sys_log/main_process_status.json`. For a first full run from a config, see the [Quick start](../../tutorials/quickstart.md); for GPU setup, [GPU acceleration](../../tutorials/gpu.md).

## Usage

Continuing from the RNA chain: reuse the preprocessed file, re-run cell typing and the sample embedding, then sample distances and clustering. Trajectory analysis and the visualization step are switched off to keep the example short:

```python
from sampledisco.wrapper.wrapper import wrapper

results = wrapper(
    output_dir="api_demo/wrapper",
    run_rna_pipeline=True,
    rna_count_data_path="data/test_RNA.h5ad",                            # required even when resuming
    rna_preprocessing=False,                                             # reuse the file below
    rna_adata_cell_path="api_demo/rna/preprocess/adata_preprocessed.h5ad",
    rna_sample_level_batch_col="batch",
    rna_n_target_cell_clusters=8,
    rna_trajectory_analysis=False,
    rna_visualize_data=False,
    rna_cluster_number=2,
)
print(sorted(results))
print(sorted(results["rna_results"]))
print(results["rna_results"]["pseudo_adata"])
```

**Writes** → under `api_demo/wrapper/`:

- `sys_log/main_process_status.json`
- `rna/preprocess/adata_preprocessed.h5ad` (re-typed copy, with `uns['X_DR_sample']`), `rna/preprocess/cell_type.csv`
- `rna/sample_embedding/sample_embedding.csv`, `sample_embedding_blocks.npz`
- `rna/Sample_distance/cosine/` and `rna/Sample_distance/correlation/` (the default methods), five files each
- `rna/sample_cluster/kmeans_clusters_sample.csv`, `kmeans_sample_embedding.png`

The input file is only read. Printed at the end (about 18 s in total, 3.4 GB):

```text
Total execution time: 18.19 seconds (0.30 minutes)
['rna_results', 'status_flags', 'system_info']
['adata', 'celltype_col', 'pseudo_adata', 'status_flags']
AnnData object with n_obs × n_vars = 8 × 7
    obs: 'celltype', 'sample', 'old', 'large', 'batch', 'sev.level', 'n_genes', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'cell_type', 'cell_type_original', 'cluster_sample_kmeans'
    uns: 'X_DR_sample', 'sample_embedding_params'
    obsm: 'X_DR_sample'
```

!!! tip "The same run from a config file"
    `sampledisco --init-config` writes the template; `sampledisco -m complex --config config.yaml` runs it. From Python, `wrapper(**yaml.safe_load(open("config.yaml")))` runs the same config (without the CLI's key check).
