# Configuration Guide

SampleDisc is designed to run from a YAML configuration file and the main wrapper pipeline.

## Run with config

=== "CLI"

    ```bash
    python /users/hjiang/GenoDistance/code/SampleDisc.py -m complex \
      --config /users/hjiang/GenoDistance/code/config/config_covid_rna.yaml
    ```

=== "Python"

    ```python
    import yaml
    from wrapper.wrapper import wrapper

    with open("/users/hjiang/GenoDistance/code/config/config_covid_rna.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    wrapper(**config)
    ```

## How to read the config

Use four blocks as your mental model:

| Block | What it controls | Example keys |
| --- | --- | --- |
| Pipeline selection | Which wrappers run | `run_rna_pipeline`, `run_atac_pipeline`, `run_multiomics_pipeline` |
| Input/output | Paths and resume files | `output_dir`, `rna_count_data_path`, `multiomics_rna_file` |
| Representation learning | Preprocessing, cell types, sample embedding | `rna_preprocessing`, `atac_leiden_cluster_resolution`, `multiomics_n_expression_components` |
| Downstream analysis | Distance, trajectory, DGE, clustering, visualization | `rna_sample_distance_calculation`, `atac_trajectory_analysis`, `multiomics_visualize_embedding` |

## Global control parameters

These parameters are shared controls and are not specific to one modality API page.

| Parameter | Description |
| --- | --- |
| `output_dir` | Root output directory. |
| `run_rna_pipeline` | Enable RNA wrapper branch. |
| `run_atac_pipeline` | Enable ATAC wrapper branch. |
| `run_multiomics_pipeline` | Enable multi-omics wrapper branch. |
| `use_gpu` | Use GPU implementations where available. |
| `initialization` | Reset existing result folders before running. |
| `verbose` | Print detailed runtime logs. |
| `save_intermediate` | Save intermediate `.h5ad` outputs. |
| `large_data_need_extra_memory` | Use memory-focused behavior for large datasets. |

## Parameter naming rule

Config keys map directly to wrapper parameters by prefix:

- `rna_*` -> `rna_wrapper(...)` parameters and RNA `downstream_analysis(...)` parameters
- `atac_*` -> `atac_wrapper(...)` parameters and ATAC `downstream_analysis(...)` parameters
- `multiomics_*` -> `multiomics_wrapper(...)` parameters and multi-omics `downstream_analysis(...)` parameters

Example:

- `rna_preprocessing` maps to `rna_wrapper(preprocessing=...)`
- `atac_sample_distance_calculation` maps to `downstream_analysis(sample_distance_calculation=...)` for ATAC
- `multiomics_color_col` maps to `downstream_analysis(multiomics_color_col=...)`

## Core config example

```yaml
output_dir: "/dcs07/hongkai/data/harry/result/test"
run_rna_pipeline: true
run_atac_pipeline: false
run_multiomics_pipeline: false
use_gpu: true

rna_count_data_path: "/path/to/rna.h5ad"
rna_sample_meta_path: "/path/to/sample_meta.csv"
rna_preprocessing: true
rna_cell_type_cluster: true
rna_derive_sample_embedding: true
rna_sample_distance_calculation: true
rna_trajectory_analysis: true
rna_trajectory_dge: true
```

## Output structure

Typical outputs appear under `{output_dir}/{modality}/`:

- `preprocess/`
- `pseudobulk/`
- `embeddings/`
- `Sample_distance/`
- `CCA/` or `TSCAN/`
- `trajectoryDEG/`
- `sample_cluster/`
- `visualization/` (multi-omics uses this for embedding visualization)

## Common pitfalls

!!! warning
    `SampleDisc.py` validates config keys against the `wrapper(...)` signature. A typo in a key causes an immediate failure.

!!! warning
    If preprocessing is disabled, resume paths (`*_adata_cell_path`, `*_adata_sample_path`, `*_pseudo_adata_path`) must point to compatible files.

