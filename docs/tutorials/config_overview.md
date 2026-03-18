# Overview: Using Config Files

SampleDisc is currently easiest to run through a YAML configuration file passed into the main pipeline wrapper.

## Where the config is used

The entry point in `code/SampleDisc.py` loads YAML with `load_config(...)`, validates keys against `wrapper(...)`, and then executes the pipeline in complex mode.

```bash
python /users/hjiang/GenoDistance/code/SampleDisc.py -m complex \
  --config /users/hjiang/GenoDistance/code/config/config_covid_rna.yaml
```

## Recommended mental model

A config file is easiest to read in four blocks:

| Section | What it controls | Representative keys |
| --- | --- | --- |
| Pipeline selection | Which modality pipelines run | `run_rna_pipeline`, `run_atac_pipeline`, `run_multiomics_pipeline` |
| Input and output | Raw files, metadata, resume paths, output folders | `output_dir`, `rna_count_data_path`, `rna_sample_meta_path` |
| Representation learning | Preprocessing, clustering, sample embedding, resolution search | `rna_min_cells`, `rna_leiden_cluster_resolution`, `rna_sample_embedding_dimension` |
| Downstream analysis | Distances, trajectories, clustering, DGE, plots | `rna_trajectory_analysis`, `rna_sample_distance_methods`, `rna_visualize_data` |

## Example excerpt from `config_covid_rna.yaml`

```yaml
output_dir: "/dcs07/hongkai/data/harry/result/long_covid"

run_rna_pipeline: true
run_atac_pipeline: false
run_multiomics_pipeline: false

rna_preprocessing: false
rna_cell_type_cluster: false
rna_derive_sample_embedding: true
rna_trajectory_analysis: true
rna_trajectory_dge: true

rna_sample_col: "sample"
rna_batch_col: "batch"
rna_celltype_col: "cell_type"

rna_sample_hvg_number: 2000
rna_sample_embedding_dimension: 10
rna_n_cca_pcs: 10
rna_trajectory_col: "month"
```

## Key parameter groups

### 1. Input and output paths

- `output_dir`: common output root for the whole run.
- `rna_count_data_path`, `atac_count_data_path`, `multiomics_rna_file`, `multiomics_atac_file`: primary data inputs.
- `rna_adata_cell_path`, `rna_adata_sample_path`, `rna_pseudo_adata_path`: resume from saved intermediate files when preprocessing is disabled.

### 2. Preprocessing options

- RNA QC is driven by keys such as `rna_min_cells`, `rna_min_genes`, and `rna_pct_mito_cutoff`.
- ATAC QC adds feature filtering and optional doublet detection through `atac_min_features`, `atac_max_features`, and `atac_doublet_detection`.
- Multi-omics preprocessing includes integration-specific gates such as `multiomics_min_cells_sample`, `multiomics_min_cell_gene`, and `multiomics_pct_mito_cutoff`.

### 3. Embedding and representation parameters

- `rna_sample_hvg_number` and `atac_sample_hvg_number` control pseudobulk feature selection.
- `rna_sample_embedding_dimension`, `atac_sample_embedding_dimension`, `multiomics_n_expression_components`, and `multiomics_n_proportion_components` define the reduced sample space.
- `*_harmony_for_proportion` toggles Harmony correction on proportion-derived embeddings.

### 4. Downstream analysis flags

- `*_sample_distance_calculation` controls distance-matrix generation.
- `*_trajectory_analysis` and `*_trajectory_dge` enable supervised CCA or unsupervised TSCAN plus GAM-based testing.
- `*_sample_cluster`, `*_proportion_test`, and `*_cluster_dge` gate clustering and differential modules.

## Output expectations

A representative RNA run can produce:

- `preprocess/` intermediates such as `adata_cell.h5ad` and `adata_sample.h5ad`
- `pseudobulk/` sample-level embedding objects
- `Sample_distance/` matrices and low-dimensional coordinates
- `CCA/` or `TSCAN/` pseudotime outputs
- `trajectoryDEG/` significant trajectory-linked features
- `sample_cluster/` clustering assignments

## Runtime guidance

!!! note
    Runtime is highly dataset- and hardware-dependent. A config that reuses preprocessed `.h5ad` files is much faster than a full raw-data run.

- RNA-only configs often complete in **minutes to under an hour** once preprocessing is skipped.
- ATAC-only configs usually require **longer QC and dimensionality reduction** than RNA.
- Multi-omics configs with GLUE training can take **hours on large cohorts**, and benefit strongly from GPU execution.

## Common pitfalls

!!! warning
    `validate_config(...)` checks config keys strictly against the `wrapper(...)` signature. A typo in a key name will fail early before any analysis begins.

!!! warning
    The repository currently mixes CLI execution and module imports rather than exposing a stable `import SampleDisc as sd` API. Prefer the documented config-driven entry point until the package namespace is finalized.
