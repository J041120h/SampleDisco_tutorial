# API reference

The API is organized by the role each function plays in a SampleDisco run. Functions shared between modalities (cell typing, sample embedding, the downstream module) are documented once and cross-linked.

!!! tip "The config-driven wrapper is the supported entry point"
    Most users do not call these functions one by one. The primary interface is the YAML-config wrapper: `sampledisco -m complex --config <config.yaml>` (or `python -m sampledisco.cli -m complex --config <config.yaml>`), which calls [`wrapper(**config)`](pipeline/wrapper.md). For each enabled modality it runs preprocess → cell typing → sample embedding (or autotune) → [`downstream_analysis`](pipeline/downstream_analysis.md). The functions below are the building blocks it calls; the [config key reference](../tutorials/config_reference.md) maps every YAML key to its function argument.

## Where to import from

The import root is `sampledisco`. Only a few names are re-exported; everything else is imported from the module that defines it.

| Import | What you get |
| --- | --- |
| `from sampledisco.sample_embedding import compute_sample_embedding` | The CPU/GPU dispatcher (re-exported by the `sample_embedding` subpackage). |
| `from sampledisco.sample_prediction import predict_phenotype, transfer_labels` | Re-exported by the `sample_prediction` subpackage. |
| `import sampledisco` → `sampledisco.compute_sample_embedding`, `sampledisco.predict_phenotype`, `sampledisco.transfer_labels` | Lazy top-level exports; `from sampledisco import predict_phenotype` works. |
| `from sampledisco.wrapper.wrapper import wrapper` | The pipeline function. **Use this spelling**: it works on every release. |
| `from sampledisco import wrapper` | A callable package: `wrapper(**cfg)` / `sampledisco.wrapper(**cfg)` runs the pipeline, and `sampledisco.wrapper.wrapper.wrapper` is still the function. (On 0.5.0 this returned a non-callable module.) |
| Everything else | The defining module, as shown on each page, e.g. `from sampledisco.preparation.rna_preprocess_cpu import preprocess`. The `__init__.py` files of `preparation`, `parameter_selection`, `sample_distance`, `sample_association`, `sample_trajectory`, `sample_clustering`, `visualization` and `utils` are empty. |

## Sections

| Section | When you reach for it | Pages |
| --- | --- | --- |
| Pipeline | Run everything from one call, or re-run only the downstream steps. | [`wrapper`](pipeline/wrapper.md), [`downstream_analysis`](pipeline/downstream_analysis.md) |
| [RNA](rna.md) | scRNA-seq preprocessing and cell typing. | [`preprocess`](rna/preprocess_linux.md), [`cell_types`](rna/cell_types_linux.md) (shared with ATAC) |
| [ATAC](atac.md) | scATAC-seq preprocessing; cell typing uses the same `cell_types`. | [`preprocess`](atac/preprocess_linux.md) (ATAC), [`cell_types`](rna/cell_types_linux.md) |
| [Multi-omics](multiomics.md) | scGLUE integration, the `Z_comp` Harmony post-pass, joint cell typing, the multi-omics sample embedding. | [`multiomics_preparation`](multiomics/multiomics_preparation.md), [`harmonize_xglue`](multiomics/harmonize_xglue.md), [`cell_types_multiomics`](multiomics/cell_types_multiomics_linux.md), [`compute_sample_embedding` (multi-omics)](multiomics/calculate_multiomics_sample_embedding.md), [`replace_optimal_dimension_reduction`](multiomics/replace_optimal_dimension_reduction.md) (legacy) |
| Shared | The sample embedding itself, used by every modality. | [`compute_sample_embedding`](shared/calculate_sample_embedding.md), [`run_autotune`](shared/run_autotune.md), [`build_sample_adata`](shared/build_sample_adata.md), [`save_embedding_to_h5ad`](shared/save_embedding_to_h5ad.md) |
| [Downstream](downstream.md) | Everything that runs on the sample embedding. | `sample_distance`, `run_dimension_association_analysis`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `raisinfit`, `run_pairwise_tests`, `predict_phenotype`, `transfer_labels`, `visualization`, `visualize_multimodal_embedding` |
| Utilities | Harmony, embedding key names, inspecting an h5ad. | [`harmonize_embedding`](utils/harmonize_embedding.md), [`migrate_obsm_keys` and the key resolvers](utils/embedding_keys.md), [`summarize_h5ad`](utils/summarize_h5ad.md) |

## Running the examples

Every **Usage** block on these pages was run as written on the COVID-19 demo (8 samples: 4 healthy, 4 severe; see [Demo data](../tutorials/demo_data.md)) with sampledisco 0.6.1 on 8 CPU cores. Put the demo files in `./data/`:

```bash
mkdir -p data
curl -fL -o data/test_RNA.h5ad "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
```

The examples write under `api_demo/` and form three chains. Run each chain in this order; within a chain, a later example reads the file an earlier one wrote.

1. **RNA:** [`preprocess`](rna/preprocess_linux.md) → [`cell_types`](rna/cell_types_linux.md) → [`compute_sample_embedding`](shared/calculate_sample_embedding.md). This leaves `api_demo/rna/preprocess/adata_preprocessed.h5ad` with cell types and `uns['X_DR_sample']`. Every RNA-based example after that (autotune, downstream, pipeline, utilities) starts by reloading it:

    ```python
    import anndata as ad
    from sampledisco.sample_embedding.sample_embedding import build_sample_adata

    adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")  # cells
    pseudo = build_sample_adata(adata, sample_col="sample")                # one row per sample
    ```

2. **ATAC:** [`preprocess` (ATAC)](atac/preprocess_linux.md) → the ATAC example on [`cell_types`](rna/cell_types_linux.md).
3. **Multi-omics** (starts from the pre-integrated file, so no scGLUE is needed): [`harmonize_xglue`](multiomics/harmonize_xglue.md) → [`cell_types_multiomics`](multiomics/cell_types_multiomics_linux.md) → [`compute_sample_embedding` (multi-omics)](multiomics/calculate_multiomics_sample_embedding.md) → [`save_embedding_to_h5ad`](shared/save_embedding_to_h5ad.md) → [`visualize_multimodal_embedding`](downstream/visualize_multimodal_embedding.md). The [key-name example](utils/embedding_keys.md) also reads the pre-integrated file.

Each example lists what it writes, what it printed, and its runtime and peak RAM. Examples that call `build_sample_adata` first print its `[build_sample_adata] dropped … obs column(s) …` line; the pages leave it out. On 8 CPU cores, running every RNA-based example in order took about 7 minutes (peak 4.6 GB, during RAISIN), the two ATAC examples about 6.5 minutes (8.7 GB), and the multi-omics chain about 3.5 minutes (3 GB).

!!! note "Your numbers will differ"
    Rerunning on the same kind of machine gives identical results, but different CPU types (e.g. AVX2 vs AVX-512) give slightly different floating-point results, and Leiden can then split the cells differently. Two runs of the RNA chain on two cluster nodes found 19 and 18 Leiden clusters before the merge to 8, and every number downstream moved with them (the CCA correlation was 0.762 and 0.774). The Harmony backend matters too: these runs used harmonypy, while a default install uses harmony-pytorch. File names, shapes and the steps that run do not change. The printed outputs on these pages come from one complete run.

!!! warning "Eight samples are a smoke test, not an analysis"
    The demo has 8 samples (16 units in multi-omics). The sample embedding is capped at 7 PCs (units − 1); a CCA on all 7 PCs has a correlation of 1.0 by construction and a meaningless permutation p-value (the examples use 2 PCs); p-values and prediction scores are noisy; and trajectory DGE has very little power. Use the examples to learn the calls and the outputs; see [Small cohorts](../tutorials/downstream/index.md#small-cohorts) for the limits.

If you are unsure which function to call next, the tutorials walk through the end-to-end order: [RNA](../tutorials/rna.md), [ATAC](../tutorials/atac.md), [Multi-omics](../tutorials/multiomics.md), [Downstream](../tutorials/downstream/index.md).
