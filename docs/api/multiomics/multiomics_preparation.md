# `multiomics_preparation`

Full scGLUE integration pipeline for unpaired (or paired) scRNA + scATAC data. Internally runs five sub-stages, each individually toggleable via `run_preprocessing`, `run_training`, `run_merge`, `run_preprocess_per_modality`, `run_visualization`:

1. **scGLUE preprocessing** — reads both modalities, merges sample metadata, selects HVGs on RNA, runs LSI on ATAC, and builds a gene-region guidance graph from an Ensembl GTF.
2. **Adversarial scGLUE training** — fits the joint embedding with a consistency-regularized adversarial loss. The primary run produces the sample-preserved embedding (the paper's `Z_rmd`). With `run_second_glue_for_sample_removal=True`, a second run with `treat_sample_as_batch=True` produces the sample-removed embedding (`Z_clust`).
3. **Cell-union merge** — builds an embedding-only union AnnData (`preprocess/adata_sample.h5ad`) carrying `obsm['Z_rmd']` (primary `X_glue`, aliased) and, when the second run is enabled, `obsm['Z_clust']`. No expression matrix is stored — see `multi_omics_merge.py`.
4. **Per-modality preprocess** — per-modality QC + normalize (RNA) and QC + TF-IDF (ATAC), writing `preprocess/adata_rna_preprocessed.h5ad` and `preprocess/adata_atac_preprocessed.h5ad` for downstream DGE / RAISIN.
5. **Visualization (optional)** — UMAPs/scatter on the joint space, colored by modality or user-specified columns.

The returned object is the merged union AnnData with the integrated embedding in `.obsm`, ready for [`cell_types_multiomics`](cell_types_multiomics_linux.md). The single sample embedding for downstream analysis is produced later by `compute_sample_embedding` and stored in `uns['X_DR_sample']`.

!!! note "Gene-activity step removed"
    Earlier versions ran a fourth sub-stage that imputed per-cell gene activity by KNN on the joint embedding. That step has been removed; the union is now built directly from the per-modality scGLUE embeddings (`build_embedding_union`). The old `run_gene_activity`, `k_neighbors`, `use_rep`, `metric`, and `use_gpu` arguments no longer exist.

**Source:** `preparation/multi_omics_glue.py:801`

## Signature

```python
def multiomics_preparation(
    # Data files
    rna_file: str,
    atac_file: str,
    rna_sample_meta_file: Optional[str] = None,
    atac_sample_meta_file: Optional[str] = None,
    additional_hvg_file: Optional[str] = None,
    # Process control flags
    run_preprocessing: bool = True,
    run_training: bool = True,
    run_merge: bool = True,
    run_preprocess_per_modality: bool = True,
    run_visualization: bool = True,
    # Preprocessing parameters
    ensembl_release: int = 98,
    species: str = "homo_sapiens",
    use_highly_variable: bool = True,
    n_top_genes: int = 2000,
    n_top_peaks: int = 50000,
    atac_min_cells_floor: int = 10,
    n_pca_comps: int = 50,
    n_lsi_comps: int = 50,
    gtf_by: str = "gene_name",
    flavor: str = "seurat_v3",
    generate_umap: bool = False,
    rna_sample_column: str = "sample",
    atac_sample_column: str = "sample",
    # Training parameters
    consistency_threshold: float = 0.05,
    treat_sample_as_batch: bool = False,
    save_prefix: str = "glue",
    batch_key: Optional[str] = None,
    sample_key: str = "sample",
    data_batch_size: int = 1024,
    max_epochs: Optional[int] = None,
    dataloader_num_workers: int = 0,
    dataloader_fetches_per_worker: int = 4,
    array_shuffle_num_workers: int = 0,
    graph_shuffle_num_workers: int = 0,
    run_second_glue_for_sample_removal: bool = False,
    second_run_save_prefix: str = "glue_no_sample",
    # Per-modality preprocess QC params
    rna_min_cells: int = 500,
    rna_min_genes: int = 500,
    rna_pct_mito_cutoff: float = 20.0,
    rna_exclude_genes: Optional[List[str]] = None,
    atac_min_cells: int = 1,
    atac_min_features: int = 2000,
    atac_max_features: int = 15000,
    atac_min_cells_per_sample: int = 1,
    atac_exclude_features: Optional[List[str]] = None,
    atac_doublet_detection: bool = True,
    atac_tfidf_scale_factor: float = 1e4,
    atac_log_transform: bool = True,
    verbose: bool = True,
    # Visualization parameters
    plot_columns: Optional[List[str]] = None,
    # Output directory
    output_dir: str = "./glue_results",
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `rna_file`, `atac_file` | str | — | Input `.h5ad` files. |
| `rna_sample_meta_file`, `atac_sample_meta_file` | str, optional | `None` | Per-modality sample metadata CSVs. |
| `additional_hvg_file` | str, optional | `None` | Plain-text list of extra HVGs forced into the RNA HVG set. |
| `run_preprocessing` / `_training` / `_merge` / `_preprocess_per_modality` / `_visualization` | bool | `True` | Toggle individual sub-stages. |
| `ensembl_release` | int | `98` | Ensembl release used for GTF retrieval. |
| `species` | str | `"homo_sapiens"` | Species key matching Ensembl. |
| `use_highly_variable` | bool | `True` | Restrict scGLUE training to HVGs. |
| `n_top_genes` | int | `2000` | HVG count for RNA. |
| `n_top_peaks` | int | `50000` | Top peak count for ATAC. |
| `atac_min_cells_floor` | int | `10` | Floor on the per-feature min-cells filter for ATAC. |
| `n_pca_comps` / `n_lsi_comps` | int | `50 / 50` | Components for RNA PCA and ATAC LSI. |
| `gtf_by` | str | `"gene_name"` | Gene identifier column used in the GTF. |
| `flavor` | str | `"seurat_v3"` | HVG selection flavor. |
| `generate_umap` | bool | `False` | UMAP after preprocessing only. |
| `rna_sample_column` / `atac_sample_column` | str | `"sample"` | Sample columns in each modality. |
| `consistency_threshold` | float | `0.05` | scGLUE consistency loss weight threshold. |
| `treat_sample_as_batch` | bool | `False` | If `True`, each sample becomes its own batch during training (removes per-sample variance). Default `False` preserves per-sample variance for the primary (RMD) embedding. |
| `save_prefix` | str | `"glue"` | Prefix for saved model and integrated files. |
| `batch_key` | str, optional | `None` | Batch column for `scglue.configure_dataset(use_batch=...)`. With `batch_key` set and `treat_sample_as_batch=False`, scGLUE removes the named batch column while preserving per-sample variance, so the primary `X_glue` is suitable as the `Z_rmd` embedding. |
| `sample_key` | str | `"sample"` | Sample identifier column used by training and the union merge. |
| `data_batch_size` | int | `1024` | scGLUE minibatch size. |
| `max_epochs` | int, optional | `None` | Cap on scGLUE training epochs (`None` = scGLUE default). |
| `dataloader_num_workers` / `dataloader_fetches_per_worker` / `array_shuffle_num_workers` / `graph_shuffle_num_workers` | int | `0 / 4 / 0 / 0` | scGLUE dataloader throughput knobs. |
| `run_second_glue_for_sample_removal` | bool | `False` | If `True`, run scGLUE a second time with `treat_sample_as_batch=True` to produce the sample-removed `Z_clust`, then merge both keys into the primary RNA + ATAC h5ads. |
| `second_run_save_prefix` | str | `"glue_no_sample"` | File prefix for the optional second (sample-removal) run. |
| `rna_min_cells` / `rna_min_genes` | int | `500 / 500` | RNA QC filters for the per-modality downstream preprocess. |
| `rna_pct_mito_cutoff` | float | `20.0` | RNA mitochondrial-percent cutoff. |
| `rna_exclude_genes` | list, optional | `None` | RNA genes to drop. |
| `atac_min_cells` | int | `1` | ATAC per-feature min-cells filter. |
| `atac_min_features` / `atac_max_features` | int | `2000 / 15000` | ATAC per-cell feature-count bounds. |
| `atac_min_cells_per_sample` | int | `1` | Minimum cells required per sample for ATAC. |
| `atac_exclude_features` | list, optional | `None` | ATAC features to drop. |
| `atac_doublet_detection` | bool | `True` | Run doublet detection on ATAC. |
| `atac_tfidf_scale_factor` | float | `1e4` | TF-IDF scale factor for ATAC. |
| `atac_log_transform` | bool | `True` | Log-transform ATAC after TF-IDF. |
| `verbose` | bool | `True` | Print progress. |
| `plot_columns` | list, optional | `None` | `.obs` columns to color UMAPs by during the visualization sub-stage. |
| `output_dir` | str | `"./glue_results"` | Writes scGLUE artifacts under `{output_dir}/integration/glue/` and union/per-modality outputs under `{output_dir}/preprocess/`. |

## Returns

The merged union AnnData, with the integrated embedding in `.obsm` (`Z_rmd`, and `Z_clust` when the second run is enabled) and modality metadata populated. Returns `None` only if `run_merge=False` and no existing `preprocess/adata_sample.h5ad` is found.

## Output files

- `{output_dir}/integration/glue/` — preprocessing artifacts, trained model(s), per-modality scGLUE embeddings.
- `{output_dir}/preprocess/adata_sample.h5ad` — canonical embedding-only union object picked up by downstream stages.
- `{output_dir}/preprocess/adata_rna_preprocessed.h5ad`, `{output_dir}/preprocess/adata_atac_preprocessed.h5ad` — per-modality QC'd objects for downstream DGE / RAISIN.

## Usage

```python
from sampledisco.preparation.multi_omics_glue import multiomics_preparation

multiomics_preparation(
    rna_file="/data/test_RNA.h5ad",
    atac_file="/data/test_ATAC.h5ad",
    additional_hvg_file="/data/unique_genes.txt",
    output_dir="/results/multiomics",
    ensembl_release=98,
    species="homo_sapiens",
    consistency_threshold=0.05,
    treat_sample_as_batch=False,
    run_second_glue_for_sample_removal=True,
)
```

!!! tip "Prefer the config-driven wrapper"
    In practice you rarely call `multiomics_preparation` directly. The supported entry point is the YAML-config wrapper, which runs preprocessing → cell embedding → cell typing → sample embedding → downstream analysis for you:

    ```bash
    sampledisco -m complex --config config.yaml
    ```

    scGLUE and `pybedtools` require the `bedtools` binary (`conda install -c bioconda bedtools`). GPU paths activate automatically when the RAPIDS stack is importable.
