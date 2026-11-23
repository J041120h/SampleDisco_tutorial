# `multiomics_preparation`

The scGLUE stage of the multi-omics pipeline, for unpaired scRNA + scATAC data. It runs up to five sub-stages, each switched by a flag:

1. **`run_preprocessing`**: reads both files, merges optional sample metadata, selects RNA HVGs and top ATAC peaks, runs PCA (RNA) and LSI (ATAC), and builds the gene–peak guidance graph from an Ensembl GTF (downloaded for `ensembl_release` / `species`). Needs `bedtools` on `PATH`.
2. **`run_training`**: trains scGLUE. The primary run gives `obsm['X_glue']`, which keeps per-sample variance (the sample-preserved view, `Z_rmd`). With `run_second_glue_for_sample_removal=True`, a second run with the sample as scGLUE batch gives the sample-removed view; both are then written into the primary per-modality files as `obsm['Z_rmd']` and `obsm['Z_comp']`. After training, a consistency score is computed and the log prints `RELIABLE` or `UNRELIABLE` against `consistency_threshold`; nothing else depends on it.
3. **`run_merge`**: builds the **embedding-only union** `preprocess/adata_sample.h5ad`: all RNA and ATAC cells, `obs['modality']`, `obs['original_barcode']`, cell ids suffixed `_RNA` / `_ATAC`, an empty `X`, and the embeddings present in both files (`X_glue`, plus `Z_comp` / `Z_clust` / `Z_rmd` after a second run, and `X_umap` if present). This is done by `build_embedding_union` in `preparation/multi_omics_merge.py`.
4. **`run_preprocess_per_modality`**: per-modality QC and normalization (RNA: the [RNA preprocess](../rna/preprocess_linux.md) filters + normalize/log1p; ATAC: the [ATAC](../atac/preprocess_linux.md) filters + TF-IDF), written to `preprocess/adata_rna_preprocessed.h5ad` and `preprocess/adata_atac_preprocessed.h5ad`. Trajectory DGE and RAISIN read the RNA file.
5. **`run_visualization`**: UMAP plots of the union under `{output_dir}/visualization/`, colored by `plot_columns`.

Without a second scGLUE run the union has only `X_glue`; the wrapper then derives `Z_comp` with [`harmonize_xglue`](harmonize_xglue.md) (which also aliases `X_glue` as `Z_rmd`). The returned union is the input to [`cell_types_multiomics`](cell_types_multiomics_linux.md).

**Import:** `from sampledisco.preparation.multi_omics_glue import multiomics_preparation` (needs `scglue`; importing the module fails without it).

**Source:** `preparation/multi_omics_glue.py:838`

!!! note "Gene-activity step removed"
    Earlier versions imputed per-cell gene activity for ATAC cells by KNN on the joint embedding. That step is gone; the old `run_gene_activity`, `k_neighbors`, `use_rep`, `metric` and `use_gpu` arguments no longer exist.

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
| `rna_file`, `atac_file` | str | — | Input `.h5ad` files (raw counts). |
| `rna_sample_meta_file`, `atac_sample_meta_file` | str or None | `None` | Per-modality sample metadata CSVs. |
| `additional_hvg_file` | str or None | `None` | Plain-text list of genes added to the RNA HVG set. |
| `run_preprocessing` / `run_training` / `run_merge` / `run_preprocess_per_modality` / `run_visualization` | bool | `True` | Sub-stage switches (see above). Without `run_merge`, an existing `preprocess/adata_sample.h5ad` is read and returned. |
| `ensembl_release` | int | `98` | Ensembl release for the GTF. |
| `species` | str | `"homo_sapiens"` | Ensembl species name. |
| `use_highly_variable` | bool | `True` | Train scGLUE on HVGs / top peaks only. |
| `n_top_genes` / `n_top_peaks` | int | `2000` / `50000` | RNA HVGs and ATAC peaks for training. |
| `atac_min_cells_floor` | int | `10` | Floor on the per-peak min-cells filter in the scGLUE preprocessing. |
| `n_pca_comps` / `n_lsi_comps` | int | `50` / `50` | RNA PCA and ATAC LSI components. |
| `gtf_by` | str | `"gene_name"` | GTF attribute matched to the RNA gene names. |
| `flavor` | str | `"seurat_v3"` | HVG flavor. |
| `generate_umap` | bool | `False` | UMAPs during the scGLUE preprocessing. |
| `rna_sample_column` / `atac_sample_column` | str | `"sample"` | Sample columns per modality. |
| `consistency_threshold` | float | `0.05` | Threshold for the post-training consistency check (log only). |
| `treat_sample_as_batch` | bool | `False` | Use the sample as scGLUE batch in the primary run (removes per-sample variance from `X_glue`; not what the pipeline wants for `Z_rmd`). |
| `save_prefix` | str | `"glue"` | Prefix of the primary run's files. |
| `batch_key` | str or None | `None` | scGLUE batch column (`use_batch`) for the primary run; removes that batch while keeping per-sample variance. Wrapper key: `multiomics_batch_col`. |
| `sample_key` | str | `"sample"` | Sample column for training and the union. |
| `data_batch_size` | int | `1024` | scGLUE minibatch size. |
| `max_epochs` | int or None | `None` | Cap on training epochs (`None`: scGLUE's default). |
| `dataloader_num_workers` / `dataloader_fetches_per_worker` / `array_shuffle_num_workers` / `graph_shuffle_num_workers` | int | `0` / `4` / `0` / `0` | scGLUE data-loading workers. The wrapper's defaults are higher (`4` / `8` / `4` / `4`). |
| `run_second_glue_for_sample_removal` | bool | `False` | Train a second scGLUE model with the sample as batch and store its embedding as `Z_comp`. Wrapper key: `multiomics_run_glue_twice_for_sample_removal`. |
| `second_run_save_prefix` | str | `"glue_no_sample"` | Prefix of the second run's files. |
| `rna_min_cells` … `atac_log_transform` | | as in the preprocess functions | QC for the per-modality files (sub-stage 4); same meaning as in the [RNA](../rna/preprocess_linux.md) and [ATAC](../atac/preprocess_linux.md) `preprocess`. |
| `verbose` | bool | `True` | Print progress. |
| `plot_columns` | list or None | `None` | `.obs` columns to color the sub-stage 5 UMAPs by. |
| `output_dir` | str | `"./glue_results"` | Root output directory. |

## Returns

The union `AnnData` from sub-stage 3 (or read from `preprocess/adata_sample.h5ad` when `run_merge=False`); `None` if `run_merge=False` and that file does not exist.

## Output files

- `{output_dir}/integration/glue/`: scGLUE preprocessing files, trained model(s) (`{save_prefix}.dill`), per-modality embeddings (`{save_prefix}-rna-emb.h5ad`, `{save_prefix}-atac-emb.h5ad`), the guidance graph, and the same for `second_run_save_prefix` when enabled.
- `{output_dir}/preprocess/adata_sample.h5ad`: the embedding-only union.
- `{output_dir}/preprocess/adata_rna_preprocessed.h5ad`, `adata_atac_preprocessed.h5ad`: per-modality files for DGE / RAISIN.
- `{output_dir}/visualization/`: sub-stage 5 plots.

## Usage

scGLUE training on the demo takes long enough that this page does not run it. [Multi-omics from a config](../../tutorials/multiomics_config.md) and the [multi-omics tutorial](../../tutorials/multiomics.md) run it, with measured runtimes; the config keys are listed in the [config key reference](../../tutorials/config_reference.md) (`multiomics_run_glue_*`, `multiomics_n_top_genes`, …). To skip scGLUE, start from the pre-integrated demo file as in the [multi-omics example chain](../index.md#running-the-examples).

!!! tip "Prefer the config-driven wrapper"
    With `run_multiomics_pipeline: true`, the wrapper runs this function, then [`harmonize_xglue`](harmonize_xglue.md), [`cell_types_multiomics`](cell_types_multiomics_linux.md), the sample embedding and the downstream steps. scGLUE needs the `bedtools` binary (`conda install -c bioconda bedtools`).
