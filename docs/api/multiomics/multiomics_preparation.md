# `multiomics_preparation`

Full GLUE integration pipeline for unpaired (or paired) scRNA + scATAC data. Internally runs four sub-stages, each individually toggleable via `run_preprocessing`, `run_training`, `run_gene_activity`, `run_visualization`:

1. **GLUE preprocessing** — reads both modalities, merges sample metadata, selects HVGs on RNA, runs LSI on ATAC, and builds a gene-region guidance graph from an Ensembl GTF.
2. **Adversarial GLUE training** — fits the joint embedding with a consistency-regularized adversarial loss; optionally treats each sample as its own batch.
3. **Gene activity computation** — for each ATAC cell, pulls nearest RNA neighbors on the joint embedding and imputes gene activity scores.
4. **Visualization (optional)** — UMAPs on the joint space, colored by modality or user-specified columns.

The returned object is a merged AnnData with `X_glue` in `.obsm` ready for [`integrate_preprocess`](integrate_preprocess.md) and [`cell_types_multiomics_linux`](cell_types_multiomics_linux.md).

**Source:** `preparation/multi_omics_glue.py:952`

## Signature

```python
def multiomics_preparation(
    rna_file: str,
    atac_file: str,
    rna_sample_meta_file: Optional[str] = None,
    atac_sample_meta_file: Optional[str] = None,
    additional_hvg_file: Optional[str] = None,
    # Process control flags
    run_preprocessing: bool = True,
    run_training: bool = True,
    run_gene_activity: bool = True,
    run_visualization: bool = True,
    # Preprocessing parameters
    ensembl_release: int = 98,
    species: str = "homo_sapiens",
    use_highly_variable: bool = True,
    n_top_genes: int = 2000,
    n_pca_comps: int = 50,
    n_lsi_comps: int = 50,
    gtf_by: str = "gene_name",
    flavor: str = "seurat_v3",
    generate_umap: bool = False,
    rna_sample_column: str = "sample",
    atac_sample_column: str = "sample",
    # Training parameters
    consistency_threshold: float = 0.05,
    treat_sample_as_batch: bool = True,
    save_prefix: str = "glue",
    # Gene activity
    k_neighbors: int = 1,
    use_rep: str = "X_glue",
    metric: str = "cosine",
    use_gpu: bool = True,
    verbose: bool = True,
    # Visualization
    plot_columns: Optional[List[str]] = None,
    # Output
    output_dir: str = "./glue_results",
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `rna_file`, `atac_file` | str | — | Input `.h5ad` files. |
| `rna_sample_meta_file`, `atac_sample_meta_file` | str, optional | `None` | Per-modality sample metadata CSVs. |
| `additional_hvg_file` | str, optional | `None` | Plain-text list of extra HVGs forced into the RNA HVG set. |
| `run_preprocessing` / `_training` / `_gene_activity` / `_visualization` | bool | `True` | Toggle individual sub-stages. |
| `ensembl_release` | int | `98` | Ensembl release used for GTF retrieval. |
| `species` | str | `"homo_sapiens"` | Species key matching Ensembl. |
| `use_highly_variable` | bool | `True` | Restrict GLUE training to HVGs. |
| `n_top_genes` | int | `2000` | HVG count for RNA. |
| `n_pca_comps` / `n_lsi_comps` | int | `50 / 50` | Components for RNA PCA and ATAC LSI. |
| `gtf_by` | str | `"gene_name"` | Gene identifier column used in the GTF. |
| `flavor` | str | `"seurat_v3"` | HVG selection flavor. |
| `generate_umap` | bool | `False` | UMAP after preprocessing only. |
| `rna_sample_column` / `atac_sample_column` | str | `"sample"` | Sample columns in each modality. |
| `consistency_threshold` | float | `0.05` | GLUE consistency loss weight threshold. |
| `treat_sample_as_batch` | bool | `True` | Each sample becomes its own batch dimension during training. |
| `save_prefix` | str | `"glue"` | Prefix for saved model and integrated files. |
| `k_neighbors` | int | `1` | Nearest RNA neighbors per ATAC cell for gene activity. |
| `use_rep` | str | `"X_glue"` | Representation to use for nearest-neighbor search. |
| `metric` | str | `"cosine"` | Distance metric for the neighbor search. |
| `use_gpu` | bool | `True` | Run the GPU-accelerated gene-activity step when available. |
| `verbose` | bool | `True` | Print progress. |
| `plot_columns` | list, optional | `None` | `.obs` columns to color UMAPs by during the visualization sub-stage. |
| `output_dir` | str | `"./glue_results"` | Writes under `{output_dir}/integration/glue/`. |

## Returns

The merged integrated AnnData, with `.obsm["X_glue"]` and modality/gene-activity metadata populated.

## Output files

- `{output_dir}/integration/glue/` — preprocessing artifacts, trained model, integrated AnnData.
- `{output_dir}/preprocess/adata_sample.h5ad` — canonical integrated object picked up by downstream stages.

## Usage

```python
from genodistance.preparation import multiomics_preparation

multiomics_preparation(
    rna_file="/data/test_RNA.h5ad",
    atac_file="/data/test_ATAC.h5ad",
    additional_hvg_file="/data/unique_genes.txt",
    output_dir="/results/multiomics",
    ensembl_release=98,
    species="homo_sapiens",
    consistency_threshold=0.05,
    treat_sample_as_batch=True,
    k_neighbors=1,
    use_rep="X_glue",
)
```
