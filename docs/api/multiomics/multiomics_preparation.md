# `multiomics_preparation`

**Module**: `code/preparation/multi_omics_glue.py`

```python
multiomics_preparation(
    rna_file,
    atac_file,
    rna_sample_meta_file=None,
    atac_sample_meta_file=None,
    additional_hvg_file=None,
    run_preprocessing=True,
    run_training=True,
    run_gene_activity=True,
    run_visualization=True,
    ensembl_release=98,
    species="homo_sapiens",
    use_highly_variable=True,
    n_top_genes=2000,
    n_pca_comps=50,
    n_lsi_comps=50,
    gtf_by="gene_name",
    flavor="seurat_v3",
    generate_umap=False,
    rna_sample_column="sample",
    atac_sample_column="sample",
    consistency_threshold=0.05,
    treat_sample_as_batch=True,
    save_prefix="glue",
    k_neighbors=1,
    use_rep="X_glue",
    metric="cosine",
    use_gpu=True,
    verbose=True,
    plot_columns=None,
    output_dir="./glue_results",
)
```

Runs GLUE preparation/training/gene activity/visualization stages.

## Parameter guide

| Parameter | Default | Controls |
| --- | --- | --- |
| `rna_file`, `atac_file` | required | RNA/ATAC input files. |
| `rna_sample_meta_file`, `atac_sample_meta_file` | `None` | Optional modality-specific sample metadata CSVs. |
| `additional_hvg_file` | `None` | Extra HVG list for preprocessing. |
| `run_preprocessing` | `True` | Run GLUE preprocessing stage. |
| `run_training` | `True` | Train GLUE model. |
| `run_gene_activity` | `True` | Compute gene activity and merged integrated object. |
| `run_visualization` | `True` | Generate integration visualization plots. |
| `ensembl_release` | `98` | Ensembl release used by annotation step. |
| `species` | `"homo_sapiens"` | Species label for gene annotation mapping. |
| `use_highly_variable` | `True` | Restrict to highly variable genes/features. |
| `n_top_genes` | `2000` | HVG target count for RNA preprocessing. |
| `n_pca_comps` | `50` | PCA dimensions for RNA branch. |
| `n_lsi_comps` | `50` | LSI dimensions for ATAC branch. |
| `gtf_by` | `"gene_name"` | Feature key used in GTF mapping. |
| `flavor` | `"seurat_v3"` | HVG selection flavor. |
| `generate_umap` | `False` | Generate UMAP during preprocess. |
| `rna_sample_column`, `atac_sample_column` | `"sample"` | Sample ID columns per modality. |
| `consistency_threshold` | `0.05` | GLUE consistency filtering threshold. |
| `treat_sample_as_batch` | `True` | Treat sample ID as batch variable during training. |
| `save_prefix` | `"glue"` | Output prefix for saved model artifacts. |
| `k_neighbors` | `1` | KNN size for gene activity transfer. |
| `use_rep` | `"X_glue"` | Embedding used for KNN transfer. |
| `metric` | `"cosine"` | Distance metric for KNN. |
| `use_gpu` | `True` | Use GPU acceleration for transfer steps. |
| `verbose` | `True` | Print detailed progress logs. |
| `plot_columns` | `None` | Metadata columns to color integration visualizations. |
| `output_dir` | `"./glue_results"` | Output root. |

