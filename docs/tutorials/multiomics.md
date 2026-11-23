# RNA + ATAC (scGLUE)

This page builds one sample embedding from scRNA-seq and scATAC-seq measured on **different cells**, here even on different donors. [scGLUE](https://www.nature.com/articles/s41587-022-01284-4) puts both modalities into one joint cell embedding. SampleDisco then types the cells jointly and gives every sample one point per modality (`<sample>_RNA`, `<sample>_ATAC`). After that, every downstream analysis works as it does for RNA or ATAC alone.

SampleDisco reads two views of the joint embedding:

- **`Z_rmd`** is scGLUE's own `X_glue`, with the differences between samples kept. It feeds the RMD displacement block.
- **`Z_comp`** is `X_glue` after one Harmony pass with `sample` as the batch key, so the sample effect is removed. It feeds the cell typing and the composition blocks.

Pick a route for step 1. Steps 2–4 are the same for both. Run everything from one empty folder.

## 1. Get the joint embedding

=== "Quick route: pre-integrated file"

    Start from the demo object that already holds the scGLUE embedding: 29,989 RNA cells from 8 donors and 29,021 ATAC cells from 8 *other* donors, half healthy and half severe COVID-19 (`obs['sev.level']` 1 or 4). This route needs only `pip install sampledisco`. Steps 1–4 take about 2 minutes on 8 CPU cores and about 3 GB of RAM.

    ```bash
    mkdir -p data
    curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
    ```

    ```python
    import anndata as ad
    from sampledisco.utils.embedding_keys import migrate_obsm_keys

    OUT = "sampledisco_demo_output/multiomics_api"

    adata = ad.read_h5ad("data/test_multiomics_integrated.h5ad")
    migrate_obsm_keys(adata)   # Z_clust -> Z_comp, Z_cmd -> Z_rmd (in memory; the file is not changed)

    # The file ships the cell types and sample embedding of an older run. Steps 2 and 3 recompute them.
    adata.obs = adata.obs.drop(columns=["cell_type"], errors="ignore")
    for key in ("X_DR_sample", "sample_embedding_params", "cell_type_colors"):
        adata.uns.pop(key, None)
    ```

    The file stores the embedding under older key names. `migrate_obsm_keys` copies them onto `Z_comp` and `Z_rmd` (see [legacy key names](../api/utils/embedding_keys.md)). The file has no expression data: its `.X` is empty.

=== "Full route: train scGLUE"

    Train scGLUE yourself from the raw counts. This needs scGLUE, PyTorch, pyensembl and the `bedtools` binary; the all-in-one environment from [Installation](../installation.md) has all of them. The whole route takes about 25–45 minutes on a GPU node, depending on the GPU, and about 40 minutes CPU-only. Ask for 20 GB of RAM and 2 GB of disk.

    ```bash
    mkdir -p data
    curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
    curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
    ```

    ```python
    import anndata as ad
    from sampledisco.preparation.multi_omics_glue import multiomics_preparation
    from sampledisco.preparation.multi_omics_batch_correction import harmonize_xglue

    OUT = "sampledisco_demo_output/multiomics_glue"

    multiomics_preparation(
        rna_file="data/test_RNA.h5ad",
        atac_file="data/test_ATAC.h5ad",
        output_dir=OUT,
        ensembl_release=98,
        species="homo_sapiens",
        n_top_genes=2000,
        n_top_peaks=50000,
        n_pca_comps=50,
        n_lsi_comps=50,
        rna_sample_column="sample",
        atac_sample_column="sample",
        data_batch_size=128,
        max_epochs=50,
        plot_columns=["modality", "sev.level"],
    )

    adata = ad.read_h5ad(f"{OUT}/preprocess/adata_sample.h5ad")
    adata = harmonize_xglue(adata, sample_col="sample", batch_col=None, max_iter=30, random_state=42)
    ```

    - `multiomics_preparation` builds the gene–peak guidance graph and trains scGLUE, on a visible CUDA GPU if there is one. It then writes the union of all cells, with `obsm['X_glue']`, to `preprocess/adata_sample.h5ad`. It also writes QC-filtered counts for each modality (`preprocess/adata_rna_preprocessed.h5ad`, `adata_atac_preprocessed.h5ad`), which DGE and RAISIN use.
    - `harmonize_xglue` adds `Z_comp` and copies `X_glue` to `Z_rmd`.
    - On first use, pyensembl downloads and indexes the Ensembl 98 annotation (about 150 MB, a few minutes) into `~/.cache/pyensembl/`. If your compute nodes have no internet, run `pyensembl install --release 98 --species homo_sapiens` once on a machine that has it and shares the same home folder.

## 2. Joint cell typing

`cell_types_multiomics` runs Leiden clustering on the RNA cells in `Z_comp` (`obs['modality'] == "RNA"`) and transfers the labels to the ATAC cells through shared nearest neighbors. It drops clusters made almost entirely of one modality and returns a **new object**: assign the result. `k_neighbors=3` is the value the config pipeline uses.

```python
from sampledisco.preparation.multi_omics_cell_type_cpu import cell_types_multiomics

adata = cell_types_multiomics(
    adata,
    use_rep="Z_comp",
    k_neighbors=3,
    save=True,
    output_dir=OUT,
    defined_output_path=f"{OUT}/preprocess/adata_sample.h5ad",
)
print(adata.n_obs, "cells kept")
```

On the quick route this finds 15 or 16 clusters, depending on the library versions, drops none and keeps all 59,010 cells.

![UMAP of the RNA and ATAC cells, colored by joint cell type](../resource/multiomics/api_umap_split_by_modality.png)
<div class="figure-caption"><code>preprocess/visualization/umap_split_by_modality.png</code>: RNA cells (left, Leiden clusters) and ATAC cells (right, transferred labels) in the same UMAP of <code>Z_comp</code>.</div>

## 3. Sample embedding

`compute_sample_embedding` is the same function as for RNA and ATAC. With `modality_col`, each sample × modality pair becomes one unit. The composition blocks use `Z_comp`, the RMD block uses `Z_rmd`, and a final Harmony pass removes the modality offset between the units. Its other options are explained in [Sample embedding in depth](sample_embedding.md).

```python
from sampledisco.sample_embedding import compute_sample_embedding
from sampledisco.sample_embedding.blocks import save_embedding_to_h5ad

adata = compute_sample_embedding(
    adata,
    output_dir=OUT,
    sample_col="sample",
    celltype_col="cell_type",
    comp_emb_key="Z_comp",      # sample removed
    rmd_emb_key="Z_rmd",        # sample differences kept (= X_glue)
    modality_col="modality",
)
save_embedding_to_h5ad(f"{OUT}/preprocess/adata_sample.h5ad", adata)   # keep the embedding with the cells
print(adata.uns["X_DR_sample"].shape)   # (16, 10): 16 units x 10 PCs
```

## 4. Plot the sample embedding

The embedding has one row per unit, while `adata.obs` has one row per cell. `build_sample_adata` builds the matching unit-level object, which the plots and every downstream module accept.

```python
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.visualization.multi_omics_visualization import visualize_multimodal_embedding

sample_adata = build_sample_adata(adata, sample_col="sample", modality_col="modality")
visualize_multimodal_embedding(adata=sample_adata, visualization_grouping_column=["sev.level", "modality"],
                               output_dir=OUT)   # the plots go to OUT/visualization/
```

![Sample embedding colored by sev.level](../resource/multiomics/api_sample_embedding_by_sev_level.png)
<div class="figure-caption">PC1 against PC2 of the 16 units: severe COVID-19 (<code>sev.level</code> 4, yellow) against healthy (1, purple). RNA and ATAC units mix, and the phenotypes mostly split along PC1: 6 of the 8 severe units have PC1 &gt; 0, and 7 of the 8 healthy units PC1 &lt; 0.</div>

With 4 + 4 donors per modality this is a picture, not a test. The full route gives the same split, although its clusters and coordinates differ, because scGLUE training is not bit-reproducible.

## Output

All under `OUT` (`sampledisco_demo_output/multiomics_api/` or `.../multiomics_glue/`):

- `preprocess/adata_sample.h5ad`: the cells with `cell_type`, `Z_comp`, `Z_rmd` and `uns['X_DR_sample']`. Reload it and call `build_sample_adata` to continue in a new session.
- `preprocess/cell_type.csv` and `preprocess/visualization/`: the cell types and the cell-type UMAPs.
- `sample_embedding/sample_embedding.csv`: 16 units × 10 PCs.
- `visualization/all_samples_*_by_{sev.level,modality}.png`: the sample-embedding plots.
- Full route only: `integration/glue/` (the scGLUE model and its inputs) and `preprocess/adata_{rna,atac}_preprocessed.h5ad`. Before DGE or RAISIN in Python, copy the new cell types onto these two files:

    ```python
    from sampledisco.preparation.multi_omics_merge import propagate_cell_type
    propagate_cell_type(union_path=f"{OUT}/preprocess/adata_sample.h5ad", celltype_col="cell_type",
                        per_modality_paths=[f"{OUT}/preprocess/adata_{m}_preprocessed.h5ad" for m in ("rna", "atac")])
    ```

## From a config

`sampledisco --config` runs the same steps, then the downstream analyses. Write the template with `sampledisco --init-config config.yaml`, set these keys, and run it. The [configuration guide](configuration.md) and the [config key reference](config_reference.md) explain every `multiomics_*` key.

```yaml
run_rna_pipeline: false
run_atac_pipeline: false
run_multiomics_pipeline: true
multiomics_integration: false   # true: train scGLUE from multiomics_rna_file / multiomics_atac_file
multiomics_integrated_h5ad_path: "./data/test_multiomics_integrated.h5ad"
```

```bash
sampledisco --config config.yaml
```

- The results go to `sampledisco_demo_output/multiomics/`; the input file is never changed. Set `multiomics_visualization_grouping_column: ["sev.level"]` to color the plots.
- Trajectory DGE and RAISIN (`multiomics_trajectory_dge`, `multiomics_cluster_dge`) need the per-modality counts that only the full route writes. With `multiomics_integration: false` they are skipped with a warning.
- `multiomics_use_gpu` takes effect only together with `use_gpu: true`. scGLUE training uses a visible GPU either way. See [GPU acceleration](gpu.md).

Continue with [Downstream analysis](downstream/index.md): every analysis there reads `uns['X_DR_sample']` the same way.
