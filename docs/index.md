# SampleDisco

<div class="hero" markdown>
**Sample-level representation learning for single-cell multi-omics.**

SampleDisco is a config-driven Python pipeline that turns single-cell RNA, ATAC, or unpaired multi-omics data into a unified sample-level embedding — combining multi-resolution cell-type composition with within-cell-type state shifts (RMD displacement) — and runs the full downstream stack of distance analysis, trajectory inference, differential testing, clustering, and visualization in a single call.

[Get started](tutorials/demo_data.md){ .md-button .md-button--primary }
[Browse the API](api/index.md){ .md-button }
</div>

## Workflow overview

![SampleDisco workflow](resource/overview.jpg)

<div class="figure-caption">From cell-level input to a sample-level embedding and downstream trajectory, clustering, and differential analyses.</div>

## Four stages

1. **Preprocessing and QC.** Filter cells and features, correct batch, build cell-level embeddings — a sample-removed view (`Z_clust`) and a sample-preserved view (`Z_rmd`) — via PCA/Harmony for RNA, TF-IDF/LSI/Harmony for ATAC, and scGLUE for unpaired multi-omics.
2. **Cell-type assignment.** Leiden clustering with optional target cluster count, or reuse of existing labels.
3. **Sample embedding.** Multi-resolution cell-type composition blocks on `Z_clust` are combined with a reference-relative mean displacement (RMD) block on `Z_rmd`; the blocks are inverse-variance weighted, stacked, PCA-reduced, and Harmony-corrected at the sample level into a single embedding stored as `adata.uns['X_DR_sample']`.
4. **Downstream analysis.** Distance, trajectory, differential genes, clustering, visualization — all running off the same sample embedding.

## Supported inputs

- `scRNA-seq` — a single `.h5ad`; sample/phenotype columns in `.obs` (or an optional metadata CSV).
- `scATAC-seq` — a peak/fragment `.h5ad`; metadata in `.obs` (or an optional CSV).
- `Unpaired multi-omics (RNA + ATAC)` — two `.h5ad` files integrated via GLUE.
- `Paired multi-omics` — same two-file entry point; GLUE still anchors the joint embedding.

!!! tip "Try it on the demo data"
    Two ready-made datasets (scRNA-seq + scATAC-seq, COVID-19 PBMC) and a one-command config are published on Zenodo — see [Demo data](tutorials/demo_data.md) to download and reproduce every tutorial.

## Quick start

**Run the full demo in three steps** (after `pip install sampledisco` — see [Installation](installation.md)):

```bash
# 1. generate a ready-to-run config (pre-wired to ./data/test_*.h5ad)
sampledisco --init-config config_demo.yaml

# 2. download the demo data into ./data  (full commands + checksums: Demo data page)
mkdir -p data
curl -L -o data/test_RNA.h5ad  "https://zenodo.org/records/21019419/files/test_RNA.h5ad?download=1"
curl -L -o data/test_ATAC.h5ad "https://zenodo.org/records/21019419/files/test_ATAC.h5ad?download=1"

# 3. run RNA + ATAC end to end
sampledisco -m complex --config config_demo.yaml
```

The general form (any config):

=== "CLI"

    ```bash
    sampledisco -m complex --config config.yaml
    # equivalent: python -m sampledisco.cli -m complex --config config.yaml
    ```

=== "Python"

    ```python
    import yaml
    from sampledisco import wrapper

    with open("config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    wrapper(**config)
    ```

Start with the [demo data](tutorials/demo_data.md), then work through the pipeline tutorials. To drive the whole pipeline from a single YAML instead, see the [Configuration guide](tutorials/configuration.md).

- [Demo data](tutorials/demo_data.md)
- [RNA pipeline](tutorials/rna.md)
- [ATAC pipeline](tutorials/atac.md)
- [Multi-omics pipeline](tutorials/multiomics.md)
- [Downstream analysis](tutorials/downstream/index.md)
- [Configuration guide](tutorials/configuration.md)

## Citation

A manuscript describing SampleDisco is in preparation. A citation block will be added here once the preprint is posted.
