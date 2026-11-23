# SampleDisco

<div class="hero" markdown>
**Sample-level representation learning for single-cell multi-omics.**

SampleDisco is a config-driven Python pipeline that turns single-cell RNA, ATAC, or unpaired multi-omics data into a unified sample-level embedding — combining multi-resolution cell-type composition with within-cell-type state shifts (RMD displacement) — and runs the full downstream stack of distance analysis, trajectory inference, differential testing, clustering, and visualization in a single call.

[Get started](tutorials/quickstart.md){ .md-button .md-button--primary }
[Browse the API](api/index.md){ .md-button }
</div>

## Workflow overview

![Overview of the SampleDisco framework](resource/figure1.jpg)

<div class="figure-caption"><strong>Figure 1 — Overview of the SampleDisco framework.</strong> From single-cell input and quality control, through cell-level integration and a unified sample-level representation, to downstream sample clustering, trajectory inference, and multilevel differential analysis.</div>

## Four stages

1. **Preprocessing and QC.** Filter cells and features, then build two cell-level embeddings: a sample-removed view (`Z_comp`, Harmony with the sample as a batch) and a sample-preserved view (`Z_rmd`). RNA uses PCA + Harmony, ATAC uses TF-IDF/LSI + Harmony, and unpaired multi-omics uses a joint embedding from scGLUE (or MaxFuse for RNA + flow). Objects written by earlier versions carry [legacy key names](legacy_keys.md), which are still read automatically.
2. **Cell-type assignment.** Leiden clustering, optionally driven to a target number of clusters, or reuse of existing labels.
3. **Sample embedding.** Multi-resolution cell-type composition blocks on `Z_comp` are combined with a reference-relative mean displacement (RMD) block on `Z_rmd`. The blocks are inverse-variance weighted, stacked and PCA-reduced (at most one PC fewer than there are samples). Harmony then corrects the embedding at the sample level, but only when the samples (or sample × modality units) number at least 8 and fall into more than one batch or modality group; otherwise the PCA is kept as is. The result is stored as `adata.uns['X_DR_sample']` ([details](tutorials/sample_embedding.md)).
4. **Downstream analysis.** From that one embedding the pipeline runs sample distances, a trajectory (supervised CCA or unsupervised TSCAN) with trajectory DGE, sample clustering, proportion tests, RAISIN DGE, dimension association and visualization. Phenotype prediction and label transfer are one Python call each.

## Supported inputs

- `scRNA-seq` — a single `.h5ad`; sample/phenotype columns in `.obs` (or an optional metadata CSV).
- `scATAC-seq` — a peak/fragment `.h5ad`; metadata in `.obs` (or an optional CSV).
- `Multi-omics (RNA + ATAC)` — paired or unpaired scRNA / scATAC, joined by scGLUE: trained by the pipeline, or loaded from a pre-integrated `.h5ad`.
- `RNA + protein / flow cytometry` — unpaired modalities matched with MaxFuse ([tutorial](tutorials/multiomics_maxfuse.md)).

!!! tip "Try it on the demo data"
    Ready-made datasets — COVID-19 PBMC scRNA-seq + scATAC-seq, and a 32-sample Sound Life scRNA-seq + flow cytometry subset — and a one-command config are published on Zenodo; see [Demo data](tutorials/demo_data.md) to download them and reproduce every tutorial.

## Quick start

After [installing](installation.md) SampleDisco, fetch the COVID-19 demo data (8 samples) and write the demo config:

```bash
mkdir -p sampledisco_quickstart/data
cd sampledisco_quickstart
curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
sampledisco --init-config config.yaml
```

Then run the whole pipeline, RNA and ATAC, from the shell or from Python:

=== "CLI"

    ```bash
    sampledisco -m complex --config config.yaml
    ```

=== "Python"

    ```python
    import yaml
    from sampledisco.wrapper.wrapper import wrapper

    with open("config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    results = wrapper(**config)
    ```

On 8 CPU cores either form takes 3 to 5.5 minutes, depending on the machine, and peaks at about 8.6 GB of RAM. Results land in `./sampledisco_demo_output/`.

The [Quick start](tutorials/quickstart.md) page walks through this run step by step: an RNA-only run first, what each stage does, and how to load the results back. [Understanding the outputs](tutorials/outputs.md) maps every output file to the stage that wrote it. From there:

- **Run from a config** — [Configuration guide](tutorials/configuration.md), [every config key](tutorials/config_reference.md), [GPU acceleration](tutorials/gpu.md), [using your own data](tutorials/own_data.md)
- **Step by step in Python** — [RNA](tutorials/rna.md), [ATAC](tutorials/atac.md), [the sample embedding in depth](tutorials/sample_embedding.md)
- **Multi-omics** — [RNA + ATAC with scGLUE](tutorials/multiomics.md), [multi-omics from a config](tutorials/multiomics_config.md), [RNA + flow cytometry with MaxFuse](tutorials/multiomics_maxfuse.md)
- **Downstream analysis** — [overview](tutorials/downstream/index.md), and a larger [32-sample walkthrough](tutorials/soundlife_rna.md) with phenotype prediction

## Citation

A manuscript describing SampleDisco is in preparation. A citation block will be added here once the preprint is posted.
