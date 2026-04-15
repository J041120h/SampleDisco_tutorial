# SampleDisc

<div class="hero" markdown>
**Sample-level representation learning for single-cell multi-omics.**

SampleDisc is a config-driven Python pipeline that turns single-cell RNA, ATAC, or unpaired multi-omics data into dual sample-level embeddings — one capturing expression signal, one capturing cell-type composition — and runs the full downstream stack of distance analysis, trajectory inference, differential testing, clustering, and visualization in a single call.

[Get started](tutorials/configuration.md){ .md-button .md-button--primary }
[Browse the API](api/index.md){ .md-button }
</div>

<div class="grid cards" markdown>

-   **Three pipelines, one wrapper**

    ---

    Point `wrapper(...)` at a YAML config and run scRNA, scATAC, or integrated multi-omics end-to-end without glue code.

-   **Dual sample representation**

    ---

    Every run produces `X_DR_expression` (pseudobulk signal) and `X_DR_proportion` (cell-type composition) embeddings, giving you two complementary views of each sample.

-   **Batch-aware, trajectory-aware**

    ---

    Harmony integration on both cell and sample levels, optional CCA-driven resolution search against a phenotype, and supervised or unsupervised pseudotime.

-   **Downstream out of the box**

    ---

    Shared downstream module runs sample distance (six metrics), CCA/TSCAN, GAM-based trajectory DGE, K-means sample clustering, proportion tests, and RAISIN cluster DGE.

</div>

## Workflow overview

![SampleDisc workflow](resource/overview.jpg)

<div class="figure-caption">From cell-level input to sample-level embeddings and downstream trajectory, clustering, and differential analyses.</div>

## Four stages

1. **Preprocessing and QC.** Filter cells and features, correct batch, build cell-level embeddings (PCA/Harmony for RNA, TF-IDF/LSI/Harmony for ATAC, GLUE for unpaired multi-omics).
2. **Cell-type assignment.** Leiden clustering with optional target cluster count, or reuse of existing labels.
3. **Sample embedding.** Pseudobulk per cell type + PCA → expression embedding. Cell-type proportions + PCA (optionally Harmony) → proportion embedding.
4. **Downstream analysis.** Distance, trajectory, differential genes, clustering, visualization — all running off the same pseudobulk object.

## Supported inputs

- `scRNA-seq` — a single `.h5ad` + sample metadata CSV.
- `scATAC-seq` — a peak/fragment `.h5ad` + sample metadata CSV.
- `Unpaired multi-omics (RNA + ATAC)` — two `.h5ad` files integrated via GLUE.
- `Paired multi-omics` — same two-file entry point; GLUE still anchors the joint embedding.

## Quick start

=== "CLI"

    ```bash
    python /users/hjiang/GenoDistance/code/SampleDisc.py -m complex \
      --config /users/hjiang/GenoDistance/code/config/config_covid_rna.yaml
    ```

=== "Python"

    ```python
    import yaml
    from wrapper.wrapper import wrapper

    with open("config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    wrapper(**config)
    ```

Continue to the [Configuration guide](tutorials/configuration.md) for a walkthrough of the YAML, or jump into one of the pipeline tutorials:

- [RNA pipeline](tutorials/rna.md)
- [ATAC pipeline](tutorials/atac.md)
- [Multi-omics pipeline](tutorials/multiomics.md)
- [Downstream analysis](tutorials/downstream.md)

## Citation

A manuscript describing SampleDisc is in preparation. A citation block will be added here once the preprint is posted.
