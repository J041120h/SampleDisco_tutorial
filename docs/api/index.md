# API reference

The API is organized by the role each function plays in a SampleDisco run. Functions shared between modalities (cell typing, sample embedding, and the whole downstream module) are documented once and cross-linked.

!!! tip "The config-driven wrapper is the supported entry point"
    Most users should not call these functions directly. The primary interface is the YAML-config wrapper — `sampledisco -m complex --config <config.yaml>` (equivalently `python -m sampledisco.cli -m complex --config <config.yaml>`), which dispatches `sampledisco.wrapper.wrapper.wrapper(**config)`. It orchestrates preprocess → cell embedding → cell typing → sample embedding → downstream for each enabled modality. The functions below are the building blocks it calls.

!!! note "Public import paths"
    The import root is `sampledisco` (NOT `genodistance`). Subpackage `__init__` files are empty, so import each public function from its concrete module — e.g. `from sampledisco.preparation.rna_preprocess_cpu import preprocess`, `from sampledisco.sample_embedding import compute_sample_embedding`. The two top-level re-exports are `sampledisco.wrapper` and `sampledisco.compute_sample_embedding`.

## Categories

| Section | When you reach for it | Key functions |
| --- | --- | --- |
| [RNA](rna.md) | scRNA-seq preprocessing and clustering. | `preprocess`, `cell_types` |
| [ATAC](atac.md) | scATAC-seq preprocessing and cell typing. | `preprocess` (ATAC build), `cell_types_atac` |
| [Multi-omics](multiomics.md) | GLUE integration, joint cell typing, multi-omics sample embedding, embedding merge. | `multiomics_preparation`, `cell_types_multiomics`, `compute_sample_embedding` |
| [Shared](shared/calculate_sample_embedding.md) | The core sample embedding — reused by all modalities. | `compute_sample_embedding` |
| [Downstream](downstream.md) | Everything that runs off the sample embedding: distance, trajectory, DGE, clustering, visualization. | `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `raisinfit`, `run_pairwise_tests`, `run_dimension_association_analysis`, `visualization`, `visualize_multimodal_embedding` |

!!! warning "Removed in the current API"
    Three functions earlier versions documented no longer exist and have no drop-in replacement: `find_optimal_cell_resolution_linux` and `find_optimal_cell_resolution_multiomics_linux` (the CCA-based resolution sweep — superseded by alpha/block-weight autotune, `sampledisco.parameter_selection.autotune.run_autotune`, wired via the `*_autotune_enable` wrapper flags), and `integrate_preprocess` (its per-modality QC + union build is now folded into `multiomics_preparation`). The two-key sample embedding (`X_DR_expression` / `X_DR_proportion`) has been replaced by a single key, `adata.uns['X_DR_sample']`.

If you are unsure which function to call next, the tutorials walk through the end-to-end order: [RNA](../tutorials/rna.md), [ATAC](../tutorials/atac.md), [Multi-omics](../tutorials/multiomics.md), [Downstream](../tutorials/downstream/index.md).
