# API reference

The API is organized by the role each function plays in a SampleDisc run. Functions shared between modalities (cell typing, sample embedding, resolution search, and the whole downstream module) are documented once and cross-linked.

## Categories

| Section | When you reach for it | Key functions |
| --- | --- | --- |
| [RNA](rna.md) | scRNA-seq preprocessing and clustering. | `preprocess_linux`, `cell_types_linux` |
| [ATAC](atac.md) | scATAC-seq preprocessing (cell typing is shared with RNA). | `preprocess_linux` (ATAC build) |
| [Multi-omics](multiomics.md) | GLUE integration, integrated QC, joint cell typing, multi-omics sample embedding, resolution search, embedding merge. | `multiomics_preparation`, `integrate_preprocess`, `cell_types_multiomics_linux`, `calculate_multiomics_sample_embedding`, `find_optimal_cell_resolution_multiomics_linux`, `replace_optimal_dimension_reduction` |
| [Shared](shared/calculate_sample_embedding.md) | Pseudobulk sample embedding and CCA-guided resolution search — reused by all modalities. | `calculate_sample_embedding`, `find_optimal_cell_resolution_linux` |
| [Downstream](downstream.md) | Everything that runs off a pseudobulk: distance, trajectory, DGE, clustering, visualization. | `sample_distance`, `CCA_Call`, `cca_pvalue_test`, `TSCAN`, `run_trajectory_gam_differential_gene_analysis`, `cluster`, `proportion_test`, `raisinfit`, `run_pairwise_tests`, `visualization`, `visualize_multimodal_embedding` |

## Page conventions

Every function page follows the same layout:

1. **Description** — 3–6 sentences on the algorithm, its inputs and outputs, and when you would call it.
2. **Signature** — copied directly from source, with full type hints and defaults.
3. **Parameters table** — name, type, default, meaning.
4. **Returns** — what the function hands back.
5. **Output files** — where artifacts land on disk.
6. **Usage** — a self-contained minimal call.

If you are unsure which function to call next, the tutorials walk through the end-to-end order: [RNA](../tutorials/rna.md), [ATAC](../tutorials/atac.md), [Multi-omics](../tutorials/multiomics.md), [Downstream](../tutorials/downstream.md).
