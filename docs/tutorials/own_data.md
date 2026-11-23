# Using your own data

This page covers what changes when you move from the demo to your own cohort: what the input files must contain, how to bring in sample metadata and existing cell-type labels, where small cohorts hit hard limits, and what the common errors mean. Start from a fresh template (`sampledisco --init-config my_config.yaml`) and the list of keys to change in the [Configuration guide](configuration.md).

## What the input must contain

Each modality is one cell-level `.h5ad` file (cells × features) holding **all samples** of the cohort.

| Part of the file | Requirement |
|---|---|
| `.X` | **Raw counts**: UMI counts for RNA, fragment or peak counts for ATAC. SampleDisco normalizes and selects features itself. Do not pass log-normalized data. |
| `.var_names` | Gene symbols for RNA (mitochondrial genes are found by the `MT-` / `mt-` prefix); peaks for ATAC. |
| `.obs[<sample column>]` | The sample (donor, specimen) of each cell. Name it **`sample`** if you can: several functions you may call directly fall back to `obs['sample']`. |
| `.obs` sample-level columns | The phenotypes to analyse (disease group, age, …) and any batch column, constant within a sample. They can also come from a CSV (below). |
| `.obs` cell-type labels | Optional. Without them, SampleDisco clusters the cells. |

## Check a file with `summarize_h5ad`

`summarize_h5ad` prints the shape, the `.obs` columns with their types and top values, and a check of `.X`:

```python
from sampledisco.utils.inspector import summarize_h5ad

summarize_h5ad("data/test_RNA.h5ad", n_examples=3, print_unique_column="sample")
```

```text
📦 Basic Info
  - Shape (cells × genes): 29989 × 14881
  - obs columns: ['celltype', 'sample', 'old', 'large', 'batch', 'sev.level']
…
🔎 Unique values for column 'sample' (found in adata.obs):
  - n_unique: 8
…
🔎 Inspecting .X matrix:
  - value type: float — values appear integer-like (no fractional parts)
  - max value: 30539
…
```

The output is shortened at `…`. What to look for:

- **`integer-like`** values and a large maximum mean `.X` holds raw counts. Fractional values or a maximum around 5–10 mean log-normalized data.
- **`print_unique_column`** lists every sample ID: check the count, and look for typos that split one sample in two.
- Each sample-level column should have a handful of unique values. On the demo, `batch` and `sev.level` have 2 each; `celltype` is a per-cell label.

The call loads the whole file into memory (5 s and 1.3 GB on the demo).

## Sample and cell metadata files

When the phenotypes are not in `.obs`, put them in a **sample-metadata CSV** and set `rna_sample_meta_path` (`atac_sample_meta_path`; for multi-omics `multiomics_rna_sample_meta_file` and `multiomics_atac_sample_meta_file`):

- One row per sample, a header row, and a column named like `*_sample_col` (normally `sample`) holding the same IDs as `.obs`.
- Comma, tab or semicolon separators are detected, and an Excel byte-order mark is handled.
- Every cell gets its sample's values. CSV columns **replace** same-named `.obs` columns. Cells whose sample is not in the CSV get missing values.

This builds such a file from the demo and merges it back with the function the pipeline uses:

```python
import anndata as ad
from sampledisco.utils.merge_sample_meta import merge_sample_metadata

adata = ad.read_h5ad("data/test_RNA.h5ad")
meta = adata.obs.groupby("sample", observed=True)[["batch", "sev.level"]].first()
meta.to_csv("sample_meta.csv")

adata.obs = adata.obs.drop(columns=["batch", "sev.level"])   # pretend they were missing
adata = merge_sample_metadata(adata, "sample_meta.csv", sample_column="sample")
```

```text
   📄 Reading CSV metadata file: sample_meta.csv
   ✅ Added 2 new metadata columns
   🔗 Matched 29989/29989 entries (100.0%)
```

In a config run, check the `Matched … entries` line in the log: below 100% means some cells' sample IDs are not in the CSV.

A **cell-metadata CSV** (`rna_cell_meta_path`, `atac_cell_meta_path`) adds per-cell columns. It needs a column named **`barcode`** matching `obs_names`. It is read before the sample CSV, so it can supply the sample column itself.

## Example: adapting the demo config

The [Sound Life subset](demo_data.md) (32 donors, `soundlife_RNA.h5ad`) has `sample`, `batch`, `age_group`, `age_at_draw`, `sex`, `cmv` and cell-type labels in `cell_type`, but no `sev.level`. In a fresh `my_config.yaml`, the RNA part needs these changes:

```yaml
output_dir: "./soundlife_output"
run_atac_pipeline: false
rna_count_data_path: "./data/soundlife_RNA.h5ad"

# keep the AIFI labels instead of clustering
rna_existing_cell_types: true
rna_celltype_col: "cell_type"
rna_n_target_cell_clusters: null       # a number would merge the labels down to that many

# replace every sev.level default
rna_grouping_columns: ["age_group"]
rna_trajectory_col: "age_at_draw"      # numeric, for the supervised (CCA) trajectory
rna_trajectory_visualization_label: ["age_group"]
rna_cluster_differential_gene_group_col: "age_group"
rna_association_categorical_cols: ["age_group", "sex", "cmv"]
rna_association_continuous_cols: ["age_at_draw"]

rna_trajectory_dge: false              # slow on 32 samples
```

```bash
sampledisco -m complex --config my_config.yaml
```

The run takes about 2 minutes and 4 GB of RAM on 8 CPU cores. It keeps the 17 labels, runs the sample-level Harmony on the 19 batches, and writes the usual folders, now named after `age_group` (`visualization/proportion_embedding_age_group.png`, `sample_cluster/proportion_test/proportion_test_Older_vs_Young.csv`). Age at draw explains about half of the variance of PC2 in `sample_association/variance_explained_sample.csv`. With only the data path changed, the same run stops in the trajectory step with `Trajectory column 'sev.level' not found in pseudo_adata.obs`.

## Existing cell-type labels

- **Use them as they are.** Set `rna_cell_type_cluster: false` and `rna_celltype_col` to your label column. Cell typing is skipped and the embedding reads your column directly. There is no `cell_type.csv` and no cell-level UMAP.
- **Pass them through cell typing.** Set `rna_cell_type_cluster: true`, `rna_existing_cell_types: true` and `rna_celltype_col` to your column, as above. Your labels are copied to `obs['cell_type']`, and the UMAP and `cell_type.csv` are written. With `rna_n_target_cell_clusters` smaller than the number of labels, they are merged into that many groups by a dendrogram of their mean embeddings, with a warning; your original column stays untouched.

Either way, once cell typing has run, the embedding and every downstream step use `obs['cell_type']`. Multi-omics ignores `multiomics_existing_cell_types`; to keep labels there, set `multiomics_cell_type_cluster: false` with the labels in `multiomics_celltype_col`.

## Small cohorts

SampleDisco runs on any number of samples, but several steps have hard limits that bite below a few dozen. The 8-sample COVID demo is at or below most of them, which is why its p-values and predictions are not meaningful.

| Step | Limit | Below it |
|---|---|---|
| Sample embedding | At most **n_units − 1** PCs. | 8 samples give 7 PCs. |
| Sample-level Harmony | Needs **≥ 8 units and > 1 batch value**. | Skipped silently: no batch correction. |
| Preprocessing | `rna_min_cells` is also the minimum number of cells **per sample**. | Smaller samples are dropped (the log prints `Samples remaining`). |
| CCA trajectory | Keep `*_n_cca_pcs` **well below n_samples − 1**. | CCA fits any labelling almost perfectly: the score is near 1 and the `*_cca_pvalue` permutation test is meaningless (it warns). |
| Trajectory DGE (GAM) | A gene is fitted only if it is non-zero in at least order + splines + 2 samples, after automatic reduction. | Below 7 samples nothing can be fitted; with 7–8, only with `*_spline_order: 2`. |
| TSCAN | `*_tscan_n_clusters` ≤ number of samples. | `ValueError: n_clusters=… cannot exceed n_samples=…` |
| Clustering, proportion test, RAISIN | Each compared group needs at least 2 samples. | Comparisons with a one-sample group are skipped. On 8 samples use `*_cluster_number` ≤ 3. |

## Common errors

| Message | Meaning and fix |
|---|---|
| `ValueError: Unexpected parameter in config: '<key>'` | A misspelled or retired key. Check it in the [Config key reference](config_reference.md). |
| `ValueError: Missing required parameter in config: '<key>'` | A key was deleted. Restore it from a fresh `--init-config` template. |
| `ValueError: Trajectory column 'sev.level' not found in pseudo_adata.obs.` | `*_trajectory_col` still has the demo value. |
| `ValueError: Missing required columns in pseudobulk AnnData: ['sev.level']` | `*_grouping_columns` names a column your data does not have (`null` also means `sev.level`). |
| `Warning: Grouping columns not found: {'sev.level'}`, then `Skipping group distance check …` | The same cause in the distance step. Only a warning: the distances are written, the group score is not. |
| `KeyError: The following variables are missing from adata.obs: [...]` | A column in `*_cell_level_batch_key` or `*_sample_level_batch_col` does not exist. |
| `ValueError: Cell type column 'cell_type' not found in adata.obs.` | Cell typing is off and the file has no labels in `*_celltype_col`. |
| `ValueError: Preprocessed data not found at …/preprocess/adata_preprocessed.h5ad.` | A resumed run (`*_preprocessing: false`) without a finished preprocessing; see [Resuming a run](outputs.md#resuming-a-run). |
| `KeyError: … not in adata.obsm` with `Z_clust` or `X_pca_harmony` listed | A file from an older release with a non-standard key passed explicitly. Leave the key `null`, or see [Legacy embedding key names](../legacy_keys.md). |
| `OSError: Unable to synchronously open file (file signature not found)` | The `.h5ad` is not an HDF5 file, typically an error page saved by `curl` without `-f`. Download it again. |

The first two come from the CLI's check, which runs before any data is read.

## Run-to-run variation { #run-to-run-variation }

- **The same machine type gives identical results.** Rerunning the same input on the same CPU type with the same number of CPUs gives bit-identical results, whether you installed from an environment file or with pip. Every install runs cell-level Harmony on harmony-pytorch, and the steps that multithreading made non-deterministic (PCA, ATAC LSI, sample-level Harmony) run single-threaded.
- **Different CPU generations can differ slightly**, for example AVX2 against AVX-512. Leiden then amplifies those last-digit differences into a slightly different cell-type count, so compare such runs by their conclusions, not their decimals.
- **GPU runs** (RAPIDS neighbours and Leiden) are not deterministic between runs.

## Memory

- **RAISIN (`*_cluster_dge`) is the heaviest step.** On the RNA demo it raised the peak memory from about 3.5 GB to 15 GB, because its ComBat step works on a dense genes × cells matrix. The template leaves it off; turn it on for RNA first, and give it at least 16 GB for a cohort the size of the demo.
- Preprocessing holds the full count matrix plus a normalized copy: plan for a few times the size of `.X`.
- On a GPU, `large_data_need_extra_memory: true` lets a dataset spill beyond GPU memory ([GPU acceleration](gpu.md)).
