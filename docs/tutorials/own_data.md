# Using your own data and troubleshooting

This page covers what changes when you move from the demo data to your own cohort: what the input files must contain, how to check them, how to bring in sample metadata and existing cell-type labels, which demo settings you must replace, where small cohorts hit hard limits, and what the common error messages mean. It assumes you have run the [Quick start](quickstart.md) and read the [Configuration guide](configuration.md).

## 1. What the input must contain

Each modality is one cell-level `.h5ad` file (cells × features) holding **all samples** of the cohort.

| Part of the file | Requirement |
|---|---|
| `.X` | **Raw counts**: UMI counts for RNA, fragment or peak counts for ATAC. SampleDisco normalizes, log-transforms and selects highly variable genes itself (Seurat v3 HVGs, which expect counts). Do not pass log-normalized data. |
| `.var_names` | Gene symbols for RNA. Mitochondrial genes are found by the `MT-` / `mt-` prefix (for `rna_pct_mito_cutoff`) and then removed. Peaks for ATAC. |
| `.obs[<sample column>]` | The sample (donor, specimen) of each cell. Name it **`sample`** if you can (see below). Without it, SampleDisco takes the part of each cell name before the first `:` as the sample ID. |
| `.obs` sample-level columns | The phenotypes you want to analyse (for example disease group, age) and any batch column, constant within a sample. They can also come from a CSV ([section 3](#3-sample-and-cell-metadata-files)). |
| `.obs` cell-type labels | Optional. Without them, SampleDisco clusters the cells ([section 5](#5-existing-cell-type-labels)). |

!!! tip "Why a column literally named `sample` is safest"
    `rna_sample_col` (and `atac_sample_col`, `multiomics_sample_col`) can name any column, but several functions you may call directly (sample distance on a cell-level object, for example) fall back to `obs['sample']`. Copying the column once, `adata.obs["sample"] = adata.obs["donor"]`, avoids surprises. With a sample-metadata CSV, 0.6.0 keeps your column and adds a `sample` copy itself. **In 0.5.0** the merge renamed your column to `sample`, and preprocessing then failed with `KeyError: '<your name>'`.

## 2. Inspect a file with `summarize_h5ad`

`summarize_h5ad` prints everything the checklist above asks about: shape, `.obs` / `.var` columns with their types and top values, the value range of `.X` and whether it is integer-like, layers and `.obsm` keys. Here it is on the RNA demo file:

```python
from sampledisco.utils.inspector import summarize_h5ad

summarize_h5ad("data/test_RNA.h5ad", n_examples=3, print_unique_column="sample")
```

```text
🔍 Loading AnnData from: data/test_RNA.h5ad

📦 Basic Info
  - Shape (cells × genes): 29989 × 14881
  - Layers: []
  - obs columns: ['celltype', 'sample', 'old', 'large', 'batch', 'sev.level']
  - var columns: []
  - obsm keys: []
  - varm keys: []

============================================================
🔎 Unique values for column 'sample' (found in adata.obs):
============================================================
  - n_unique: 8
     1. 'CoV-3-Wilk'
     2. 'CoV-34-Aruna'
     3. 'CoV-36-Aruna'
     4. 'CoV-4-Wilk'
     5. 'HD-1-Wilk'
     6. 'HD-20-Aruna'
     7. 'HD-21-Aruna'
     8. 'HD-4-Wilk'

🔎 Inspecting .X matrix:
  - storage type: sparse (csc_matrix)
  - dtype: float64
  - .X nonzeros: 43474639 / 446266309 (9.7419% nonzero, 90.2581% zero)
  - value type: float — values appear integer-like (no fractional parts)
  - contains NaN: False
  - contains Inf: False
  - min value: 1
  - max value: 30539
…
============================================================
📋 Detailed obs columns inspection:
============================================================

    [celltype]
      dtype: category
      non-NA: 29989/29989 (100.0%)
      unique: 26
      top values:
        'T_c2': 4820 (16.1%)
        …
    [batch]
      dtype: category
      non-NA: 29989/29989 (100.0%)
      unique: 2
      top values:
        'Wilk': 17925 (59.8%)
        'Aruna': 12064 (40.2%)

    [sev.level]
      dtype: int64
      non-NA: 29989/29989 (100.0%)
      unique: 2
      min: 1, max: 4, mean: 2.39081
      example values: [np.int64(4), np.int64(1)]
…
```

The output is shortened at `…`: the full printout also shows example `.X` values, cell names, the other `.obs` columns (`sample`, `old`, `large`) and `.var`.

What to look for:

- **`value type: … integer-like`** and a large maximum (here 30,539): `.X` holds raw counts. Fractional values or a maximum around 5–10 mean log-normalized data. The `min value` line ignores the zeros of a sparse matrix.
- **`print_unique_column="sample"`** lists every sample ID: check the count and look for typos that split one sample in two.
- In the `.obs` inspection, each candidate sample-level column should have a handful of unique values. On the demo, `batch` (2 values) and `sev.level` (2 values) are sample-level; `celltype` is a cell-level label.
- `obsm keys` is empty for a raw file. A file that already has `Z_comp` / `Z_rmd` (or the older `X_pca_harmony` / `Z_clust` names) has been preprocessed by SampleDisco before.

The call loads the whole file into memory and returns nothing. On the 30,000-cell demo file it took 5 s and 1.3 GB of RAM.

## 3. Sample and cell metadata files

When the phenotypes are not in `.obs`, put them in a **sample-metadata CSV** and set `rna_sample_meta_path` (`atac_sample_meta_path`; for multi-omics `multiomics_rna_sample_meta_file` and `multiomics_atac_sample_meta_file`).

- One row per sample, a header row, and one column named like `*_sample_col` (normally `sample`) holding the same IDs as `.obs`.
- The separator is detected automatically (comma, tab or semicolon), and a UTF-8 byte-order mark from Excel is handled.
- It is a left join on the sample column: every cell gets its sample's values. Columns that already exist in `.obs` are **replaced** by the CSV version. Empty text fields in the CSV become `"Unknown"`; cells whose sample is not in the CSV get missing values.
- If the sample column is not called `sample`, it is kept and a copy named `sample` is added.

This builds such a file from the demo and merges it back with the function the pipeline uses:

```python
import anndata as ad
from sampledisco.utils.merge_sample_meta import merge_sample_metadata

adata = ad.read_h5ad("data/test_RNA.h5ad")
meta = adata.obs.groupby("sample", observed=True)[["batch", "sev.level"]].first()
meta.to_csv("sample_meta.csv")
print(open("sample_meta.csv").read())

adata.obs = adata.obs.drop(columns=["batch", "sev.level"])   # pretend they were missing
adata = merge_sample_metadata(adata, "sample_meta.csv", sample_column="sample")
print(adata.obs[["sample", "batch", "sev.level"]].drop_duplicates().to_string(index=False))
```

```text
sample,batch,sev.level
CoV-3-Wilk,Wilk,4
CoV-4-Wilk,Wilk,4
CoV-34-Aruna,Aruna,4
CoV-36-Aruna,Aruna,4
HD-1-Wilk,Wilk,1
HD-4-Wilk,Wilk,1
HD-20-Aruna,Aruna,1
HD-21-Aruna,Aruna,1

   📄 Reading CSV metadata file: sample_meta.csv
   ✅ Added 2 new metadata columns
   🔗 Matched 29989/29989 entries (100.0%)
      sample batch  sev.level
  CoV-3-Wilk  Wilk          4
  CoV-4-Wilk  Wilk          4
   HD-1-Wilk  Wilk          1
   HD-4-Wilk  Wilk          1
CoV-36-Aruna Aruna          4
 HD-20-Aruna Aruna          1
 HD-21-Aruna Aruna          1
CoV-34-Aruna Aruna          4
```

In a config, `rna_sample_meta_path: "sample_meta.csv"` does the same merge at the start of preprocessing. Check the `Matched … entries` line in the log: anything below 100% means some cells' sample IDs are not in the CSV, and those cells get empty metadata.

A **cell-metadata CSV** (`rna_cell_meta_path`, `atac_cell_meta_path`) adds per-cell columns. It needs a column named **`barcode`** whose values match `obs_names`; the other columns are joined onto the matching cells. It is read before the sample CSV, so it can supply the sample column itself.

## 4. Adapting the demo config

Start from a fresh template (`sampledisco --init-config my_config.yaml`). The keys you must look at are the data paths, the batch columns, the QC thresholds and the six `sev.level` keys listed in the [Configuration guide](configuration.md#configuration-guide). As an example, here is the change for the RNA part of the [Sound Life subset](demo_data.md) (32 donors, 16 young and 16 older adults, `soundlife_RNA.h5ad`), whose `.obs` has `sample`, `batch`, `age_group`, `age_at_draw`, `sex`, `cmv` and cell-type labels in `cell_type`, but no `sev.level`:

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

rna_trajectory_dge: false              # long on 32 samples; see the Sound Life tutorial
```

The batch keys stay as they are: the template already uses `batch` for both cell-level and sample-level Harmony, and Sound Life has a `batch` column. Then run it:

```bash
sampledisco -m complex --config my_config.yaml
```

On 8 CPU cores this took 2 min 31 s with a peak of 3.7 GB RAM.

Writes → `soundlife_output/rna/` with the same folders as the demo run (`preprocess/`, `sample_embedding/`, `Sample_distance/`, `CCA/`, `sample_cluster/`, `visualization/`, `sample_association/`), now named after `age_group`: `visualization/proportion_embedding_age_group.png`, `sample_cluster/proportion_test/proportion_test_Older_vs_Young.csv`.

- Cell typing kept the 17 labels of `cell_type` (`[cell_types] Found existing cell type annotation.`) and only added the UMAP.
- Sample-level Harmony ran on the 19 processing batches (`[batch correction] 19 groups`).
- Age at draw explained 48% of the variance of PC2 (Pearson r = −0.69, permutation p = 0.01 with the template's 99 permutations), in `sample_association/variance_explained_sample.csv`.
- The proportion test compared `Older` with `Young` (not k-means clusters, because `rna_cluster_differential_gene_group_col` is set) and found one cell type at FDR < 0.01: CD8aa, lower in older donors (logFC −1.24).

The [Sound Life RNA tutorial](soundlife_rna.md) runs this cohort in full, including trajectory DGE and phenotype prediction.

With only `output_dir`, `run_atac_pipeline` and the data path changed, the same run gets through preprocessing, cell typing and the embedding, warns `Grouping columns not found: {'sev.level'}` and `Skipping group distance check …` in the distance step, and then stops in the trajectory step after about 2 min 40 s:

```text
RNA pipeline failed: Trajectory column 'sev.level' not found in pseudo_adata.obs.
```

## 5. Existing cell-type labels

There are two ways to use labels you already have.

- **Use them as they are.** Set `rna_cell_type_cluster: false` and `rna_celltype_col` to your label column. Cell typing is skipped and the sample embedding reads your column directly. You get no `preprocess/cell_type.csv` and no cell-level UMAP.
- **Pass them through cell typing.** Set `rna_cell_type_cluster: true`, `rna_existing_cell_types: true` and `rna_celltype_col` to your column, as in the example above. Your labels are copied to `obs['cell_type']`, optionally merged down to `rna_n_target_cell_clusters` groups by a dendrogram of their mean embeddings, and the UMAP and `cell_type.csv` are written.

!!! note "Cell typing always writes to `obs['cell_type']`"
    Whatever `*_celltype_col` says, the labels from cell typing (copied, merged or clustered) land in `obs['cell_type']`, and once cell typing has run the embedding and every downstream step use that column. On the demo with `rna_existing_cell_types: true`, `rna_celltype_col: "celltype"` (26 labels) and `rna_n_target_cell_clusters: 8`, `obs['cell_type']` held the 8 merged groups and the embedding was built on them (`uns['sample_embedding_params']['K_c'] == 8`); the original labels stay in `celltype`. With `rna_cell_type_cluster: false` and the same column, the embedding used all 26 labels directly.

    **In 0.5.0** the embedding kept reading `*_celltype_col`, so with labels in another column the merge had no effect: the same setting used all 26 original labels. On 0.5.0, put your labels in a column named `cell_type` before merging them.

Multi-omics ignores `multiomics_existing_cell_types` and `multiomics_n_target_clusters`. To keep labels there, set `multiomics_cell_type_cluster: false` with the labels already in `multiomics_celltype_col`.

## 6. Small cohorts

SampleDisco works on any number of samples, but several steps have hard limits that bite below a few dozen. The COVID demo (8 samples) is at or below most of them, which is why its p-values and predictions are not meaningful; the 32-donor [Sound Life RNA tutorial](soundlife_rna.md) is the better reference for results.

| Step | Limit | What happens below it |
|---|---|---|
| Sample embedding | At most **n_units − 1** PCs, whatever `*_sample_embedding_pca_components` says. | 8 samples give 7 PCs. |
| Sample-level Harmony | Needs **≥ 8 units and > 1 batch value**. | Skipped silently: the embedding is not batch-corrected. |
| Preprocessing | `rna_min_cells` is also the minimum number of cells **per sample**. | Samples with fewer cells are dropped (the log prints `Samples remaining`). |
| CCA trajectory | Fits `*_n_cca_pcs` PCs (capped at the number available) to one label per sample. Keep it small and **well below n_samples − 1**. | With PCs close to n_samples − 1, CCA fits any labelling almost perfectly: the score is near 1 and every permuted score is too. |
| CCA p-value (`*_cca_pvalue`) | Same PCs as the observed score. | 0.6.0 warns that the permutation null is degenerate and still writes a p-value, but with `*_n_cca_pcs ≥ n_samples − 1` that p-value only reflects floating-point noise. **In 0.5.0** the step crashes with `ValueError: Too many bins for data range. Cannot create 30 finite-sized bins`. |
| Trajectory DGE (GAM) | A gene is fitted only if it is non-zero in at least *order + splines + 2* samples, after the spline settings are reduced for small cohorts (table below). | Fewer than 7 samples: no gene can be fitted. 7–8 samples: only with `*_spline_order: 2`. |
| TSCAN | `*_tscan_n_clusters` ≤ number of samples. | Raises `ValueError: n_clusters=… cannot exceed n_samples=…`. |
| Sample k-means, proportion test, RAISIN | Each compared group needs at least 2 samples. | Comparisons involving a single-sample cluster are skipped. On 8 samples use `*_cluster_number` ≤ 3. |

Minimum number of samples in which a gene must be non-zero to be fitted by trajectory DGE, with the default `*_num_splines: 5` (computed with `calculate_optimal_spline_parameters`, the function the pipeline calls):

| Samples with a pseudotime | `*_spline_order: 2` | `*_spline_order: 3` |
|---|---|---|
| 4–6 | 7: nothing fits | 7 (reduced to order 2): nothing fits |
| 7–8 | 7 | 9: nothing fits |
| 9 | 7 | 9 |
| 10 | 8 | 9 |
| 11 | 9 | 9 |
| 12–20 | 9 | 10 |

## 7. Common errors

| Message | Meaning and fix |
|---|---|
| `ValueError: Unexpected parameter in config: '<key>'` | A misspelled or retired key. Compare it with the [Config key reference](config_reference.md). |
| `ValueError: Missing required parameter in config: '<key>'` | A key was deleted. Every key must stay in the file (except the three `*_save_cell_adata_after_embedding`); restore it from a fresh `--init-config` template. |
| `ValueError: Trajectory column 'sev.level' not found in pseudo_adata.obs.` | `*_trajectory_col` still has the demo value. See [section 4](#4-adapting-the-demo-config). |
| `ValueError: Missing required columns in pseudobulk AnnData: ['sev.level']` | `*_grouping_columns` names a column your data does not have (a `null` value also means `sev.level`). Raised by the visualization step. |
| `Warning: Grouping columns not found: {'sev.level'}`, then `Skipping group distance check for sample_DR: no grouping column available in adata.obs …` | The same cause in the sample-distance step. Only a warning: the distances are still written, the group-separation score is not. In 0.5.0 the second line reads `Warning: Distance check failed for sample_DR: 'NoneType' object is not iterable`. |
| `KeyError: The following variables are missing from adata.obs: [...]` | A batch column in `*_cell_level_batch_key` or `*_sample_level_batch_col` does not exist. |
| `ValueError: Cell type column 'cell_type' not found in adata.obs.` | Cell typing is off and the file has no labels in `*_celltype_col`. |
| `ValueError: Preprocessed data not found at …/preprocess/adata_preprocessed.h5ad.` | A resumed run (`*_preprocessing: false`) without a finished preprocessing; see [Resuming a run](configuration.md#9-resuming-a-run). |
| `KeyError: comp embedding key 'Z_comp' not in adata.obsm (available: [..., 'Z_clust', ...])` | A file written by an older release, with an embedding key passed explicitly, on 0.5.0. 0.6.0 translates an explicit key to its legacy alias with a `FutureWarning`. Leave the key `null` so it resolves automatically, or rename the keys in memory with `migrate_obsm_keys` (below). |
| `TypeError: 'module' object is not callable` | `from sampledisco import wrapper` on 0.5.0. Use `from sampledisco.wrapper.wrapper import wrapper`, which works on every release. |
| `KeyError: '<your sample column>'` during preprocessing | 0.5.0 only: a sample-metadata CSV with `*_sample_col` other than `sample` ([section 1](#1-what-the-input-must-contain)). Upgrade, or name the column `sample`. |
| `ValueError: Too many bins for data range. Cannot create 30 finite-sized bins` | 0.5.0 only: CCA p-value on a cohort too small for `*_n_cca_pcs` ([section 6](#6-small-cohorts)). Set `*_cca_pvalue: false` or fewer PCs. |
| `OSError: Unable to synchronously open file (file signature not found)` | The `.h5ad` is not an HDF5 file, typically an HTML error page saved by `curl` without `-f`. Download again. |

The first two come from the CLI's check, which runs before any data is read:

```bash
sed 's/^rna_min_genes:/rna_min_gene:/' my_config.yaml > typo.yaml
sampledisco --config typo.yaml
grep -v '^rna_umap:' my_config.yaml > missing.yaml
sampledisco --config missing.yaml
```

```text
Error: Unexpected parameter in config: 'rna_min_gene'
Traceback (most recent call last):
  ...
ValueError: Unexpected parameter in config: 'rna_min_gene'
Error: Missing required parameter in config: 'rna_umap'
Traceback (most recent call last):
  ...
ValueError: Missing required parameter in config: 'rna_umap'
```

Both commands exit with status 1 before any data is loaded (about 15 s each, mostly importing the pipeline).

Old embedding key names can be renamed in memory. `migrate_obsm_keys` copies each legacy key to its current name, keeps the old one, and writes no file:

```python
import anndata as ad
from sampledisco.utils.embedding_keys import migrate_obsm_keys

adata = ad.read_h5ad("data/test_multiomics_integrated.h5ad")
print(list(adata.obsm))
migrate_obsm_keys(adata)
print(list(adata.obsm))
```

```text
['X_glue', 'X_umap', 'Z_clust', 'Z_cmd']
[sampledisco] migrated obsm['Z_clust'] → obsm['Z_comp']
[sampledisco] migrated obsm['Z_cmd'] → obsm['Z_rmd']
['X_glue', 'X_umap', 'Z_clust', 'Z_cmd', 'Z_comp', 'Z_rmd']
```

`test_multiomics_integrated.h5ad` uses the pre-0.3.0 names (see [Legacy embedding key names](../legacy_keys.md)). Save the object (`adata.write_h5ad(...)`) to keep the new keys.

## 8. Memory and runtime

Measured on 8 CPU cores (wall time; peak resident memory of the largest process):

| Run | Data | Time | Peak RAM |
|---|---|---|---|
| `summarize_h5ad` | COVID RNA, 30,000 cells × 14,881 genes | 5 s | 1.3 GB |
| Demo template, RNA + ATAC | COVID RNA + ATAC (29,000 cells × 230,000 peaks) | 3 min – 8 min 40 s, depending on the node | 7.1 GB |
| RNA only, TSCAN ([configuration guide, section 8](configuration.md#8-a-worked-change-unsupervised-trajectory-cluster-based-tests)) | COVID RNA | 1 min 33 s | 3.3 GB |
| Resumed RNA run **with RAISIN** ([section 9](configuration.md#9-resuming-a-run)) | COVID RNA | 1 min 20 s | **15.0 GB** |
| Sound Life RNA, [section 4](#4-adapting-the-demo-config) | 32 donors, 32,000 cells × 18,456 genes | 2 min 31 s | 3.7 GB |

- **RAISIN (`*_cluster_dge`) is the heaviest step.** On the RNA demo it raised the peak memory from 3.3 to 15.0 GB. It runs a ComBat correction and then fits every gene on the cell-level data in parallel workers. 0.6.0 caps the workers at the Slurm allocation (`SLURM_CPUS_PER_TASK`); in 0.5.0 it started one per CPU the job could see, which was 16 on our 8-core allocation (peak 14.8 GB, so most of the peak is in the main process, not the workers). The demo template leaves RAISIN off for all three modalities; its comment notes that ATAC RAISIN (29,000 cells × 230,000 peaks) was killed for lack of memory with 48 cores. Turn it on for RNA first, and give it at least 16 GB for a cohort the size of the demo.
- Preprocessing holds the full count matrix plus a normalized copy; plan for a few times the size of `.X`.
- On a GPU, `use_gpu: true` moves Harmony, the neighbor graph, Leiden and UMAP onto the GPU (the embedding's k-means stays on the CPU), and `large_data_need_extra_memory: true` lets a dataset spill beyond GPU memory. See [GPU acceleration](gpu.md).
