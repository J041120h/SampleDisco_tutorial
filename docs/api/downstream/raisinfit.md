# `raisinfit`

Python port of the RAISIN hierarchical model for differential expression between groups of samples, from single-cell data. For each gene it estimates the per-sample means and two variance components: between samples within a group (`sigma2`) and between cells within a sample (`omega2`), with empirical-Bayes shrinkage. With `batch_col`, ComBat is applied first (preserving the group effect). The returned fit is tested with [`run_pairwise_tests`](run_pairwise_tests.md). Designs: `"unpaired"` (groups of samples), `"paired"` (with `individual_col`), `"continuous"` and `"custom"`.

Expression is read from `adata.raw.X` when present, otherwise from `adata.X` (log-normalized after SampleDisco preprocessing). With `batch_col`, ComBat corrects that same matrix (the `.raw` genes when `.raw` exists), and the model is fitted on the corrected values. The matrix is densified to genes × cells for the fit, so memory grows with genes × cells.

**Import:** `from sampledisco.sample_clustering.RAISIN import raisinfit`

**Source:** `sample_clustering/RAISIN.py:219`

## Signature

```python
def raisinfit(
    adata,
    sample_col,
    testtype="unpaired",
    group_col=None,
    individual_col=None,
    batch_col=None,
    sample_to_clade=None,
    custom_design=None,
    intercept=True,
    filtergene=False,
    filtergenequantile=0.5,
    max_features=None,
    n_jobs=None,
    verbose=True,
    seed: int = 42,
)
```

## Parameters

| Name | Type | Default | Description |
| --- | --- | --- | --- |
| `adata` | AnnData | — | Cell-level (or pseudobulk) AnnData. |
| `sample_col` | str | — | Sample column. |
| `testtype` | str | `"unpaired"` | `"unpaired"`, `"paired"`, `"continuous"` or `"custom"`. |
| `group_col` | str or None | `None` | Group (or, for `"continuous"`, covariate) column. Used when present in `.obs`; then `sample_to_clade` is ignored. |
| `individual_col` | str or None | `None` | Subject column for `"paired"`. |
| `batch_col` | str or None | `None` | Run ComBat on this column first (not for `"custom"`). |
| `sample_to_clade` | dict or None | `None` | `{sample_id: group}`, used when `group_col` is not given, e.g. from [`cluster`](cluster.md). |
| `custom_design` | dict or None | `None` | For `"custom"`: keys `"X"`, `"Z"`, `"group"`. |
| `intercept` | bool | `True` | Intercept in the fixed-effect design. |
| `filtergene` | bool | `False` | Drop lowly expressed genes first. |
| `filtergenequantile` | float | `0.5` | Quantile used by `filtergene`. |
| `max_features` | int or None | `None` | Keep only the most variable genes, to bound memory. |
| `n_jobs` | int or None | `None` | Processes; `None` or `-1` = the CPUs available to the process (its CPU-affinity mask, capped at `SLURM_CPUS_PER_TASK` when set). |
| `verbose` | bool | `True` | Print progress. |
| `seed` | int | `42` | Random seed. |

## Returns

`dict`:

| Key | Meaning |
| --- | --- |
| `"mean"` | genes × samples mean expression |
| `"sigma2"` | genes × variance groups, between-sample variance |
| `"omega2"` | genes × samples, within-sample (cell-level) variance |
| `"X"`, `"Z"` | fixed- and random-effect design matrices |
| `"group"` | variance group of each random effect |
| `"failgroup"` | groups whose variance could not be estimated (no residual degrees of freedom, e.g. a group with a single sample); their variance is set to 0, and `run_pairwise_tests` skips every contrast that uses them |
| `"eb_fallback"` | per group, genes whose shrunken variance was not finite and was set to 1 (as in R) |
| `"sample_names"` | sample order |
| `"batch_corrected"` | whether ComBat ran |

## Usage

Continuing from the RNA chain, healthy vs severe:

```python
import anndata as ad
from sampledisco.sample_clustering.RAISIN import raisinfit

adata = ad.read_h5ad("api_demo/rna/preprocess/adata_preprocessed.h5ad")

fit = raisinfit(adata, sample_col="sample", group_col="sev.level", testtype="unpaired", n_jobs=8)
print(sorted(fit))
print(fit["mean"].shape, fit["failgroup"])
```

Output, after the fitting log (`Expression matrix: 9936 genes x 29680 cells`, one variance estimate per group `1` and `4`):

```text
['X', 'Z', 'batch_corrected', 'eb_fallback', 'failgroup', 'group', 'mean', 'omega2', 'sample_names', 'sigma2']
(9936, 8) []
```

About 75 s and 4.5 GB including the worker processes. The log line `Using counts from adata.X` means "using `adata.X`"; here that is the log-normalized matrix from preprocessing.

The next step is [`run_pairwise_tests`](run_pairwise_tests.md).
