# GPU acceleration

SampleDisco can run its cell-level steps on an NVIDIA GPU through [RAPIDS](https://rapids.ai) and PyTorch: Harmony, the neighbor graph, Leiden and UMAP. Parts of the sample embedding can run there too. GPU use is off by default. This page covers:

- what moves to the GPU;
- how to install the GPU environment and switch the GPU on;
- how to tell which backend a run actually used;
- how much a GPU run differs from a CPU run on the COVID-19 demo, in results and in run time.

!!! info "In short"
    - **Switching on:** set `use_gpu: true`. Installing RAPIDS is not enough: with `use_gpu: false`, SampleDisco's GPU code stays off. (Two PyTorch steps can still use a visible GPU; see section 1.)
    - **Multi-omics:** it uses the GPU only when **both** `multiomics_use_gpu` and `use_gpu` are true.
    - **Fallback:** if the GPU stack is missing or broken, the run falls back to the CPU and keeps going. Check `gpu_available` in `sys_log/main_process_status.json` and `uns['sample_embedding_params']['backend']`.
    - **Results:** GPU and CPU results agree in their broad structure but are not identical. Two GPU runs also differ from each other, while CPU runs repeat exactly. Pick one backend per project; use the CPU when results must repeat exactly.
    - **Speed:** on the 8-sample demo the GPU makes cell typing 2× faster for RNA and 6× faster for multi-omics, but the whole RNA run is no faster.

## 1. What runs on the GPU

| Stage | On the GPU with `use_gpu: true` | Still on the CPU | Code |
| --- | --- | --- | --- |
| RNA preprocessing | both Harmony passes (`Z_comp`, `Z_rmd`), through harmony-pytorch on CUDA | QC, HVG selection, normalization, PCA (scanpy) | `preparation/rna_preprocess_gpu.py` |
| ATAC preprocessing | both Harmony passes | TF-IDF, HVF selection, LSI, doublet removal (muon / scanpy) | `preparation/atac_preprocess_gpu.py` |
| Cell typing (RNA and ATAC) | neighbor graph, Leiden, UMAP (rapids-singlecell) | dendrogram merge to `n_target_cell_clusters`, plots | `preparation/cell_type_gpu.py` |
| Multi-omics Harmony post-pass | `harmonize_xglue(use_gpu=True)` (harmony-pytorch on CUDA) | – | `preparation/multi_omics_batch_correction.py` |
| Multi-omics cell typing | neighbor graph + Leiden on the RNA cells (rapids-singlecell), kNN label transfer to ATAC (cuML), UMAP | plots | `preparation/multi_omics_cell_type_gpu.py` |
| Sample embedding | RBF soft assignment of cells to the k-means centres (CuPy), final PCA (cuML) | the k-means fits (see below), the RMD block | `sample_embedding/sample_embedding_gpu.py` |
| Downstream analyses | – | all of them: distances, CCA / TSCAN, GAM DGE, clustering, proportion tests, RAISIN, plots | – |

!!! note "Where the k-means runs"
    The GPU sample embedding asks cuML for a `MiniBatchKMeans`. cuML 24.12, the version in the GPU environment, has no `MiniBatchKMeans`. The code then silently uses the same scikit-learn `MiniBatchKMeans` as the CPU path, with the same seeds. So the k-means centres are computed on the CPU on both backends.

!!! note "Two PyTorch steps choose their own device"
    `use_gpu` does not control these two steps; both use the GPU whenever PyTorch can see one:

    - **The sample-level Harmony** at the end of the sample embedding, on both the GPU and the CPU path. It runs through harmonypy 0.2.0, which picks the device itself.
    - **scGLUE training** (`multiomics_integration: true`): scGLUE also picks its own device.

When `use_gpu: true` and the GPU stack loads, the wrapper also sets up the [RMM](https://docs.rapids.ai/api/rmm/stable/) memory pool for CuPy and cuML (section 9).

## 2. Install the GPU environment

The GPU environment is the CPU environment plus RAPIDS 24.12 and PyTorch 2.5.1 (the CUDA 12.1 wheel). This one environment runs every part of SampleDisco 0.6.0: RNA, ATAC, scGLUE, MaxFuse and the GPU paths. Build it from `environment-gpu.yml` (the same file as in the SampleDisco repository):

```bash
curl -fLO https://j041120h.github.io/SampleDisco_tutorial/assets/environment-gpu.yml
conda env create -f environment-gpu.yml
conda activate sampledisco-gpu
pip install --no-deps rapids-singlecell==0.13.1 sampledisco==0.6.0
```

- **Platform:** Linux x86-64 with an NVIDIA GPU only. macOS and CPU-only machines use the CPU environment (see [Installation](../installation.md)).
- **Driver:** RAPIDS 24.12 is built for CUDA 12.0–12.5 and runs on a driver that supports CUDA 12.5 (we used driver 555.42.06). RAPIDS 25.04 and later need a driver for CUDA 12.6 or newer; on an older driver they fail at import with `cudaErrorInsufficientDriver`. If your driver supports CUDA 12.6 or newer, you can move the RAPIDS pins to 25.x.
- **No GPU on the build machine** (for example a login node): prefix the first command with `CONDA_OVERRIDE_CUDA=12.5`.
- **Keep `--no-deps` on the `pip` line.** The env file already installs every dependency. A dependency-resolving `pip install sampledisco` would replace conda's `pyarrow` with a version that breaks cuDF (RAPIDS 24.12 needs `pyarrow < 18`).
- **Build time:** about 32 minutes on our cluster (4 cores).

Check the stack on a GPU node:

```bash
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
python -c "import cupy, cuml, rapids_singlecell, torch; print('cuML', cuml.__version__, '| rapids-singlecell', rapids_singlecell.__version__, '| torch sees a GPU:', torch.cuda.is_available())"
```

```text
name, driver_version, memory.total [MiB]
NVIDIA L40S, 555.42.06, 46068 MiB
cuML 24.12.00 | rapids-singlecell 0.13.1 | torch sees a GPU: True
```

If the import fails, fix the environment first. SampleDisco never installs packages at run time. With a broken stack it runs on the CPU and prints a warning (section 9).

## 3. Switch it on

Three config keys control the GPU. Their values in the demo template:

```yaml
use_gpu: false                        # master switch; set true to use the GPU
multiomics_use_gpu: true              # multi-omics only, and only if use_gpu is also true
large_data_need_extra_memory: false   # true: CUDA managed memory instead of the RMM pool (section 9)
```

- **`use_gpu: true` is required.** With `use_gpu: false`, SampleDisco's own steps stay on the CPU even when RAPIDS is installed and a GPU is present.
- **Multi-omics needs both switches.** The multi-omics pipeline runs on the GPU only when `multiomics_use_gpu` **and** `use_gpu` are true. `multiomics_use_gpu: false` keeps multi-omics on the CPU while RNA and ATAC use the GPU.
- **From Python:** pass the same keys to `wrapper()`. Or pass `use_gpu=True` to `compute_sample_embedding` and `harmonize_xglue` directly (sections 7 and 8).
- **macOS:** `use_gpu: true` is ignored there (the wrapper checks for Linux).

## 4. Run the demo on the GPU

The rest of the page runs the RNA pipeline of the COVID-19 demo (8 samples) twice, once on the GPU and once on the CPU, and compares the two runs. Work on a GPU node with the GPU environment active, from a new folder:

```bash
mkdir -p sampledisco_gpu/data
cd sampledisco_gpu
curl -fL -o data/test_RNA.h5ad "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
sampledisco --init-config config.yaml
sed -i.bak \
    -e 's/^run_atac_pipeline: true/run_atac_pipeline: false/' \
    -e 's/^rna_trajectory_dge: true/rna_trajectory_dge: false/' \
    -e 's/^use_gpu: false/use_gpu: true/' \
    -e 's|^output_dir: .*|output_dir: "./gpu_output"|' \
    config.yaml
grep -E '^(output_dir|run_[a-z]+_pipeline|rna_trajectory_dge|use_gpu):' config.yaml
```

```text
Wrote the demo config to 'config.yaml'.
Edit the data paths / options, then run:
  sampledisco -m complex --config config.yaml
WARNING: 'config.yaml' is the config used for the Zenodo demo dataset — its QC thresholds and trajectory/grouping columns (e.g. 'sev.level', '*_supervised: true') are tuned for that demo data, not yours. Review and adjust every field before running on your own data.
output_dir: "./gpu_output"
run_rna_pipeline: true
run_atac_pipeline: false
run_multiomics_pipeline: false   # set true to run multi-omics; pre-wired to skip scGLUE training (see multiomics_integration below)
use_gpu: true   # set true to use the RAPIDS GPU paths if installed; auto-falls back to CPU otherwise
rna_trajectory_dge: false
```

The edits do four things:

- switch the ATAC pipeline off;
- switch off the trajectory GAM (`rna_trajectory_dge`), which runs on the CPU either way and would take most of the run time;
- switch the GPU on;
- give the run its own output folder.

Now run the pipeline and keep a copy of the log:

```bash
sampledisco --config config.yaml 2>&1 | tee gpu_run.log
```

Every fallback to the CPU prints a warning. None should appear:

```bash
grep -iE "DISABLED|falling back|using the CPU implementation" gpu_run.log || echo "no CPU fallback"
```

```text
no CPU fallback
```

**Writes** →

- `gpu_output/rna/`: the same files as a CPU run (see [Quick start](quickstart.md));
- `gpu_output/sys_log/main_process_status.json`;
- possibly a small `rmm_log.txt` (RMM's log) in the working folder.

## 5. The same run on the CPU

Copy the config with the GPU switched off and a second output folder. Then run it with the GPU hidden:

```bash
sed -e 's/^use_gpu: true/use_gpu: false/' \
    -e 's|^output_dir: .*|output_dir: "./cpu_output"|' \
    config.yaml > config_cpu.yaml
CUDA_VISIBLE_DEVICES="" sampledisco --config config_cpu.yaml 2>&1 | tee cpu_run.log
```

!!! note "Why hide the GPU"
    `use_gpu: false` switches off SampleDisco's GPU code, but not PyTorch's own device choice. harmonypy 0.2.0 runs the sample-level Harmony on both paths, and it moves to any GPU that PyTorch can see; the log then says `Running Harmony (PyTorch on cuda)`. `CUDA_VISIBLE_DEVICES=""` hides the GPU, so this run computes exactly what a machine without a GPU computes. On a CPU-only machine you don't need it.

    It matters: we also ran the CPU config with the GPU visible. The cell-level results were identical, but the sample embedding changed (the correlation of sample distances fell to 0.979). That came from harmonypy's device alone.

## 6. Compare the two runs

### Which backend ran

```python
import json
import anndata as ad

runs = {}
for name in ["gpu", "cpu"]:
    with open(f"{name}_output/sys_log/main_process_status.json") as fh:
        info = json.load(fh)["system_info"]
    adata = ad.read_h5ad(f"{name}_output/rna/preprocess/adata_preprocessed.h5ad")
    backend = adata.uns["sample_embedding_params"]["backend"]
    print(f"{name}: use_gpu={info['use_gpu']}  gpu_available={info['gpu_available']}  backend={backend}")
    runs[name] = adata
```

```text
gpu: use_gpu=True  gpu_available=True  backend=gpu
cpu: use_gpu=False  gpu_available=False  backend=cpu
```

- **`gpu_available`** in `sys_log/main_process_status.json` is the wrapper's decision for the whole run. It is true only if `use_gpu` was true, the machine runs Linux and the RAPIDS stack initialized.
- **`uns['sample_embedding_params']['backend']`** (`"gpu"` or `"cpu"`) records which sample-embedding code produced `uns['X_DR_sample']`.

### Cell types

```python
import pandas as pd
from sklearn.metrics import adjusted_rand_score

gpu, cpu = runs["gpu"], runs["cpu"]
print("same cells:", gpu.obs_names.equals(cpu.obs_names))
sizes = pd.concat({name: a.obs["cell_type"].value_counts().reset_index(drop=True)
                   for name, a in runs.items()}, axis=1)
print(sizes.T)   # cluster sizes, largest first: Leiden numbers clusters differently in each run
print("ARI:", round(adjusted_rand_score(gpu.obs["cell_type"], cpu.obs["cell_type"]), 3))
```

```text
same cells: True
        0     1     2     3     4     5     6    7
gpu  9192  7797  4037  3946  1995  1234  1203  584
cpu  9604  7991  7105  1993  1189  1178   855   73
ARI: 0.833
```

Both runs keep the same 29,988 cells and merge the Leiden clusters down to `rna_n_target_cell_clusters: 8`, but the partitions differ (adjusted Rand index 0.83). Two things contribute:

- **Different input.** The PCA is computed on the CPU in both runs and is identical. The Harmony output `Z_comp` that the clustering uses differs a little between the GPU and the CPU implementation (per-dimension correlation ≥ 0.9997).
- **Different, non-deterministic algorithm.** The neighbor graph and Leiden run in rapids-singlecell instead of scanpy / igraph, and the GPU versions are not deterministic. A second GPU run of the same config gave the same `Z_comp` but different clusters: ARI 0.92 against the first GPU run. Two CPU runs gave identical clusters.

### Sample embedding

```python
import numpy as np
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist

E_gpu = gpu.uns["X_DR_sample"]
E_cpu = cpu.uns["X_DR_sample"].loc[E_gpu.index]
_, _, disparity = procrustes(E_cpu.to_numpy(), E_gpu.to_numpy())
r_dist = np.corrcoef(pdist(E_gpu), pdist(E_cpu))[0, 1]
r_pc = [abs(np.corrcoef(E_gpu.iloc[:, k], E_cpu.iloc[:, k])[0, 1]) for k in range(3)]
print(f"Procrustes disparity: {disparity:.3f}")
print(f"correlation of sample distances: {r_dist:.3f}")
print("|r| of PC1-PC3:", np.round(r_pc, 3))
```

```text
Procrustes disparity: 0.030
correlation of sample distances: 0.963
|r| of PC1-PC3: [0.91  0.908 0.968]
```

- **Procrustes disparity** is 0 for embeddings that are identical up to rotation, scaling and shift, and 1 for unrelated ones. 0.03 means the eight samples sit in nearly the same configuration.
- **Sample distances** correlate at 0.963 between the two runs, and the leading PCs at 0.91–0.97. PCA signs are arbitrary, so the correlations are taken as absolute values.
- **Downstream results** follow the embedding: the CCA pseudotime orders the samples almost the same way (Spearman |ρ| = 0.976).
- **Two GPU runs** differ from each other by about as much as GPU from CPU: sample-distance correlation 0.975, pseudotime |ρ| = 0.88. Two CPU runs were bit-identical.

With 8 samples these numbers are noisy. For a cohort large enough for stable sample-level statistics, see the [Sound Life tutorial](soundlife_rna.md).

### Run time

We measured both runs on the same node: an NVIDIA L40S (46,068 MiB) and 8 CPU cores. The stage times come from time-stamped logs, and peak RAM from `/usr/bin/time -v`. The "repeat" columns rerun the same configs.

| Stage | GPU | GPU, repeat | CPU (GPU hidden) | CPU, repeat |
| --- | --- | --- | --- | --- |
| Preprocessing (QC, HVGs, PCA, 2 × Harmony) | 56 s | 37 s | 60 s | 51 s |
| Cell typing (neighbors, Leiden, UMAP, plots) | 24 s | 22 s | 47 s | 47 s |
| Sample embedding | 27 s | 20 s | 3 s | 2 s |
| Downstream (on the CPU in both runs) | 12 s | 10 s | 10 s | 10 s |
| **Whole command (wall clock)** | **2 min 14 s** | **1 min 40 s** | **2 min 11 s** | **2 min 1 s** |
| Peak RAM | 4.3 GB | 4.2 GB | 3.4 GB | 3.4 GB |
| Peak GPU memory (`nvidia-smi`) | 23,561 MiB | 23,397 MiB | – | – |

On this 30,000-cell demo the GPU halves cell typing, but the whole run is no faster:

- **Preprocessing** is mostly CPU work (QC, HVG selection, PCA) on both paths.
- **The GPU sample embedding is slower:** 20–27 s against 2–3 s. Starting CuPy and cuML costs more than the GPU saves on 30,000 cells, and the k-means runs on the CPU anyway (section 1).
- **Run-to-run variation** of 10–20 s comes from file reads on our shared file system.
- **GPU memory** is mostly the RMM pool, reserved at start-up (section 9).

The multi-omics demo (section 8), with twice as many cells, gains more.

## 7. GPU code from Python

`compute_sample_embedding` takes `use_gpu` directly. Running it both ways on the same input separates the sample-embedding code from everything upstream. Continue in the same Python session:

```python
from sampledisco.sample_embedding import compute_sample_embedding

def embed(use_gpu, batch_method):
    a = runs["cpu"].copy()   # the same cells, cell types, Z_comp and Z_rmd every time
    compute_sample_embedding(a, "embedding_check", use_gpu=use_gpu, sample_col="sample",
                             celltype_col="cell_type", batch_col="batch",
                             batch_method=batch_method, save=False, verbose=False)
    return a.uns["X_DR_sample"].to_numpy(), a.uns["sample_embedding_params"]["backend"]

for batch_method in ["none", "harmony"]:
    (G, g), (C, c) = embed(True, batch_method), embed(False, batch_method)
    G = G * np.sign((G * C).sum(axis=0))   # PCA signs are arbitrary: align each PC
    print(f"batch_method={batch_method!r}: backends {g}/{c}, largest |GPU - CPU| = "
          f"{np.abs(G - C).max():.1e} (largest value {np.abs(C).max():.1f})")
```

```text
batch_method='none': backends gpu/cpu, largest |GPU - CPU| = 1.3e-05 (largest value 12.0)
batch_method='harmony': backends gpu/cpu, largest |GPU - CPU| = 1.1e+00 (largest value 12.0)
```

With `save=False` nothing is written, and the output folder argument is not used.

On identical input, the GPU and CPU sample-embedding code agree to rounding error (1e-5) up to the final PCA. The sample-level Harmony then turns that rounding difference into a visible one: up to 1.1 on values of up to 12. Harmony on 8 samples is that sensitive to its input. This is also why hiding the GPU from harmonypy (section 5) changes the result. So most of the GPU-versus-CPU difference in section 6 comes from upstream (Harmony on the cells, GPU cell typing), and the sample-level Harmony enlarges it.

## 8. Multi-omics on the GPU

The multi-omics pipeline can start from the demo's pre-integrated scGLUE file. Download it into the same folder, switch RNA off and multi-omics on, and run:

```bash
curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
sed -e 's/^run_rna_pipeline: true/run_rna_pipeline: false/' \
    -e 's/^run_multiomics_pipeline: false/run_multiomics_pipeline: true/' \
    config.yaml > config_mo.yaml
sampledisco --config config_mo.yaml 2>&1 | tee mo_gpu_run.log
grep -E "GPU mode|sample_embedding_gpu\] done" mo_gpu_run.log
```

```text
GPU mode: Enabled
[sample_embedding_gpu] done in 8.99s; shape=(16, 10)
```

`GPU mode: Enabled` appears only when `multiomics_use_gpu` and `use_gpu` are both true. The sample embedding has 16 units: each sample's RNA and ATAC cells. The demo file already contains a composition embedding, so the pipeline skips the Harmony post-pass and goes straight to cell typing and the sample embedding.

**Writes** → `gpu_output/multiomics/`. SampleDisco 0.6.0 writes the union with the new cell types and `uns['X_DR_sample']` to `gpu_output/multiomics/preprocess/adata_sample.h5ad` and leaves the downloaded file unchanged.

We ran the same config with `use_gpu: false` (GPU hidden) for comparison, on the same L40S node with 8 CPU cores:

| Multi-omics demo (59,010 cells) | GPU | CPU |
| --- | --- | --- |
| Cell typing (Leiden on 29,989 RNA cells, label transfer to 29,021 ATAC cells, UMAP, plots) | 13 s | 84 s |
| Sample embedding | 11 s | 5 s |
| Downstream | 12 s | 18 s |
| **Whole command (wall clock)** | **56 s** | **2 min 2 s** |
| Peak RAM | 2.5 GB | 2.6 GB |

Results differ more here than for RNA:

- **Cell types:** the GPU Leiden found 14 RNA clusters and the CPU Leiden found 16, at the same resolution (0.8). The partitions agree with an adjusted Rand index of 0.75.
- **Sample embeddings:** the sample distances correlate at 0.94.

`harmonize_xglue` is the Harmony pass that turns the scGLUE embedding `X_glue` into the sample-removed `Z_comp`. The pipeline skips it for the demo file, as above. Called directly, it runs on the GPU with `use_gpu=True`:

```python
import time
import anndata as ad
import numpy as np
from sampledisco.preparation.multi_omics_batch_correction import harmonize_xglue

mo = ad.read_h5ad("data/test_multiomics_integrated.h5ad")
for key, use_gpu in [("Z_comp_gpu", True), ("Z_comp_cpu", False)]:
    t0 = time.time()
    harmonize_xglue(mo, sample_col="sample", in_key="X_glue", out_key=key,
                    use_gpu=use_gpu, max_iter=50, random_state=42, verbose=False)
    print(f"use_gpu={use_gpu}: {time.time() - t0:.1f} s")
r = [np.corrcoef(mo.obsm["Z_comp_gpu"][:, k], mo.obsm["Z_comp_cpu"][:, k])[0, 1]
     for k in range(mo.obsm["X_glue"].shape[1])]
print(f"{mo.n_obs} cells x {len(r)} dims; GPU vs CPU correlation per dim: min {min(r):.3f}, median {np.median(r):.3f}")
```

```text
Use GPU mode.
	Initialization is completed.
	Completed 1 / 50 iteration(s).
	...
Reach convergence after 7 iteration(s).
use_gpu=True: 7.0 s
	Initialization is completed.
	...
Reach convergence after 7 iteration(s).
use_gpu=False: 22.3 s
59010 cells x 50 dims; GPU vs CPU correlation per dim: min 1.000, median 1.000
```

The iteration lines come from harmony-pytorch. On 59,010 cells the GPU Harmony is 3× faster than the CPU one, and the two corrected embeddings agree (every dimension correlates at 1.000 to three decimals). The GPU time includes starting CUDA in this Python session.

## 9. Fallbacks and GPU memory

### When the run falls back to the CPU

A missing or broken GPU stack does not stop a run. Each stage falls back to its CPU code, and every case below except the first and the last prints a warning. The checks, in the order a run meets them:

| Where | Condition | What runs instead | Printed |
| --- | --- | --- | --- |
| Start of `wrapper()` | `use_gpu: true` on macOS | everything on the CPU; `gpu_available: false` | nothing |
| Start of `wrapper()` | `rmm`, `cupy` or `rapids_singlecell` fails to import, or the RMM set-up fails (no GPU visible, driver too old) | everything on the CPU; `gpu_available: false` | `GPU acceleration DISABLED — running on CPU. The RAPIDS stack failed to initialize (...)` |
| RNA / ATAC pipeline | the `*_gpu` modules fail to import | that pipeline on the CPU | `[sampledisco] WARNING: GPU acceleration DISABLED — running on CPU.` |
| RNA / ATAC cell typing | any error inside the GPU cell typing | the CPU cell typing, from the start | `[cell_types] WARNING: GPU cell-typing failed (...); falling back to CPU implementation.` |
| Sample embedding | `cupy` / `cuml` missing, or any error in the GPU code | the CPU code; `backend: "cpu"` | `... using the CPU implementation.` or `GPU sample embedding unavailable (...); falling back to the CPU implementation.` |
| Multi-omics cell typing | the GPU module fails to import | the CPU cell typing | `Warning: GPU cell-typing unavailable (...). Falling back to CPU.` |
| Cell-level Harmony passes | harmony-pytorch not installed | harmonypy, on whatever device PyTorch sees | nothing |

We tested the second row by running the section 4 config (`use_gpu: true`) with the GPU environment on a node without a GPU. The `grep` from section 4 then prints the warning instead of `no CPU fallback`:

```text
2026-11-23 00:24:59 | [WARNING] GPU acceleration DISABLED — running on CPU. The RAPIDS stack failed to initialize (RuntimeError: Function "cuInit" not found). Verify with `python -c "import cuml, cupy, rapids_singlecell"`; see the installation guide. (This run's backend is recorded as gpu_available in sys_log/main_process_status.json.)
```

The run finished on the CPU (5 min 45 s on that CPU-only node). `main_process_status.json` recorded `use_gpu: true, gpu_available: false`, and the embedding recorded `backend: "cpu"`.

Errors that happen once the GPU code is running are different. The GPU preprocessing (RNA and ATAC) and the multi-omics GPU cell typing do not catch them. For example, running out of GPU memory there stops the run with the CUDA error. Rerun with `use_gpu: false`. If the error comes from a RAPIDS step (neighbor graph, Leiden, UMAP, cuML), `large_data_need_extra_memory: true` may help (next section).

So after a GPU run, check three things:

- no warning in the log (the `grep` in section 4);
- `gpu_available: true` in `sys_log/main_process_status.json`;
- `backend: "gpu"` in `uns['sample_embedding_params']`.

### GPU memory

When the GPU stack loads, `wrapper()` routes CuPy's and cuML's GPU memory through RMM, in one of two modes set by `large_data_need_extra_memory`:

- **`false` (default): an RMM memory pool.** The pool reserves half of the GPU's memory as soon as the run starts. On our L40S (46,068 MiB), 22,965 MiB were in use within seconds of start-up, and the RNA run peaked at 23,561 MiB. Other jobs sharing the GPU see that memory as taken.
- **`true`: CUDA managed (unified) memory, without a pool.** Allocations can grow beyond the GPU's memory and page out to host RAM, at some cost in speed. Use it when a large cohort runs out of GPU memory in the RAPIDS steps. On the RNA demo, the section 4 config with `large_data_need_extra_memory: true` peaked at 1,159 MiB of GPU memory and finished in 1 min 39 s. On data this small, managed memory costs no time.

PyTorch (harmony-pytorch, harmonypy) keeps its own memory outside RMM. Calling the GPU functions directly from Python (sections 7 and 8) sets up no pool; those calls used under 1 GB of GPU memory.
