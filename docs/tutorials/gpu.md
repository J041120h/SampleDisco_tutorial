# GPU acceleration

SampleDisco can run its cell-level steps (Harmony, the neighbor graph, Leiden and UMAP) and parts of the sample embedding on an NVIDIA GPU, through RAPIDS and PyTorch. It is off by default. Install the GPU environment as described in [Installation](../installation.md), then switch it on in the config.

## Switch it on

```yaml
use_gpu: true                         # master switch; the template has false
multiomics_use_gpu: true              # multi-omics only, and only together with use_gpu
large_data_need_extra_memory: false   # true: CUDA managed memory instead of the RMM pool
```

- **`use_gpu: true` is required.** Installing RAPIDS is not enough: with `use_gpu: false`, everything runs on the CPU.
- **Multi-omics needs both switches.** `multiomics_use_gpu: false` keeps multi-omics on the CPU while RNA and ATAC use the GPU.
- **From Python**, pass the same keys to `wrapper()`, or `use_gpu=True` to `compute_sample_embedding` and `harmonize_xglue`.
- **Linux only.** On macOS, `use_gpu: true` is ignored.

## What runs on the GPU

| Stage | On the GPU | Stays on the CPU |
|---|---|---|
| RNA and ATAC preprocessing | both cell-level Harmony passes (harmony-pytorch) | QC, feature selection, normalization, PCA or LSI |
| Cell typing (RNA and ATAC) | neighbor graph, Leiden, UMAP (rapids-singlecell) | the merge to `*_n_target_cell_clusters`, plots |
| Multi-omics | the `harmonize_xglue` pass; neighbors and Leiden on the RNA cells, kNN label transfer to ATAC (cuML), UMAP | plots |
| Sample embedding | soft assignment of cells to the anchors (CuPy), the final PCA (cuML) | the k-means fits, the RMD block, the sample-level Harmony |
| Downstream analyses | none | all of them |

scGLUE training (`multiomics_integration: true`) is separate: it uses a visible CUDA GPU on its own, whatever `use_gpu` says.

## Run and check

The RNA pipeline of the demo on the GPU, from a folder with `data/test_RNA.h5ad` ([Demo data](demo_data.md)). The edits switch ATAC off, switch off the trajectory GAM (CPU-only, and most of the run time), switch the GPU on and give the run its own output folder:

```bash
sampledisco --init-config config.yaml
sed -i.bak \
    -e 's/^run_atac_pipeline: true/run_atac_pipeline: false/' \
    -e 's/^rna_trajectory_dge: true/rna_trajectory_dge: false/' \
    -e 's/^use_gpu: false/use_gpu: true/' \
    -e 's|^output_dir: .*|output_dir: "./gpu_output"|' \
    config.yaml
sampledisco --config config.yaml 2>&1 | tee gpu_run.log
grep -iE "DISABLED|falling back|using the CPU implementation" gpu_run.log || echo "no CPU fallback"
```

Every fallback to the CPU prints a warning, so the last command should print `no CPU fallback`. Then check which backend the run used:

```python
import json
import anndata as ad

with open("gpu_output/sys_log/main_process_status.json") as fh:
    info = json.load(fh)["system_info"]
adata = ad.read_h5ad("gpu_output/rna/preprocess/adata_preprocessed.h5ad")
print("gpu_available:", info["gpu_available"], "| embedding backend:", adata.uns["sample_embedding_params"]["backend"])
```

```text
gpu_available: True | embedding backend: gpu
```

- `gpu_available` in `sys_log/main_process_status.json` is `True` only if `use_gpu` was true, the machine runs Linux and the RAPIDS stack initialized.
- `backend` in `uns['sample_embedding_params']` records which code computed the sample embedding.
- The run writes the same files as a CPU run (see [Outputs](outputs.md)), plus possibly a small `rmm_log.txt` in the working folder.

## When it falls back to the CPU

A missing or broken GPU stack does not stop a run. Each stage falls back to its CPU code:

| Condition | What runs instead | Printed |
|---|---|---|
| `use_gpu: true` on macOS | everything on the CPU | nothing |
| RAPIDS fails to import or initialize (no GPU visible, driver too old) | everything on the CPU; `gpu_available: false` | `GPU acceleration DISABLED — running on CPU. The RAPIDS stack failed to initialize (...)` |
| An error inside GPU cell typing (RNA, ATAC) | CPU cell typing, from the start | `[cell_types] WARNING: GPU cell-typing failed (...); falling back to CPU implementation.` |
| `cupy` or `cuml` missing, or an error in the GPU sample embedding | the CPU embedding; `backend: "cpu"` | `... using the CPU implementation.` or `... falling back to the CPU implementation.` |
| The multi-omics GPU cell typing fails to import | CPU cell typing | `Warning: GPU cell-typing unavailable (...). Falling back to CPU.` |

Errors inside the GPU preprocessing and the multi-omics GPU cell typing are not caught. Running out of GPU memory there, for example, stops the run with the CUDA error. Rerun with `use_gpu: false`, or, if the error comes from a RAPIDS step, try `large_data_need_extra_memory: true`.

## GPU memory

When the GPU stack loads, `wrapper()` routes CuPy's and cuML's memory through RMM, in one of two modes:

- **`large_data_need_extra_memory: false` (default): an RMM memory pool.** It reserves half of the GPU's memory at start-up, so other jobs on the same GPU see that memory as taken.
- **`true`: CUDA managed memory.** Allocations can grow beyond the GPU's memory and page out to host RAM, at some cost in speed. Use it when a large cohort runs out of GPU memory in the RAPIDS steps.

PyTorch (harmony-pytorch, scGLUE) keeps its own memory outside RMM.

## Speed and results

Measured on one NVIDIA V100S (32 GB) with 8 CPU cores, the same configs with `use_gpu` on and off:

| Demo run | GPU | CPU only |
|---|---|---|
| RNA (29,988 cells, trajectory GAM off) | 2 min 12 s | 3 min 4 s |
| Multi-omics from the pre-integrated file (59,010 cells) | 50 s | 2 min 37 s |

- **Speed.** The GPU mostly speeds up cell typing (neighbor graph, Leiden, UMAP), so the gain grows with the number of cells. Preprocessing is mostly CPU work, and on 30,000 cells starting CuPy and cuML costs the GPU sample embedding more than it saves. The first GPU run on a node is slower (4.6 minutes here for the RNA run).
- **Memory.** The RNA run peaked at 4.5 GB of RAM (3.7 GB on the CPU) and 16.7 GB of GPU memory, mostly the RMM pool.
- **Results.** GPU and CPU runs agree in their broad structure but are not identical. On the RNA demo the cell types agreed with an adjusted Rand index of 0.81–0.89 and the sample distances correlated at 0.98; for multi-omics, 0.77 and 0.94 (15 vs 16 RNA clusters). The GPU neighbor graph and Leiden are not deterministic, so two GPU runs also differ (ARI 0.85), while two CPU runs agreed almost exactly (ARI 0.995). Pick one backend per study, and use the CPU when results must repeat.
