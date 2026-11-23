# Installation

SampleDisco runs in **one conda environment** that covers every feature on this site:
- the RNA, ATAC and multi-omics pipelines;
- training scGLUE from scratch (with the `bedtools` binary it needs);
- MaxFuse for RNA + protein / flow;
- every downstream analysis.

Pick the **CPU** environment, or the **GPU** environment, which is the same plus NVIDIA RAPIDS. Both install SampleDisco itself, so there is no separate `pip install` to remember.

!!! info "Requirements"
    - **Linux x86-64.** The environment files are validated on Linux. GPU acceleration needs an NVIDIA GPU, and the CPU must support AVX2; see [Troubleshooting](#troubleshooting).
    - **conda** (Miniforge or Miniconda).
    - **RAM, measured on the demo data with 8 cores:**

        | Run | Peak RAM |
        |---|---|
        | RNA pipeline | about 4 GB |
        | RNA + ATAC | about 8.5 GB |
        | RAISIN cluster DGE | about 16 GB |

## 1. CPU environment (recommended)

```bash
curl -fLO https://raw.githubusercontent.com/J041120h/SampleDisco/main/environment-cpu.yml
conda env create -f environment-cpu.yml
conda activate sampledisco-cpu
```

The build takes about 15–25 minutes. Check it:

```bash
python -c "import sampledisco, scglue, maxfuse; print('sampledisco', sampledisco.__version__)"
sampledisco --help
bedtools --version
```

The first line should print `sampledisco 0.6.0`. Next: the [Quick start](tutorials/quickstart.md).

## 2. GPU environment (Linux + NVIDIA)

```bash
curl -fLO https://raw.githubusercontent.com/J041120h/SampleDisco/main/environment-gpu.yml
conda env create -f environment-gpu.yml      # building on a machine without a GPU? prefix CONDA_OVERRIDE_CUDA=12.5
conda activate sampledisco-gpu
pip install --no-deps rapids-singlecell==0.13.1 sampledisco==0.6.0
python -c "import cuml, cupy, rapids_singlecell, torch; print('GPU OK', torch.cuda.is_available())"
```

- The build takes about 30 minutes.
- This is everything in the CPU environment, plus RAPIDS 24.12 and torch 2.5.1+cu121.
- **Keep the `--no-deps` in the last step.** The env file already installs every dependency. A dependency-resolving `pip install sampledisco` would replace the `pyarrow` that RAPIDS needs and break `cudf`.
- RAPIDS 24.12 targets CUDA 12.0–12.5 drivers. With a driver for CUDA ≥ 12.6 you may move the RAPIDS pins to 25.x.
- Installing RAPIDS is not enough to use the GPU: set `use_gpu: true` in the config. See [GPU acceleration](tutorials/gpu.md) for what runs on the GPU and how the CPU fallback works.

## 3. pip only

If you would rather not use an environment file:

```bash
conda create -n sampledisco python=3.10 pip && conda activate sampledisco
pip install "sampledisco[all]"
conda install -c bioconda -c conda-forge bedtools   # binary that scGLUE / pybedtools call
```

| Install | What you get |
|---|---|
| `pip install sampledisco` | The core pipeline: RNA, ATAC, multi-omics from a pre-integrated file, and all downstream analyses |
| `[glue]` extra | scGLUE 0.3.2 and harmony-pytorch |
| `[maxfuse]` extra | MaxFuse 0.0.2 |
| `[all]` extra | Both of the above |

This route uses PyPI's default torch wheel and does not pin numpy or scanpy to the validated versions. It is tested, but use an environment file when you want the exact validated stack.

!!! note "Core-only installs use a different Harmony implementation"
    Cell-level Harmony runs on **harmony-pytorch** whenever it is installed: in both environment files and with `[glue]` or `[all]`. This is the original code path. A core-only `pip install sampledisco` has no harmony-pytorch and falls back to **harmonypy**. Both are valid Harmony implementations, but they do not give identical cell embeddings, so cell types and downstream numbers differ between the two installs. Use one install type throughout a project.

## What the environment contains

| Feature | Packages | Tutorial |
|---|---|---|
| Core pipeline (RNA / ATAC / sample embedding / downstream) | scanpy 1.11.4, anndata 0.10.9, harmonypy 0.2.0, numpy 2.0.2 | [Quick start](tutorials/quickstart.md) |
| scGLUE integration (RNA + ATAC from scratch) | scglue 0.3.2, torch 2.5.1, bedtools 2.31, pyensembl | [Multi-omics](tutorials/multiomics.md) |
| MaxFuse (RNA + protein / flow) | maxfuse 0.0.2 | [RNA + flow cytometry](tutorials/multiomics_maxfuse.md) |
| GPU acceleration (GPU env only) | RAPIDS 24.12 (cuml, cudf, cugraph), rapids-singlecell 0.13.1 | [GPU acceleration](tutorials/gpu.md) |

??? note "Why these version pins"
    - **scGLUE** must come from PyPI, never `conda install -c bioconda scglue`: the bioconda recipe caps `numpy<1.22`.
    - scGLUE 0.3.2 needs **`anndata<0.11`** (it imports a class anndata 0.11 removed) and **`torch<2.7`** (its training loop uses an argument torch 2.7 removed). Its own package metadata declares neither cap.
    - **MaxFuse 0.0.2** needs **`anndata<0.13`**.
    - **`pynndescent<0.7`** (and so `umap-learn<0.6`): pynndescent 0.7 breaks `scanpy.pp.neighbors` on more than about 8,000 cells (`'csr_array' object has no attribute 'getnnz'`). sampledisco 0.6.0 declares this pin itself.

## Upgrading from 0.5.x

- **In an environment built from an env file:** run `pip install --no-deps sampledisco==0.6.0`.
- **In a pip-only environment:** run `pip install -U "sampledisco[all]"`.

The [changelog](https://github.com/J041120h/SampleDisco/blob/main/CHANGELOG.md) lists what changed. If you stay on 0.5.0 with a fresh pip install, add `pip install "pynndescent<0.7"`.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Illegal instruction (core dumped)` when importing scglue or pyensembl | The CPU lacks AVX2, which the conda `polars` build needs. Use a newer machine or node. |
| `libtorch_cuda.so: undefined symbol: ncclCommResume` when importing anndata or scanpy | A system CUDA directory on `LD_LIBRARY_PATH` shadows the env's own libraries. Run `unset LD_LIBRARY_PATH` (or remove the cluster CUDA module) before activating the env. |
| `cudaErrorInsufficientDriver` when importing cuml | The GPU driver is older than the RAPIDS build. RAPIDS 24.12 needs a CUDA ≥ 12.0 driver. |
| `'csr_array' object has no attribute 'getnnz'` in `sc.pp.neighbors` | pynndescent 0.7 is installed. Run `pip install "pynndescent<0.7"`. |
| `TypeError: 'module' object is not callable` after `from sampledisco import wrapper` | sampledisco 0.5.0. Upgrade, or use `from sampledisco.wrapper.wrapper import wrapper`. |
