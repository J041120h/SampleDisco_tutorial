# Installation

`pip install sampledisco` installs the full **CPU** pipeline (RNA, ATAC, and multi-omics from a pre-integrated file). GPU acceleration and training GLUE from scratch are optional, environment-specific add-ons.

!!! info "Requirements"
    - **Python ≥ 3.10**, **macOS or Linux**. GPU acceleration is **Linux + NVIDIA only**.
    - **RAM:** ≈ 8 GiB for RNA; **≈ 16 GiB** for the full RNA+ATAC demo (ATAC's peak matrix is large).

## 1. Core install (CPU)

```bash
conda create -n sampledisco python=3.10 && conda activate sampledisco
pip install sampledisco
```

Verify, then drive the pipeline from a YAML config (see the [Configuration guide](tutorials/configuration.md)):

```bash
python -c "import sampledisco; print(sampledisco.__version__)"
sampledisco --help
```

## 2. GPU acceleration (optional — Linux + NVIDIA)

Build the GPU environment from the repo's validated env file. (Layering `conda install` of RAPIDS onto the CPU env does **not** solve — a channel clash.)

```bash
git clone https://github.com/J041120h/SampleDisco.git && cd SampleDisco
conda env create -f environment-gpu.yml         # RAPIDS 24.12 + torch 2.5.1+cu121
conda activate sampledisco-gpu
pip install rapids-singlecell==0.13.1 --no-deps
pip install sampledisco --no-deps
```

Verify the stack, then set `use_gpu: true` in your config:

```bash
python -c "import cuml, cupy, rapids_singlecell; print('GPU OK')"
```

If the stack is missing, runs fall back to CPU (the backend used is recorded as `gpu_available` in `sys_log/main_process_status.json`). GPU and CPU embeddings are **not identical** — the k-means backends differ — so stay on one backend for reproducible results; on small data the GPU may be no faster than CPU. The 24.12 pins target a CUDA 12.0–12.5 driver; bump to RAPIDS 25.x on a newer driver.

## 3. Training GLUE from scratch (optional)

Only needed to build the multi-omics integration yourself (the demo starts from a pre-integrated file). scGLUE 0.3.2 needs `anndata < 0.11`, so use a dedicated env:

```bash
conda create -n sampledisco-glue python=3.10 && conda activate sampledisco-glue
conda install -c bioconda bedtools
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121   # cu12 build first
pip install "anndata<0.11" scglue==0.3.2 harmony-pytorch
pip install sampledisco --no-deps
```

Install scGLUE from PyPI, never `conda install -c bioconda scglue` (that recipe caps `numpy<1.22`).

## Reproducible env files

The repo ships exact-version `environment-cpu.yml` and `environment-gpu.yml`. Create with `conda env create -f <file>`, then `pip install -e . --no-deps`.
