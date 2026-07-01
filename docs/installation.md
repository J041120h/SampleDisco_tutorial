# Installation

Install SampleDisco from PyPI with `pip`. It runs on CPU by default; GPU acceleration is optional and turns on automatically when the NVIDIA RAPIDS stack is importable — no separate build and no reinstall — otherwise it falls back to CPU.

!!! info "Requirements"
    - **Python ≥ 3.10** (developed and tested on 3.10). On a fresh machine or HPC login node, create a matching environment first — e.g. `conda create -n sampledisco python=3.10 && conda activate sampledisco`.
    - **OS:** macOS or Linux. GPU acceleration is **Linux + NVIDIA only** (RAPIDS needs a CUDA GPU); on macOS everything runs on CPU.
    - **Core dependencies** — installed automatically by `pip`: NumPy ≥ 1.23, SciPy ≥ 1.9, scikit-learn ≥ 1.2, scanpy ≥ 1.11, anndata ≥ 0.10, harmonypy, leidenalg ≥ 0.9, numba ≥ 0.57 (and others). **No PyTorch** — Harmony runs on CPU via harmonypy.
    - **Multi-omics extra** — `pip install sampledisco[multiomics]` adds PyTorch ≥ 2.0 + scGLUE ≥ 0.3 (+ harmony-pytorch). Needed for scGLUE integration and for the GPU/torch-based Harmony path.
    - **`bedtools`** (system binary, [step 2](#2-optional-bedtools-only-for-scglue-training)) — needed **only** to train scGLUE from scratch. The multi-omics demo skips it via the [pre-integrated file](tutorials/multiomics.md#1-load-the-integrated-data).
    - **Resources** — the CPU demo runs comfortably in **≈4–8 GiB RAM**; the RAPIDS GPU install needs **≥8 GiB RAM and ~10 GiB disk** to unpack the CUDA libraries. (`rna_cluster_dge` / RAISIN is multithreaded and memory-heavy — see the demo config note.)

## Install SampleDisco

### 1. Core install (CPU) — from PyPI

```bash
pip install sampledisco                 # core: CPU, single-omics (no PyTorch)
pip install "sampledisco[multiomics]"   # adds PyTorch + scGLUE for multi-omics
```

The core install has **no PyTorch** and runs Harmony via harmonypy — enough for the full RNA-only and ATAC-only pipelines. Add the `multiomics` extra only if you need scGLUE integration (it also switches Harmony to the faster torch backend automatically).

### 2. Optional — bedtools (only for scGLUE training)

The `bedtools` binary is required **only if you train multi-omics scGLUE integration from scratch** (scGLUE / `pybedtools` call it). RNA-only, ATAC-only, and the multi-omics demo that resumes from the [pre-integrated file](tutorials/multiomics.md#1-load-the-integrated-data) do not need it.

```bash
conda install -c bioconda bedtools     # macOS: also available via `brew install bedtools`
```

### 3. GPU acceleration (optional — Linux + NVIDIA)

The GPU code paths (RAPIDS-accelerated normalization, Harmony, k-means / PCA, Leiden, scGLUE training) activate **only when the RAPIDS stack is importable**. RAPIDS requires an NVIDIA CUDA GPU and is **Linux-only** (skip this section on macOS). It is CUDA-driver-specific and conda-only, so you install it separately, matching your driver. The pins below target a CUDA 12.5 driver:

```bash
conda install -c rapidsai -c conda-forge -c nvidia \
    cuml=24.12 cudf=24.12 cugraph=24.12 rmm=24.12 cuvs=24.12 cupy=13 cuda-version=12.5
pip install rapids-singlecell==0.13.1 --no-deps
pip install docrep scikit-image        # rapids-singlecell runtime deps that --no-deps skips
```

For the multi-omics / torch-based Harmony GPU path, also install the extra with a **CUDA-12 torch build** — the default PyPI torch may be a CUDA-13 wheel that reports `torch.cuda.is_available() == False` on a CUDA-12 driver and silently runs on CPU:

```bash
pip install "sampledisco[multiomics]"
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121   # match your driver's CUDA
```

Then set `use_gpu: true` in your config. **You do not reinstall SampleDisco** — once those libraries are importable the GPU paths turn on automatically; if any part of the GPU stack is missing or the driver is too old, SampleDisco falls back cleanly to the CPU implementations (harmonypy Harmony, scikit-learn k-means).

!!! warning "RAPIDS version is driver-specific"
    RAPIDS is pinned to **24.12** on purpose: it is built for CUDA 12.0–12.5 and runs on a 12.5 driver. RAPIDS 25.04+ needs a newer driver (≥ CUDA 12.6) and fails at import with `cudaErrorInsufficientDriver` on a 12.5 driver. On machines with driver ≥ CUDA 12.6 you may bump the pins to 25.x.

## Alternative: reproducible conda environment

For an exactly-pinned stack (the most reproducible option), use the conda environment files that ship in the repository. **Clone the repo first** — that's where the environment files and the package source live:

```bash
git clone https://github.com/J041120h/SampleDisco.git
cd SampleDisco
```

Then create the environment from the matching file and install the package itself from the clone:

=== "CPU (macOS / Linux)"

    ```bash
    conda env create -f environment-cpu.yml
    conda activate sampledisco-cpu
    pip install -e . --no-deps        # the env already provides every dependency
    ```

=== "GPU (Linux + NVIDIA, RAPIDS 24.12)"

    ```bash
    conda env create -f environment-gpu.yml
    conda activate sampledisco-gpu
    pip install rapids-singlecell==0.13.1 --no-deps   # pip-only; --no-deps preserves the pinned RAPIDS
    pip install -e . --no-deps
    ```

This pins the whole stack (NumPy 2.0.2, scanpy 1.11.4, anndata 0.10.9, scGLUE 0.3.2; the GPU env adds RAPIDS 24.12 + torch 2.5.1+cu121). With the environment active, run the pipeline from a YAML config — see the [Configuration guide](tutorials/configuration.md), e.g. on the [demo data](tutorials/demo_data.md):

```bash
sampledisco -m complex --config your_config.yaml
```

!!! warning "scGLUE must come from PyPI"
    Install scGLUE with `pip install scglue==0.3.2`, never `conda install -c bioconda scglue` — the bioconda recipe caps `numpy<1.22` and breaks the rest of the stack. The env files already do this.

## Verify the install

Confirm the package imports and the console script is on your `PATH`:

```bash
python -c "import sampledisco; print(sampledisco.__version__)"
sampledisco --help
```

The primary interface is the config-driven wrapper, invoked through the `sampledisco` console script:

```bash
sampledisco -m complex --config config.yaml
# equivalently:
python -m sampledisco.cli -m complex --config config.yaml
```
