# Installation

Install SampleDisco from PyPI with `pip`. It runs on CPU by default; GPU acceleration is optional and turns on automatically when the NVIDIA RAPIDS stack is importable — no separate build and no reinstall — otherwise it falls back to CPU.

!!! info "Requirements"
    - **Python ≥ 3.10** (developed and tested on 3.10).
    - **OS:** macOS or Linux. GPU acceleration is **Linux + NVIDIA only** (RAPIDS needs a CUDA GPU); on macOS everything runs on CPU.
    - **Core dependencies** — installed automatically by `pip`: NumPy ≥ 1.23, SciPy ≥ 1.9, scikit-learn ≥ 1.2, scanpy ≥ 1.11, anndata ≥ 0.10, PyTorch ≥ 2.0, scGLUE ≥ 0.3, leidenalg ≥ 0.9, numba ≥ 0.57 (and others).
    - **`bedtools`** (system binary, [step 2](#2-optional-bedtools-only-for-scglue-training)) — needed **only** to train scGLUE from scratch. The multi-omics demo skips it via the [pre-integrated file](tutorials/multiomics.md#1-load-the-integrated-data).

## Install SampleDisco

### 1. Core install (CPU) — from PyPI

```bash
pip install sampledisco
```

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
```

Then set `use_gpu: true` in your config. **You do not reinstall SampleDisco** — once those libraries are importable the GPU paths turn on automatically; if they are missing or the driver is too old, SampleDisco falls back cleanly to CPU equivalents (harmonypy, scikit-learn k-means, PyTorch CPU).

!!! warning "RAPIDS version is driver-specific"
    RAPIDS is pinned to **24.12** on purpose: it is built for CUDA 12.0–12.5 and runs on a 12.5 driver. RAPIDS 25.04+ needs a newer driver (≥ CUDA 12.6) and fails at import with `cudaErrorInsufficientDriver` on a 12.5 driver. On machines with driver ≥ CUDA 12.6 you may bump the pins to 25.x.

## Alternative: one-command environment

For a fully reproducible environment — including bedtools, scGLUE 0.3.2, and (GPU) the full RAPIDS 24.12 stack — use the provided conda files instead of the steps above. This is the **most reproducible option**:

```bash
# CPU
conda env create -f environment-cpu.yml
conda activate sampledisco-cpu
pip install sampledisco          # or: pip install -e . --no-deps  from a clone

# GPU (RAPIDS 24.12)
conda env create -f environment-gpu.yml
conda activate sampledisco-gpu
pip install rapids-singlecell==0.13.1 --no-deps   # rapids-singlecell is pip-only
pip install sampledisco          # or: pip install -e . --no-deps  from a clone
```

The CPU stack is cross-platform (macOS + Linux); the GPU environment is Linux + NVIDIA only. Pinned versions: numpy 2.0.2 + scGLUE 0.3.2 + scanpy 1.11.4 + anndata 0.10.9; the GPU environment adds RAPIDS 24.12 + torch 2.5.1+cu121.

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
