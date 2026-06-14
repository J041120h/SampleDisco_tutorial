# Installation

SampleDisco is **one package** (`sampledisco`). The core install is pip-only and runs on CPU; **GPU acceleration lights up automatically when the NVIDIA RAPIDS stack is present** in your environment — there is no separate "GPU build" and you never reinstall SampleDisco.

## 1. Core install (CPU) — from PyPI

```bash
pip install sampledisco
```

!!! tip "Phylogenetic tree clustering (optional)"
    For phylogenetic tree-based sample clustering, install the optional extra:

    ```bash
    pip install 'sampledisco[trees]'
    ```

## 2. System prerequisite — bedtools

scGLUE (the multi-omics integrator) and `pybedtools` call the `bedtools` binary, which pip cannot provide. Install it via conda:

```bash
conda install -c bioconda bedtools
```

## 3. GPU acceleration (optional — install RAPIDS yourself)

The GPU code paths (RAPIDS-accelerated normalization, Harmony, k-means / PCA, Leiden, scGLUE training) activate **only when the RAPIDS stack is importable**. RAPIDS is CUDA-driver-specific and conda-only, so you install it separately, matching your driver. The pins below target a CUDA-12.5 driver (e.g. the cluster's GPU nodes):

```bash
conda install -c rapidsai -c conda-forge -c nvidia \
    cuml=24.12 cudf=24.12 cugraph=24.12 rmm=24.12 cuvs=24.12 cupy=13 cuda-version=12.5
pip install rapids-singlecell==0.13.1 --no-deps
```

Then set `use_gpu: true` in your config. **You do not reinstall SampleDisco** — once those libraries are importable the GPU paths turn on automatically; if they are missing or the driver is too old, SampleDisco falls back cleanly to CPU equivalents (harmonypy, scikit-learn k-means, PyTorch CPU).

!!! warning "RAPIDS version is driver-specific"
    RAPIDS is pinned to **24.12** on purpose: it is built for CUDA 12.0–12.5 and runs on the 12.5 driver. RAPIDS 25.04+ needs a newer driver (≥ CUDA 12.6) and fails at import with `cudaErrorInsufficientDriver` on 12.5 nodes. On nodes with driver ≥ CUDA 12.6 you may bump the pins to 25.x.

## One-command environments (recommended)

For a fully reproducible environment — including bedtools, scGLUE 0.3.2, and (GPU) the full RAPIDS 24.12 stack — use the provided conda files instead of the manual steps above:

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

Both environments were built and smoke-tested end-to-end on JHPCE (numpy 2.0.2 + scGLUE 0.3.2 + scanpy 1.11.4 + anndata 0.10.9; the GPU env adds RAPIDS 24.12 + torch 2.5.1+cu121, all coexisting on an A100).

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
