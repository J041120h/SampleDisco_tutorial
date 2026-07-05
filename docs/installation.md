# Installation

`pip install sampledisco` gives you the full **CPU** pipeline — RNA, ATAC, and multi-omics from a pre-integrated file. GPU acceleration and *training GLUE from scratch* are optional add-ons you install yourself (they're environment-specific and can't be set up reliably via a pip extra).

!!! info "Requirements"
    - **Python ≥ 3.10** — on a fresh machine, create an env first: `conda create -n sampledisco python=3.10 && conda activate sampledisco`.
    - **OS:** macOS or Linux. The CPU core has **no PyTorch** (Harmony runs via harmonypy); GPU acceleration is **Linux + NVIDIA only**.
    - **Resources:** CPU demo ≈ 4–8 GiB RAM; the full RAPIDS GPU env is ≈ 16 GiB installed (allow ~30 GiB transient for the conda package cache during the solve) and needs ≥ 8 GiB RAM. RAISIN (`rna_cluster_dge`) is multithreaded and memory-heavy.

## 1. Core install (CPU)

```bash
pip install sampledisco
```

Verify:

```bash
python -c "import sampledisco; print(sampledisco.__version__)"
sampledisco --help
```

The supported interface is the config-driven CLI — `sampledisco -m complex --config config.yaml` (see the [Configuration guide](tutorials/configuration.md)).

## 2. GPU acceleration & training GLUE from scratch (optional)

These aren't bundled — they're CUDA-driver- and OS-specific. If any part is missing, SampleDisco falls back cleanly to CPU.

!!! note "GPU pays off only on large data"
    On the small demo (8 samples, ~30 k cells) the import and host↔device transfer overhead outweighs the speedup — a GPU run may be **no faster, or even slightly slower**, than CPU. The GPU paths matter on large datasets. Also note the GPU and CPU embeddings are **not equivalent, and can differ materially**: the composition k-means uses different backends (cuML vs scikit-learn), and on the demo the CPU-vs-GPU sample-distance-matrix correlation was only ≈0.4–0.6 — enough to change some downstream sample clusters. (The sample-level batch correction itself is the same on both — `harmonypy`.) **For reproducible or publication results, pick one backend and stay on it.**

**RAPIDS (GPU — Linux + NVIDIA).** No single command works everywhere; get the line for *your* driver from the [RAPIDS install selector](https://docs.rapids.ai/install) and [rapids-singlecell docs](https://rapids-singlecell.readthedocs.io/). The known-good stack we test (CUDA 12.5 driver):

```bash
conda install -c rapidsai -c conda-forge -c nvidia \
    cuml=24.12 cudf=24.12 cugraph=24.12 rmm=24.12 cuvs=24.12 cupy=13 cuda-version=12.5
pip install rapids-singlecell==0.13.1 --no-deps
pip install docrep scikit-image          # rapids-singlecell deps that --no-deps skips
```

For GPU-accelerated Harmony during **preprocessing** (not just the sample-embedding step), also install `harmony-pytorch` from the block below — without it, preprocessing Harmony silently runs on CPU (`harmonypy`).

**Training GLUE from scratch / torch-based Harmony.** Install `bedtools`, a **CUDA-12** PyTorch, scGLUE, and optionally harmony-pytorch. Use a cu12 torch build — the default PyPI torch may be a CUDA-13 wheel that silently runs on CPU on a CUDA-12 driver:

```bash
conda install -c bioconda bedtools                                           # scGLUE needs it (macOS: brew install bedtools)
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121   # match your driver's CUDA
pip install scglue==0.3.2 harmony-pytorch                                     # scGLUE + faster Harmony
```

See scGLUE's [install guide](https://scglue.readthedocs.io/en/latest/install.html) — install it from PyPI, never `conda install -c bioconda scglue` (that recipe caps `numpy<1.22`). scGLUE/harmony-pytorch also run on CPU (slower); without harmony-pytorch, Harmony uses harmonypy.

!!! warning "scGLUE 0.3.2 needs `anndata < 0.11`"
    `scglue==0.3.2` imports `anndata._core.sparse_dataset.SparseDataset`, which was **removed in anndata 0.11**. A plain `pip install sampledisco` pulls `anndata 0.11.x`, so `import scglue` then fails with `AttributeError: … has no attribute 'SparseDataset'`. This only affects *training GLUE from scratch* (the demo starts from the pre-integrated file and never imports scGLUE). If you do need from-scratch GLUE, train it in a **dedicated env** that pins `anndata<0.11`:

    ```bash
    conda create -n sampledisco-glue python=3.10
    conda activate sampledisco-glue
    pip install "anndata<0.11" scglue==0.3.2 harmony-pytorch
    pip install sampledisco --no-deps      # avoid re-pulling anndata 0.11
    ```

!!! warning "RAPIDS is driver-specific"
    The 24.12 pins target a CUDA 12.0–12.5 driver. RAPIDS 25.04+ needs ≥ CUDA 12.6 and fails to import (`cudaErrorInsufficientDriver`) on older drivers; bump the pins on newer drivers.

## Alternative: pinned conda environment (most reproducible)

The repo ships exact-version env files. Clone it, create the env, and install the package from the clone:

```bash
git clone https://github.com/J041120h/SampleDisco.git && cd SampleDisco
conda env create -f environment-cpu.yml      # CPU (macOS / Linux)
conda activate sampledisco-cpu
pip install -e . --no-deps
```

For GPU, use `environment-gpu.yml` (RAPIDS 24.12 + torch 2.5.1+cu121) and add `pip install rapids-singlecell==0.13.1 --no-deps` after activating.
