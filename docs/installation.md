# Installation

One conda environment runs everything on this site: the RNA, ATAC and multi-omics pipelines, scGLUE training (with the `bedtools` binary it needs), MaxFuse and every downstream analysis. Use the **CPU** environment, or the **GPU** environment, which adds NVIDIA RAPIDS. Both install SampleDisco itself. A pip-only route is also available.

**Requirements:** Linux x86-64 with an AVX2 CPU, and conda (Miniforge or Miniconda); the GPU environment also needs an NVIDIA GPU. On the demo data the RNA pipeline needs about 4 GB of RAM, RNA + ATAC about 8.5 GB, and the optional RAISIN DGE 16–18 GB.

## 1. CPU environment (recommended)

```bash
curl -fLO https://j041120h.github.io/SampleDisco_tutorial/assets/environment-cpu.yml
conda env create -f environment-cpu.yml
conda activate sampledisco-cpu
```

The build takes about 15–25 minutes. Check it:

```bash
python -c "import sampledisco, scglue, maxfuse; print('sampledisco', sampledisco.__version__)"
sampledisco --help
bedtools --version
```

Next: [download the demo data](tutorials/demo_data.md) and run the [Quick start](tutorials/quickstart.md).

## 2. GPU environment (Linux + NVIDIA)

```bash
curl -fLO https://j041120h.github.io/SampleDisco_tutorial/assets/environment-gpu.yml
conda env create -f environment-gpu.yml      # building on a machine without a GPU? prefix CONDA_OVERRIDE_CUDA=12.5
conda activate sampledisco-gpu
pip install --no-deps rapids-singlecell==0.13.1 sampledisco==0.10.0
python -c "import cuml, cupy, rapids_singlecell, torch; print('GPU OK', torch.cuda.is_available())"
```

- The build takes about 30 minutes. You get the CPU environment plus RAPIDS 24.12 and torch 2.5.1+cu121.
- **Keep `--no-deps`.** The env file already installs every dependency. A dependency-resolving `pip install sampledisco` would replace the `pyarrow` that RAPIDS needs and break `cudf`.
- RAPIDS 24.12 needs a driver for CUDA 12.0–12.5. With a driver for CUDA 12.6 or later you may move the RAPIDS pins to 25.x.
- Set `use_gpu: true` in the config to use the GPU; see [GPU acceleration](tutorials/gpu.md).

## 3. pip only

```bash
conda create -n sampledisco -c conda-forge -c bioconda python=3.10 pip bedtools   # bedtools: binary that scGLUE / pybedtools call
conda activate sampledisco
pip install "sampledisco[all]"
```

- Create the environment with Python, pip and bedtools in one step, as above. Adding bedtools to an existing Miniconda environment later can replace its Python with an incompatible build, and then nothing imports.
- `pip install sampledisco` alone gives the core pipeline: RNA, ATAC, multi-omics from a pre-integrated file, and every downstream analysis. The `[glue]` extra adds scGLUE, `[maxfuse]` adds MaxFuse, and `[all]` adds both.
- This route does not pin numpy, scanpy or torch to the validated versions. Use an environment file when you want the exact validated stack.

**Upgrading:** in an environment built from an env file, run `pip install --no-deps -U sampledisco`; in a pip-only environment, run `pip install -U "sampledisco[all]"`. The [changelog](https://github.com/J041120h/SampleDisco/blob/main/CHANGELOG.md) lists what changed.

??? note "What the environments contain, and why the pins"
    | Feature | Packages |
    |---|---|
    | Core pipeline | scanpy 1.11.4, anndata 0.10.9, numpy 2.0.2, harmonypy 0.2.0, harmony-pytorch 0.1.8 |
    | scGLUE (RNA + ATAC from scratch) | scglue 0.3.2, torch 2.5.1, bedtools 2.31, pyensembl |
    | MaxFuse (RNA + protein / flow) | maxfuse 0.0.2 |
    | GPU (GPU environment only) | RAPIDS 24.12 (cuml, cudf, cugraph), rapids-singlecell 0.13.1, torch 2.5.1+cu121 |

    - scGLUE must come from PyPI, never `conda install -c bioconda scglue`: the bioconda recipe caps `numpy<1.22`.
    - scGLUE 0.3.2 needs `anndata<0.11` and `torch<2.7`, although its own package metadata declares neither.
    - MaxFuse 0.0.2 needs `anndata<0.13`.
    - `pynndescent<0.7` (and so `umap-learn<0.6`): pynndescent 0.7 breaks `scanpy.pp.neighbors` on more than about 8,000 cells. SampleDisco declares this pin itself.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Illegal instruction (core dumped)` when importing scglue or pyensembl | The CPU lacks AVX2, which the conda `polars` build needs. Use a newer machine or node. |
| `libtorch_cuda.so: undefined symbol: ncclCommResume` when importing anndata or scanpy | A system CUDA directory on `LD_LIBRARY_PATH` shadows the env's own libraries. Run `unset LD_LIBRARY_PATH` (or unload the cluster CUDA module) before activating the env. |
| `cudaErrorInsufficientDriver` when importing cuml | The GPU driver is older than RAPIDS 24.12 needs (a CUDA 12.0 driver or later). |
| `'csr_array' object has no attribute 'getnnz'` in `sc.pp.neighbors` | pynndescent 0.7 is installed. Run `pip install "pynndescent<0.7"`. |
