# Demo data

Every tutorial on this site runs on the same public datasets, published on Zenodo — **[zenodo.org/records/21019419](https://zenodo.org/records/21019419)** (CC-BY-4.0). Download them below, then follow the pipeline tutorials.

## What's in the box

Three single-cell AnnData (`.h5ad`) files — unpaired scRNA-seq and scATAC-seq from COVID-19 and healthy PBMC donors, plus their pre-computed scGLUE-integrated object.

| File | Contents | Shape | Units | Size | MD5 |
| --- | --- | --- | --- | --- | --- |
| `test_RNA.h5ad` | scRNA-seq counts | 29,989 cells × 14,881 genes | 8 samples (CoV / HD donors) | 94 MB | `69f00904fb9ff21e1d8e83afd63aa441` |
| `test_ATAC.h5ad` | scATAC-seq counts | 29,021 cells × 230,356 peaks | 8 samples (SRA accessions) | 260 MB | `9053f14b396ae2553677d00b69ba18a8` |
| `test_multiomics_integrated.h5ad` | scGLUE-integrated RNA + ATAC | 59,010 cells × 13,820 genes | 29,989 RNA + 29,021 ATAC cells | 54 MB | `3d474460e651fe61cc0751a3911512d1` |

The third file is **optional** — the [multi-omics tutorial](multiomics.md) starts from it to skip scGLUE training. It already carries the GLUE joint embedding (`obsm['X_glue']`), the sample-removed view (`obsm['Z_clust']`), and joint cell-type labels.

!!! tip "No metadata CSV needed"
    All per-cell metadata is already inside each file's `.obs`: the sample identifier (`sample`), a severity phenotype (`sev.level`; `1` = healthy, `4` = severe COVID-19), and — for RNA — a `batch` column (`Aruna` vs `Wilk`) and pre-computed `celltype` labels. You do **not** need a separate `sample_meta.csv` or HVG list.

The RNA and ATAC modalities are **unpaired** (they share no donors), which is exactly the setting the [multi-omics pipeline](multiomics.md) is built for: GLUE anchors RNA and ATAC into a joint embedding without any cross-modality cell pairing.

## Download

Grab the files into a local `data/` folder:

=== "wget"

    ```bash
    mkdir -p data
    wget -O data/test_RNA.h5ad  "https://zenodo.org/records/21019419/files/test_RNA.h5ad?download=1"
    wget -O data/test_ATAC.h5ad "https://zenodo.org/records/21019419/files/test_ATAC.h5ad?download=1"
    # optional — used by the multi-omics tutorial:
    wget -O data/test_multiomics_integrated.h5ad "https://zenodo.org/records/21019419/files/test_multiomics_integrated.h5ad?download=1"
    ```

=== "curl"

    ```bash
    mkdir -p data
    curl -L -o data/test_RNA.h5ad  "https://zenodo.org/records/21019419/files/test_RNA.h5ad?download=1"
    curl -L -o data/test_ATAC.h5ad "https://zenodo.org/records/21019419/files/test_ATAC.h5ad?download=1"
    # optional — used by the multi-omics tutorial:
    curl -L -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/21019419/files/test_multiomics_integrated.h5ad?download=1"
    ```

=== "Python"

    ```python
    import urllib.request, os
    os.makedirs("data", exist_ok=True)
    base = "https://zenodo.org/records/21019419/files"
    files = ["test_RNA.h5ad", "test_ATAC.h5ad", "test_multiomics_integrated.h5ad"]
    for f in files:
        urllib.request.urlretrieve(f"{base}/{f}?download=1", f"data/{f}")
    ```

## Citation

> Jiang, H., & Ji, H. *SampleDisco tutorial demo data: scRNA-seq and scATAC-seq COVID-19 PBMC subsets.* Zenodo. https://doi.org/10.5281/zenodo.20988712
