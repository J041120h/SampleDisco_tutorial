# Demo data

Every tutorial on this site runs on the same public datasets, published on Zenodo under DOI **[10.5281/zenodo.22549610](https://doi.org/10.5281/zenodo.22549610)** (**CC BY-NC 4.0 — noncommercial use**) — this DOI always resolves to the **latest version**, so it's the link to share. The commands below fetch v2.0 (Zenodo record 24451067), which is version-pinned so the MD5s below always match.

## What's in the box

Five single-cell AnnData (`.h5ad`) files from two cohorts: unpaired scRNA-seq and scATAC-seq from COVID-19 and healthy PBMC donors (plus their pre-computed scGLUE-integrated object), and unpaired scRNA-seq and spectral flow cytometry from the Sound Life cohort of young and older adults.

| File | Contents | Shape | Units | Size | MD5 |
| --- | --- | --- | --- | --- | --- |
| `test_RNA.h5ad` | scRNA-seq counts | 29,989 cells × 14,881 genes | 8 samples (CoV / HD donors) | 94 MB | `69f00904fb9ff21e1d8e83afd63aa441` |
| `test_ATAC.h5ad` | scATAC-seq counts | 29,021 cells × 230,356 peaks | 8 samples (SRA accessions) | 260 MB | `9053f14b396ae2553677d00b69ba18a8` |
| `test_multiomics_integrated.h5ad` | scGLUE-integrated RNA + ATAC | 59,010 cells × 13,820 genes | 29,989 RNA + 29,021 ATAC cells | 54 MB | `3d474460e651fe61cc0751a3911512d1` |
| `soundlife_RNA.h5ad` | Sound Life scRNA-seq counts | 32,000 cells × 18,456 genes | 32 samples (16 young / 16 older) | 137 MB | `0446c8728d187cb95f4e55cc5c6f1339` |
| `soundlife_FLOW.h5ad` | Sound Life spectral flow (PS1) | 32,000 events × 24 markers | the same 32 samples | 7 MB | `b57f57b32d39161ac8081a5f1287f33e` |

The third file is **optional** — the [multi-omics tutorial](multiomics.md) starts from it to skip scGLUE training. It already carries the GLUE joint embedding, stored under the pre-0.3 key names: the sample-**preserved** view `obsm['X_glue']` (with a copy under the legacy name `Z_cmd`, today's `Z_rmd`) and the sample-**removed** view `obsm['Z_clust']` (today's `Z_comp`). SampleDisco resolves these names automatically; `migrate_obsm_keys(adata)` from `sampledisco.utils.embedding_keys` copies them onto the current names in memory (see [legacy key names](../legacy_keys.md)). The file also holds the joint cell-type labels (`obs['cell_type']`) and a sample embedding (`uns['X_DR_sample']`) from the original run — the tutorials recompute both.

!!! tip "No metadata CSV needed"
    All per-cell metadata is already inside each file's `.obs`: the sample identifier (`sample`), a severity phenotype (`sev.level`; `1` = healthy, `4` = severe COVID-19), and — for RNA — a `batch` column (`Aruna` vs `Wilk`) and pre-computed `celltype` labels. You do **not** need a separate `sample_meta.csv` or HVG list.

The RNA and ATAC modalities are **unpaired** (they share no donors), which is exactly the setting the [multi-omics pipeline](multiomics.md) is built for: GLUE anchors RNA and ATAC into a joint embedding without any cross-modality cell pairing.

The two **Sound Life** files are used by the [RNA + flow cytometry (MaxFuse) tutorial](multiomics_maxfuse.md). They hold one baseline visit (Flu Year 1 Day 0) per donor — 16 young and 16 older adults, four per sex × CMV group — with 1,000 cells per sample in each modality. Here the modalities share **samples** but not cells. `.obs` carries `sample`, `subject`, `age_group`, `age_at_draw`, `sex`, `cmv` and the processing `batch`; the RNA file adds AIFI cell types (`cell_type`, `cell_type_l1`) and the flow file the marker-gated `lineage`. The flow `.X` is `arcsinh(unmixed / 3000)`, with the unmixed values in `layers['unmixed']`.

!!! warning "Noncommercial use"
    The Sound Life data come from the Allen Institute for Immunology and are redistributed under the [Allen Institute Terms of Use](https://alleninstitute.org/terms-of-use/), which permit research and other noncommercial use. From v2.0 the whole record is licensed CC BY-NC 4.0. If you use the Sound Life files, cite the source publication (below).

## Download

Grab the files into a local `data/` folder:

=== "wget"

    ```bash
    mkdir -p data
    wget -O data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
    wget -O data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
    # optional — used by the multi-omics tutorial:
    wget -O data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
    # used by the RNA + flow cytometry (MaxFuse) tutorial:
    wget -O data/soundlife_RNA.h5ad  "https://zenodo.org/records/24451067/files/soundlife_RNA.h5ad?download=1"
    wget -O data/soundlife_FLOW.h5ad "https://zenodo.org/records/24451067/files/soundlife_FLOW.h5ad?download=1"
    ```

=== "curl"

    ```bash
    mkdir -p data
    curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
    curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
    # optional — used by the multi-omics tutorial:
    curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
    # used by the RNA + flow cytometry (MaxFuse) tutorial:
    curl -fL -o data/soundlife_RNA.h5ad  "https://zenodo.org/records/24451067/files/soundlife_RNA.h5ad?download=1"
    curl -fL -o data/soundlife_FLOW.h5ad "https://zenodo.org/records/24451067/files/soundlife_FLOW.h5ad?download=1"
    ```

=== "Python"

    ```python
    import urllib.request, os
    os.makedirs("data", exist_ok=True)
    base = "https://zenodo.org/records/24451067/files"
    files = ["test_RNA.h5ad", "test_ATAC.h5ad", "test_multiomics_integrated.h5ad",
             "soundlife_RNA.h5ad", "soundlife_FLOW.h5ad"]
    for f in files:
        urllib.request.urlretrieve(f"{base}/{f}?download=1", f"data/{f}")
    ```

## Verify your download (optional)

A corrupted or partial download can surface later as a cryptic HDF5 error (e.g. "bad symbol table node signature"). To catch it early, check the MD5s against the [table above](#whats-in-the-box):

=== "Linux"

    ```bash
    md5sum data/*.h5ad
    ```

=== "macOS"

    ```bash
    md5 data/*.h5ad
    ```

Expected (tool-independent `<file>: <hash>` pairs — `md5sum` prepends the `data/` path and prints `hash  file`, while macOS `md5` prints `MD5 (file) = hash`, so match on the hash regardless of layout):

```
test_RNA.h5ad:                     69f00904fb9ff21e1d8e83afd63aa441
test_ATAC.h5ad:                    9053f14b396ae2553677d00b69ba18a8
test_multiomics_integrated.h5ad:   3d474460e651fe61cc0751a3911512d1
soundlife_RNA.h5ad:                0446c8728d187cb95f4e55cc5c6f1339
soundlife_FLOW.h5ad:               b57f57b32d39161ac8081a5f1287f33e
```

## Citation

> Jiang, H., & Ji, H. *SampleDisco tutorial demo data: COVID-19 scRNA-seq/scATAC-seq and Sound Life scRNA-seq/flow cytometry subsets.* Zenodo. https://doi.org/10.5281/zenodo.22549610

For the Sound Life files, also cite the source:

> Gong, Q. et al. Longitudinal multi-omic immune profiling reveals age-related immune cell dynamics in healthy adults. *Nature* (2025). https://doi.org/10.1038/s41586-025-09686-5 — data: Allen Institute for Immunology, Human Immune System Explorer, https://doi.org/10.57785/bn17-3h16
