# Demo data

Every tutorial on this site is built from the same public datasets. They are hosted on Zenodo so you can download them and reproduce the entire walkthrough — RNA, ATAC, and multi-omics — on your own machine.

!!! abstract "Zenodo record"
    **SampleDisco tutorial demo data: scRNA-seq and scATAC-seq COVID-19 PBMC subsets**<br/>
    DOI (this version, 1.1): [10.5281/zenodo.20990399](https://doi.org/10.5281/zenodo.20990399) · License: **CC-BY-4.0** · Open access<br/>
    Always-latest (concept) DOI: [10.5281/zenodo.20988712](https://doi.org/10.5281/zenodo.20988712)

## What's in the box

Three single-cell AnnData (`.h5ad`) files — unpaired scRNA-seq and scATAC-seq from COVID-19 and healthy PBMC donors, plus their pre-computed scGLUE-integrated object.

| File | Contents | Shape | Units | Size | MD5 |
| --- | --- | --- | --- | --- | --- |
| `test_RNA.h5ad` | scRNA-seq counts | 29,989 cells × 14,881 genes | 8 samples (CoV / HD donors) | 501 MB | `d2f7f9a84a8082d09c18d102669966e3` |
| `test_ATAC.h5ad` | scATAC-seq counts | 29,021 cells × 230,356 peaks | 8 samples (SRA accessions) | 2.7 GB | `70df48c7f3007651ac23c2644e49e829` |
| `test_multiomics_integrated.h5ad` | scGLUE-integrated RNA + ATAC | 59,010 cells × 13,820 genes | 29,989 RNA + 29,021 ATAC cells | 54 MB | `3d474460e651fe61cc0751a3911512d1` |

The third file is **optional** — you only need it for the [multi-omics tutorial](multiomics.md) if you want to **skip scGLUE training** (see [below](#skip-scglue-training-multi-omics)). It already carries the GLUE joint embedding (`obsm['X_glue']`), the sample-removed view (`obsm['Z_clust']`), and joint cell-type labels.

!!! tip "No metadata CSV needed"
    All per-cell metadata is already inside each file's `.obs`: the sample identifier (`sample`), a severity phenotype (`sev.level`; `1` = healthy, `4` = severe COVID-19) used throughout the trajectory/association steps, and — for RNA — a `batch` column (`Aruna` vs `Wilk`) and pre-computed `celltype` labels. Because of this, every `*_sample_meta_path` is left empty in the demo config; you do **not** need a separate `sample_meta.csv` or HVG list to run the tutorials.

The RNA and ATAC modalities are **unpaired** (they share no donors), which is exactly the setting the [multi-omics pipeline](multiomics.md) is built for: GLUE anchors RNA and ATAC into a joint embedding without any cross-modality cell pairing.

## Download

Grab the files into a local `data/` folder:

=== "wget"

    ```bash
    mkdir -p data
    wget -O data/test_RNA.h5ad  "https://zenodo.org/records/20990399/files/test_RNA.h5ad?download=1"
    wget -O data/test_ATAC.h5ad "https://zenodo.org/records/20990399/files/test_ATAC.h5ad?download=1"
    # optional — only for the "skip scGLUE training" path:
    wget -O data/test_multiomics_integrated.h5ad "https://zenodo.org/records/20990399/files/test_multiomics_integrated.h5ad?download=1"
    ```

=== "curl"

    ```bash
    mkdir -p data
    curl -L -o data/test_RNA.h5ad  "https://zenodo.org/records/20990399/files/test_RNA.h5ad?download=1"
    curl -L -o data/test_ATAC.h5ad "https://zenodo.org/records/20990399/files/test_ATAC.h5ad?download=1"
    # optional — only for the "skip scGLUE training" path:
    curl -L -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/20990399/files/test_multiomics_integrated.h5ad?download=1"
    ```

=== "Python"

    ```python
    import urllib.request, os
    os.makedirs("data", exist_ok=True)
    base = "https://zenodo.org/records/20990399/files"
    files = ["test_RNA.h5ad", "test_ATAC.h5ad", "test_multiomics_integrated.h5ad"]
    for f in files:
        urllib.request.urlretrieve(f"{base}/{f}?download=1", f"data/{f}")
    ```

Verify the downloads (optional):

```bash
md5sum data/*.h5ad          # macOS: use `md5 data/*.h5ad`
# d2f7f9a84a8082d09c18d102669966e3  data/test_RNA.h5ad
# 70df48c7f3007651ac23c2644e49e829  data/test_ATAC.h5ad
# 3d474460e651fe61cc0751a3911512d1  data/test_multiomics_integrated.h5ad
```

## Run the demo

Once the package is [installed](../installation.md), there are two ways to reproduce the tutorials.

### Option A — one command with the demo config

Download the ready-to-run config — [**`config_demo.yaml`**](../assets/config_demo.yaml) — into the **same folder** that contains your `data/` directory, then:

```bash
sampledisco -m complex --config config_demo.yaml
```

The config already points at `./data/test_RNA.h5ad` and `./data/test_ATAC.h5ad` (relative paths, so run it from that folder), writes everything under `./sampledisco_demo_output/{rna,atac,multiomics}/`, and runs on CPU by default (`use_gpu: false`).

- **RNA and ATAC** pipelines are enabled out of the box and run end to end (preprocess → embedding → all downstream).
- **Multi-omics** is off by default (`run_multiomics_pipeline: false`); flip it to `true` to add it. The config is pre-wired to **skip scGLUE training** by loading the integrated file (see below), so the multi-omics branch stays light and needs no `bedtools`.

!!! warning "Use the full config — it is validated exactly"
    In `complex` mode the CLI checks the YAML **key-for-key** against `wrapper()`: every parameter must be present and no extras are allowed. `config_demo.yaml` is a complete, validated config — edit values, but don't delete keys. See the [Configuration guide](configuration.md) for what each block means.

### Skip scGLUE training (multi-omics)

scGLUE integration (training the cross-modality model + downloading genome annotation) is the slowest part of the multi-omics pipeline and needs the `bedtools` binary. To skip it, use the pre-computed `test_multiomics_integrated.h5ad` and resume from it:

```yaml
run_multiomics_pipeline: true
multiomics_integration: false                                         # do NOT train scGLUE
multiomics_integrated_h5ad_path: "./data/test_multiomics_integrated.h5ad"   # load the integrated object instead
```

These three keys are already set this way in [`config_demo.yaml`](../assets/config_demo.yaml). With them, the multi-omics branch loads the joint embedding directly and proceeds to cell typing → sample embedding → downstream — no GLUE training, no `bedtools`. To instead train scGLUE yourself from the raw counts, set `multiomics_integration: true` and `multiomics_integrated_h5ad_path: null` (and install `bedtools`).

### Option B — step through the pipeline tutorials

Prefer to see each function and its output? The pipeline tutorials call the building blocks directly on the downloaded files:

- [RNA pipeline](rna.md) — uses `data/test_RNA.h5ad`
- [ATAC pipeline](atac.md) — uses `data/test_ATAC.h5ad`
- [Multi-omics pipeline](multiomics.md) — uses **both** raw files (train GLUE), or `test_multiomics_integrated.h5ad` to skip it
- [Downstream analysis](downstream/index.md) — runs off the resulting sample embedding

## Citation

If you use these data, please cite the Zenodo record:

> Jiang, H., & Ji, H. *SampleDisco tutorial demo data: scRNA-seq and scATAC-seq COVID-19 PBMC subsets.* Zenodo. https://doi.org/10.5281/zenodo.20988712
