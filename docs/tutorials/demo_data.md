# Download the demo data

Every tutorial runs on two public datasets in one Zenodo record, [10.5281/zenodo.22549610](https://doi.org/10.5281/zenodo.22549610). That DOI always resolves to the latest version. The commands below fetch version 2.0 (record 24451067), so the checksums below always match.

The **COVID-19 PBMC** files hold unpaired scRNA-seq and scATAC-seq, each with 8 samples (4 healthy and 4 severe COVID-19, `sev.level` 1 and 4) from different donors; they feed the quick start and the RNA, ATAC, RNA + ATAC and downstream pages. The third file holds the same RNA and ATAC cells already integrated by scGLUE, so the [RNA + ATAC page](multiomics.md) can skip training. The two **Sound Life** files, used only by the [RNA + flow cytometry page](multiomics_maxfuse.md), hold scRNA-seq and spectral flow cytometry from one baseline visit of 32 healthy adults (16 young, 16 older), with 1,000 cells per sample in each modality; the two modalities share samples but not cells. All sample metadata is already in each file's `.obs`, so you need no separate metadata file.

| File | Contents | Size | MD5 |
| --- | --- | --- | --- |
| `test_RNA.h5ad` | COVID scRNA-seq counts: 29,989 cells × 14,881 genes | 94 MB | `69f00904fb9ff21e1d8e83afd63aa441` |
| `test_ATAC.h5ad` | COVID scATAC-seq counts: 29,021 cells × 230,356 peaks | 260 MB | `9053f14b396ae2553677d00b69ba18a8` |
| `test_multiomics_integrated.h5ad` | COVID RNA + ATAC integrated by scGLUE: 59,010 cells | 54 MB | `3d474460e651fe61cc0751a3911512d1` |
| `soundlife_RNA.h5ad` | Sound Life scRNA-seq counts: 32,000 cells × 18,456 genes | 137 MB | `0446c8728d187cb95f4e55cc5c6f1339` |
| `soundlife_FLOW.h5ad` | Sound Life spectral flow: 32,000 events × 24 markers | 7 MB | `b57f57b32d39161ac8081a5f1287f33e` |

## Download

Fetch the files into a `data/` folder. Keep the `-f` in `curl`: without it, a failed download is saved as an error page under the `.h5ad` name.

=== "curl"

    ```bash
    mkdir -p data
    curl -fL -o data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
    curl -fL -o data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
    # optional — used by the RNA + ATAC (scGLUE) tutorial:
    curl -fL -o data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
    # used by the RNA + flow cytometry (MaxFuse) tutorial:
    curl -fL -o data/soundlife_RNA.h5ad  "https://zenodo.org/records/24451067/files/soundlife_RNA.h5ad?download=1"
    curl -fL -o data/soundlife_FLOW.h5ad "https://zenodo.org/records/24451067/files/soundlife_FLOW.h5ad?download=1"
    ```

=== "wget"

    ```bash
    mkdir -p data
    wget -O data/test_RNA.h5ad  "https://zenodo.org/records/24451067/files/test_RNA.h5ad?download=1"
    wget -O data/test_ATAC.h5ad "https://zenodo.org/records/24451067/files/test_ATAC.h5ad?download=1"
    # optional — used by the RNA + ATAC (scGLUE) tutorial:
    wget -O data/test_multiomics_integrated.h5ad "https://zenodo.org/records/24451067/files/test_multiomics_integrated.h5ad?download=1"
    # used by the RNA + flow cytometry (MaxFuse) tutorial:
    wget -O data/soundlife_RNA.h5ad  "https://zenodo.org/records/24451067/files/soundlife_RNA.h5ad?download=1"
    wget -O data/soundlife_FLOW.h5ad "https://zenodo.org/records/24451067/files/soundlife_FLOW.h5ad?download=1"
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

To check the files, compare their MD5s with the table. A partial download otherwise shows up later as an HDF5 error such as "bad symbol table node signature".

```bash
md5sum data/*.h5ad     # on macOS: md5 data/*.h5ad
```

## License

The record is licensed [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/), for noncommercial use. The Sound Life data come from the Allen Institute for Immunology and are redistributed under the [Allen Institute Terms of Use](https://alleninstitute.org/terms-of-use/), which permit research and other noncommercial use.

## Citation

> Jiang, H., & Ji, H. *SampleDisco tutorial demo data: COVID-19 scRNA-seq/scATAC-seq and Sound Life scRNA-seq/flow cytometry subsets.* Zenodo. https://doi.org/10.5281/zenodo.22549610

If you use the Sound Life files, also cite the source:

> Gong, Q. et al. Longitudinal multi-omic immune profiling reveals age-related immune cell dynamics in healthy adults. *Nature* (2025). https://doi.org/10.1038/s41586-025-09686-5 — data: Allen Institute for Immunology, Human Immune System Explorer, https://doi.org/10.57785/bn17-3h16
