# About and Citation

## What SampleDisc is for

SampleDisc targets **sample-level representation learning** from single-cell RNA, ATAC, and combined multi-omics datasets. The design is especially useful when your main scientific unit is the **sample or cohort**, not the individual cell.

## Current public surface

The repository currently exposes its workflow through:

- `code/SampleDisc.py` for CLI-based execution
- `code/wrapper/wrapper.py` for config-driven orchestration
- modality wrappers in `code/wrapper/`

!!! note
    The repository name (`GenoDistance`) and the user-facing method name (`SampleDisc`) are not yet perfectly aligned in code layout. These docs therefore focus on the actual wrapper functions and CLI entry points available today.

## Documentation version

This site currently labels itself as **v0.1-dev** to indicate that the documentation tracks the present development snapshot rather than a packaged release.

## Citation

```text
[Paper citation placeholder]
```

## Recommended next improvements

- Finalize the Python import namespace so examples can use a stable `import SampleDisc as sd` pattern.
- Publish installation instructions once environment and packaging are stabilized.
- Add rendered figures from benchmark runs to complement the copied CSV and TSV artifacts already included in the docs site.
