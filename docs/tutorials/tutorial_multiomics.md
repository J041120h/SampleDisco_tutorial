# Tutorial 3: Multi-omics Integration

The multi-omics branch combines RNA and ATAC samples through GLUE-based integration, followed by a shared sample-level embedding and downstream analyses.

## High-level workflow

1. Run GLUE preprocessing and training on RNA and ATAC inputs.
2. Integrate the modalities into a common cell representation.
3. Preprocess the integrated object for sample-level analysis.
4. Optionally search for the optimal cell-type resolution.
5. Build expression and proportion sample embeddings, then run downstream analysis.

## Important wrapper behavior

The multi-omics wrapper documents one major sequencing rule:

```python
# find_optimal_resolution runs BEFORE dimensionality_reduction
```

That ordering matters when you compare candidate resolutions across RNA and ATAC structure.

## Representative outputs

### Integrated sample coordinates

The copied distance-reduction coordinates file mixes RNA and ATAC samples in one table:

| Sample | PC1 | PC2 | PC3 |
| --- | ---: | ---: | ---: |
| `SRR14466474_ATAC` | -1.7112 | -3.4762 | 1.5893 |
| `CoV-34-Aruna_RNA` | 7.5848 | 19.6423 | 11.0849 |
| `CoV-3-Wilk_RNA` | 9.6717 | 4.4960 | -10.9849 |
| `HD-4-Wilk_RNA` | -0.0195 | -9.7368 | 8.5743 |

Full artifact: [expression_dr_coordinates.csv](../resource/data/multiomics/expression_dr_coordinates.csv)

### Modality alignment and resolution summaries

- [resolution_results_rna_expression.csv](../resource/data/multiomics/resolution_results_rna_expression.csv)
- [celltype_modality_distribution.csv](../resource/data/multiomics/celltype_modality_distribution.csv)

### Trajectory-linked differential signals

The copied GAM results include p-values, FDR, effect size, direction, and a `pseudoDEG` flag.

- [trajectory_deg_significant.tsv](../resource/data/multiomics/trajectory_deg_significant.tsv)

## Practical notes

!!! note
    Multi-omics is the best tutorial for understanding how SampleDisc uses both expression-derived and composition-derived sample spaces to compare modalities.

!!! warning
    GLUE integration and gene-activity generation are typically the longest-running steps in the whole project. Plan for GPU access when possible and test on a smaller subset before scaling to a full study collection.

## Runtime estimate

- Small demonstration runs: **roughly an hour** when major intermediates are reused.
- Larger cross-study integrations: **multiple hours**, especially with GLUE training, gene activity, and full downstream analysis enabled.

## Related references

- [Overview: Using Config Files](config_overview.md)
- [Sample Embedding API](../api/embedding.md)
- [Downstream Analysis API](../api/downstream.md)
