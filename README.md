# PLP1/PMD Digital Twin (Module A)

Mechanism-to-phenotype model for PLP1 variants. Builds a mouse oligodendrocyte-lineage reference atlas (Marques et al. 2016), defines interpretable cell-state modules, and reports measured disease-state shifts relative to a healthy baseline.

Part of the extended PLP1/PMD project:

- Part 1: [plp1-pmd-mechanism-mapping](https://github.com/muneeb-dotcom/plp1-pmd-mechanism-mapping)
- Module B (therapeutic design): [module-b-therapeutics](https://github.com/muneeb-dotcom/plp1-pmd-mechanism-mapping/tree/main/module-b-therapeutics)
- Module C (validation): [module-c-validation](https://github.com/muneeb-dotcom/plp1-pmd-mechanism-mapping/tree/main/module-c-validation)

## Honest scope

This is **not** a mechanistic simulation. It reports measured module scores from real reference and disease data. It does not extrapolate to mechanism classes with no available data (see `results/severity_analysis.md`).

## Structure

- `scripts/`: pipeline, phases A1-A8 (see filenames a01-a18)
- `twin/pmd_twin.py`: the packaged, callable twin
- `results/`: module definitions, state vectors, analysis writeups
- `figures/`: atlas UMAP, marker validation, maturation axis

## Key findings

- Reference atlas validated against canonical oligodendrocyte markers (Pdgfra, Mog, Plp1)
- Mouse-derived modules validated on independent human data (GSE118257)
- Cross-mechanism expression data unavailable for loss-of-function and duplication variants (see `results/severity_analysis.md`), so the twin reports `no_measured_data` rather than extrapolating

## Reproduce

Large data files are gitignored (~2 GB) and can be regenerated with `scripts/a01` to `scripts/a03`. Run the scripts in order, a01 through a18.
