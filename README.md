\# PLP1/PMD Digital Twin (Module A)



Mechanism-to-phenotype model for PLP1 variants. Builds a mouse

oligodendrocyte-lineage reference atlas (Marques et al. 2016),

defines interpretable cell-state modules, and reports measured

disease-state shifts relative to a healthy baseline.



Part of the extended PLP1/PMD project:

\- Part 1: https://github.com/muneeb-dotcom/plp1-pmd-mechanism-mapping

\- Module B: (therapeutic design)

\- Module C: https://github.com/muneeb-dotcom/plp1-pmd-validation



\## Honest scope

NOT a mechanistic simulation. Reports measured module scores from

real reference/disease data. Does not extrapolate to mechanism

classes with no available data (see results/severity\_analysis.md).



\## Structure

\- `scripts/` — pipeline, phases A1-A8 (see filenames a01-a18)

\- `twin/pmd\_twin.py` — the packaged, callable twin

\- `results/` — module definitions, state vectors, analysis writeups

\- `figures/` — atlas UMAP, marker validation, maturation axis



\## Key findings

\- Reference atlas validated against canonical OL markers (Pdgfra, Mog, Plp1)

\- Mouse-derived modules validated on independent human data (GSE118257)

\- Cross-mechanism expression data unavailable for LOF/duplication (see

&#x20; results/severity\_analysis.md) — twin honestly reports "no\_measured\_data"

&#x20; rather than extrapolating



\## Reproduce

Large data files are gitignored (regenerable via scripts/a01-a03,

\~2GB). See scripts in order a01 through a18.

