\# Phase A6 — Severity Analysis (Digital Twin, Module A)



\## Objective

Originally scoped as: fit a model mapping variant features (mechanism,

ddG, dosage) to a measured cell-state vector, using expression data

across multiple mechanism classes.



\## Data availability finding

A targeted search found no deposited transcriptomic (RNA-seq/microarray)

dataset in GEO for either a Plp1-null (loss-of-function) or Plp1-duplication

mouse model. Extensive protein/histology-level literature exists for both

(Karim et al. 2007; Readhead et al. 1994; Cai et al. 2013's Plp1dup mouse),

but none of these studies deposited expression data. GSE111605 (jimpy,

misfolding) appears to be the only mechanism class with public

transcriptomic data available.



\*\*Consequence:\*\* a cross-mechanism expression-based model (as originally

scoped) cannot be fit with real data. Rather than fabricate synthetic

training points, Phase A6 was redirected to a data source that is real

and underused: ClinVar's own phenotype annotations, tested across the

full 371-variant set from Part 1 (14x larger than the Nevin et al. 2017

iPSC panel of 12 lines).



\## Method

Pulled `germline\_classification.trait\_set` for all 375 pathogenic ClinVar

PLP1 records via Entrez esummary. Assigned an ordinal severity score:

\- 0 = Hereditary spastic paraplegia 2 (SPG2, mild end of spectrum)

\- 1 = Pelizaeus-Merzbacher disease (classic/undifferentiated)

\- 2 = Pelizaeus-Merzbacher disease, connatal (severe end)



Merged against Part 1's mechanism classification and Phase 2's FoldX ddG

scores. Tested Spearman correlation within the misfolding class (ddG vs

severity) and attempted the same for duplication (copy number vs severity).



\## Results



\*\*Misfolding class (n=73 variants with both ddG and a severity label):\*\*

Spearman rho = -0.051, p = 0.6696 — not significant.



\*\*Duplication class:\*\* 0 of 63 variants had a usable phenotype label (all

fell into "See cases" or "not provided" trait categories). This test could

not be run.



\## Interpretation



The non-significant misfolding-class correlation should NOT be read as

evidence that structural destabilization is unrelated to clinical

severity. 166 of 172 matched variants (96.5%) fall into a single severity

bucket ("Pelizaeus-Merzbacher disease," undifferentiated) -- ClinVar's

phenotype field functions largely as a diagnosis label, not a graded

clinical outcome measure, and has essentially no variance to correlate

against.



\*\*Honest conclusion:\*\* ClinVar-derived phenotype labels are too coarse to

test a ddG-severity dose-response relationship at this scale. This is

itself a citable methodological finding for anyone attempting a similar

variant-to-severity analysis from public data alone.



\*\*A secondary, speculative hypothesis worth stating:\*\* if a true

biological relationship between destabilization magnitude and severity

is weak or absent even with better-graded data, this would be consistent

with UPR activation behaving as a threshold/switch-like response (known

PERK/eIF2-alpha pathway dynamics) rather than a linear dose-response to

the degree of misfolding. This is not tested here and should be flagged

as a hypothesis for future work, not a finding.



\*\*Duplication class:\*\* the copy-number/severity relationship is

well-documented in the literature (Wolf et al. 2005) but is not testable

from ClinVar's structured fields. A different data source -- the original

patient case tables from primary literature -- would be needed.



\## Status

Phase A6 is closed on this basis. No predictive variant-to-state-vector

model was fit, for the data-availability reasons stated above. The

digital twin (Phase A8) therefore reports only measured, empirical states

per mechanism class, and does not extrapolate to unmeasured mechanism

classes or attempt severity prediction.

