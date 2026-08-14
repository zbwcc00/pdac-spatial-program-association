# Analysis timeline and confirmatory-status record

This record distinguishes the locked primary operational analysis from later
robustness analyses. It is not a prospective preregistration and must not be
read as one.

| Version/date | Change | Rationale and evidentiary status |
|---|---|---|
| Exploratory phase, before v2 | Public-data discovery, program detectability, coordinate/package audits, and exploratory spatial analysis | Feasibility work; not confirmatory. |
| v2 unified primary rerun | Fixed mregDC-like and Tfh-like lists, covariates, six-neighbour focal-spot-excluded fields, DC-core-above-mean mask, GSE278687 patient aggregation, and GSE277116 sample-level role | Locked before the unified rerun, but not prospectively preregistered. This defines the primary operational estimand. |
| v3 | Increased hierarchy-matched block-null draws to 999 and global-Moran graph-sensitivity draws to 1,000 in four deterministic batches | Increased Monte Carlo resolution and corrected the patient-level null hierarchy. Sensitivity/randomization refinement. |
| v4 | Added DC-core-mask-Moran graph sensitivity and continuous primary-program scoreability thresholds | Added after identification of global-to-mask Moran mismatch and platform sparsity concerns. Post-hoc sensitivity work. |
| v6 | Consolidated the manuscript, release record, and figures | Reporting and reproducibility consolidation. |
| v7 | Added whole-tissue/DC-core-conditioning, DC-core-threshold, k=4/6/12, coordinate-distance, VIF, and all-patient-display analyses | Added in response to methodological review. Post-hoc sensitivity/descriptive analyses; they neither replace nor confirm the locked primary endpoint. |
| 2026-08-14 single-cell attribution addendum | Audited the complete GSE202051 author-annotated integrated object and scored unchanged frozen programs using patient-level author-label comparisons | Added after the spatial analyses to resolve a cell-identity evidence gap. It is program-attribution/non-specificity QC, not a prospective confirmation, validation of the spatial association, or a cell-interaction analysis. |

All scripts, deterministic seeds, output paths, and SHA-256 records in the local
release manifest must be reconciled before public release. A public repository,
environment file, release tag, and Zenodo DOI remain required before submission.
