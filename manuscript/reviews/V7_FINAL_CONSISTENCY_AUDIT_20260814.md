# V7 Final Consistency Audit

**Audit date:** 2026-08-14  
**Scope:** historical internal consistency audit completed before the public
release of the v7 manuscript source/Word render. It is superseded for release
status by the public v1.0.5 package and remains neither an external peer-review
decision nor a claim of prospective confirmation.

## Result

**PASS after the corrections documented below.** No unresolved discrepancy was
found between the locked primary analysis, randomization summaries, v7 text,
main figures, supplementary package, or the then-current verification record.

## Numerical Cross-checks

- GSE278687: 21 sections from 18 patients; patient-level median local r =
  0.4183273, reported as 0.418; bootstrap 95% CI 0.312-0.516; 17/18 positive.
- Principal array-block null: 0/999 at-least-as-extreme draws; two-sided
  add-one Monte Carlo P = 0.001. The manuscript and Figure 1 report this as
  the principal randomization analysis.
- Global- and DC-core-mask-Moran graph sensitivities: 0/1,000 extreme draws;
  add-one P = 0.000999. The manuscript distinguishes their separate
  calibration targets and alternative-scale audits.
- GSE277116: 18 scoreable tumour packages; sample-level median local r =
  0.4136335, reported as 0.414; 95% CI 0.287-0.470; 18/18 positive. It is
  consistently described as external **sample-level technical replication**,
  not patient-level validation.
- GSE202051: the sole designated primary comparison is matched-control
  activated-DC versus cDC2 score among 23 paired patients; mregDC-like median
  difference = 0.4785579, reported as 0.479 (bootstrap 95% CI 0.424-0.668),
  two-sided paired Wilcoxon P = 2.38e-6. Exploratory minimum-cell sensitivities
  retain 9 patients at >=3 cells per label (difference 0.718, 95% CI
  0.359-0.933, FDR q=0.00641) and 7 at >=5 (difference 0.718, 95% CI
  0.359-0.747, q=0.0169). The author reference has no author-defined Tfh
  label, which is consistently retained as a non-specificity boundary.

## Claims and Methods

- The local field definition consistently excludes the focal spot and uses six
  nearest neighbours.
- The manuscript explicitly records double conditioning on DC-core, the
  outcome-only nature of graph nulls, incomplete within-block dependence
  preservation, non-comparable physical coordinate scales, and lack of spatial
  residual modelling in the joint model.
- Searches of the active v7 manuscript found no stale `mregDC-Tfh niche`,
  one-sided Wilcoxon, single-cell-spatial-validation, or Tfh-specific-spatial
  claim. TLS, interaction, causal, prognostic and treatment-prediction claims
  are consistently excluded.
- All non-primary single-cell comparisons, score variants, leave-one-gene-out
  tests and minimum-cell sensitivities are explicitly exploratory and have
  Benjamini-Hochberg FDR values in Supplementary Table S31. Supplementary
  Table S33 contains the activated-DC/cDC2 >=3 and >=5 cell-per-label results.
- Supplementary Table S21 contains all six manuscript GEO accessions and a
  cited/verified provenance status for each.

## Reproducibility Package Corrections

- `REPRODUCE_ALL_V7.ps1` now regenerates holdout results, Figure 2/S1-S5,
  Figure 5, and S6-S8 before building the Word manuscript and manifest.
- All 20 scripts invoked by the runbook are content-addressed. The v7 Word
  builder's inherited v5 builder and all five main figure PDF/PNG assets are
  also hashed.
- The verification manifest now recursively hashes every active v7 supplementary file.
  Superseded generated v1/v3 supplemental manifest indexes are removed by the
  v7 synchronizer; the package has one active `supplementary_pack_manifest_v7.json`.
- Figure 5 was relabelled from “holdout validation” to “holdout sensitivity”
  to match its non-inferential role in the Methods and Results.

## Final Integrity Checks

- `local_release_manifest_v7.json`: 111 content-addressed files after the
  single-cell threshold, environment and wording corrections, including
  Supplementary Table S33.
- `release_verification_v7.json`: PASS; 112 checks after the same rebuild,
  including the 2.59-GB compressed GSE202051 archive hash.
- The Word file opens as a valid DOCX and contains 95 paragraphs, 5 tables and
  5 embedded main figures. The five main figures were regenerated during the
  final build and visually inspected for readable labels and scoped claims.

## Deliberately Outstanding Before Submission

GitHub repository creation, a release tag and Zenodo DOI were intentionally
deferred by instruction and are not claimed in the manuscript or manifest.
Before external submission, publish the code/nonrestricted derivative assets,
archive the tagged release, then replace the local-only availability language
with the resulting persistent links and DOI.
