# Independent PDAC single-cell program-attribution audit (v3)

## Scope

This audit tests only whether the fixed spatial-study programs are detectably expressed at higher levels in author-annotated cell states. It does **not** validate spatial proximity, demonstrate a cell-cell interaction or TLS, establish mechanism or causality, or evaluate prognosis or treatment response. The locked spatial estimand is unchanged.

## Reference qualification

The exact-length GEO integrated object from GSE202051 (Hwang et al., *Nature Genetics*, 2022; PMID 35902743) passed local integrity checks: compressed size 2,589,715,190 bytes; gzip SHA-256 `e9b43be8b5bf8d7a606b9cb3c972b1bee93826d2c48f7ac565aea0fe57bb1a43`; decompressed H5AD 8,321,198,912 bytes. It contains 224,988 cells, 22,164 genes, unique cell identifiers, 43 patients, 44 samples, treatment fields, and author fields `Level 1 Annotation`, `Level 2 Annotation`, and `Level 3 Annotation`.

The author level-3 annotation contains activated DC (91 cells from 26 patients), cDC2 (1,008 cells from 37 patients), cDC1, pDC, macrophage, CD4+ T, Treg and CD8-state labels. It contains **no author-defined Tfh label**. This asymmetry determines the permissible interpretation below.

## Locked scoring and inferential unit

Frozen mregDC-like and Tfh-like lists were transferred unchanged from the spatial analysis. The primary score was the mean author-provided log-normalized `X` expression of the program genes minus an expression-bin-matched control-gene mean (25 bins, up to 50 controls per target gene, seed 20260814). The sole designated primary test in this module was the two-sided paired Wilcoxon comparison of this mregDC-like score between author-labelled activated DC and cDC2, using each patient's within-label median. It was designated for this module but was not prospectively preregistered. All other label comparisons, unadjusted-score and gene-deletion variants, and minimum-cell sensitivities are exploratory and have Benjamini-Hochberg FDR values reported in Supplementary Table S31. `X` and the raw-UMI `counts` layer had identical sparse zero patterns, allowing gene detection to be reported against raw-count detection. No cell-level P value was used.

## Findings that can be used

### mregDC-like attribution

The mregDC-like score was higher in author-labelled activated DC than cDC2 in 23 paired patients: median paired difference 0.479 (bootstrap 95% CI 0.424-0.668); two-sided paired Wilcoxon P=2.38e-6. With both labels required to have at least 3 or 5 cells per patient, the result remained positive in 9 patients (difference 0.718, 95% CI 0.359-0.933; P=0.00391; exploratory FDR q=0.00641) and 7 patients (difference 0.718, 95% CI 0.359-0.747; P=0.0156; exploratory FDR q=0.0169), respectively. The result also held against cDC1, pDC and macrophages. It persisted under the unadjusted score (difference 0.473; exploratory P=2.38e-6). All six mregDC-like genes were available.

`CCL19` was detected in only 1/91 activated DC cells (1.1%). Removing it increased the activated-DC versus cDC2 paired difference to 0.573 (exploratory two-sided P=2.38e-6; FDR q=1.40e-5). Leave-one-gene-out remained directionally positive; deletion of LAMP3 produced the greatest attenuation (difference 0.280; exploratory two-sided P=1.53e-4). These score sensitivities do not establish molecular contribution or mechanism. They support **higher mregDC-like scores in the author-defined activated DC compartment**, not a spatial interaction or recruitment mechanism.

### Tfh-like specificity boundary

Because no author Tfh category is available, this dataset cannot validate Tfh identity. Exploratorily, the Tfh-like score was higher in CD4+ T cells than Tregs (25 paired patients; median difference 0.061, two-sided P=0.0051; FDR q=0.00770) and CD8+ T cells (28 paired patients; difference 0.064, two-sided P=0.00058; FDR q=0.00114), but the CD4-versus-Treg effect was small. It was retained with the unadjusted score (difference 0.066; exploratory two-sided P=0.0125).

Crucially, this contrast is not stable as a Tfh-specific feature: deleting `IL7R` reversed the CD4-versus-Treg median difference (-0.083); this exploratory result is recorded in the accompanying output table and does not support a directional Tfh-specific interpretation. Canonical discriminating components were sparse in author CD4+ T cells: CXCL13 0.22%, PDCD1 1.31%, and BCL6 1.74% detection; IL7R was detected in 72.1%. Therefore, the appropriate interpretation is a **broader T-cell/lymphoid program with incomplete Tfh specificity**, not a Tfh identity assignment.

## Required manuscript wording

Use: "In an independent author-annotated PDAC single-cell reference, the fixed mregDC-like score was higher in author-labelled activated dendritic cells, whereas the operational Tfh-like score was distributed across broader CD4/T-cell states and lacked an author-defined Tfh identity reference."

Do not use: "single-cell validation of the spatial association," "mregDC-Tfh interaction," "Tfh niche," "TLS mechanism," or any causal language.

## Files

- Script: `scripts/run_gse202051_author_annotated_attribution.py`
- Results: `data/03_results/GSE202051_author_annotated_attribution_v1/`
- Supplementary Figure S9 source: `figures/GSE202051_author_annotated_attribution_v1/Figure_S9_author_annotated_program_attribution.png`
- Structural audit: `data/03_results/GSE202051_author_reference_audit/gse202051_author_reference_audit.md`
