# Supplementary methods and data dictionary (v7, single-cell attribution update)

## Monte Carlo reporting

The GSE278687 array-block randomization test contains 999 draws and is aggregated section -> patient -> cohort median at every draw. The two-sided add-one P value is `(1 + number of draws with abs(null) >= abs(observed))/(B + 1)`. It had 0/999 extreme draws and P=0.001.

Global-Moran and mask-Moran graph sensitivities each concatenate four separately seeded deterministic 250-draw batches (seeds 20360813, 20460813, 20560813 and 20660813). Every graph result reports 0/1,000 extreme draws and P=0.000999 where applicable. These deterministic seed streams should not be described as a claim of philosophical or hardware-level independence.

## Complementary Moran targets

S8-S9 calibrate whole-section outcome Moran's I and audit DC-core-mask Moran's I. S24-S25 instead calibrate outcome Moran's I on the induced DC-core-mask subgraph and audit whole-section Moran's I. A single one-parameter graph filter cannot generally match both targets simultaneously. Both graph sensitivities hold the observed mregDC-like field fixed and generate only Tfh-like outcome surrogates; they do not represent a joint bivariate spatial null. Thus both graph procedures are conditional sensitivity evidence, while the hierarchy-matched patient block null is the principal randomization evidence.

## Scoreability and reference programs

S17/S20 report gene coverage and score zero fractions. S22/S23 use post-hoc scoreability thresholds specified before that sensitivity rerun (joint primary-program nonzero-score proportions 0.20, 0.30, 0.40 and 0.50). Reference-program sparsity does not justify interpreting weak reference effects as negative cell-abundance controls.

## Conditional, scale, and joint-model sensitivity analyses

S26 reports post-hoc whole-tissue and DC-core-threshold analyses with and without DC-core residualization, plus k=4/6/12 local fields. The specified primary estimand deliberately conditions twice on DC-core: an above-mean DC-core mask and DC-core residualization of mregDC-like. It targets residual mregDC-like variation within DC-core-enriched tissue, not a tumour-wide association. S27 reports section-specific raw-coordinate neighbour-distance distributions. GSE278687 coordinates are archive image pixels and GSE277116 coordinates are archive array units; they are not treated as common physical micrometre distances. S28-S29 report descriptive joint-model variance-inflation factors and standardized coefficients. These coefficients do not model spatial residual dependence and are not independent-effect tests.

## Analysis evolution

Supplementary Data S2 records when and why the v2 primary operational analysis and later v3-v7 refinements were introduced. No prospective preregistration is claimed.

## Single-cell QC boundary

Supplementary Figure S1 and its supporting tables are marker-gated gene-detectability/gating-feasibility QC only. Gates use score genes, so expected score enrichment is circular and is not independent cellular validation, localization or intercellular-association evidence.

GSE202051 is a separate, author-annotated PDAC single-cell reference. The exact GEO-length compressed archive passed gzip and SHA-256 checks before use. Author level-1 to level-3 labels were retained without relabelling. The fixed spatial-study gene lists were scored on author-provided log-normalized `X` as a program mean minus expression-bin-matched control mean (25 bins, up to 50 control genes per target gene; seed 20260814); unadjusted program means are reported as a sensitivity. Patient-level within-label medians are the inferential units. The sole designated primary single-cell comparison is the two-sided paired Wilcoxon test of matched-control mregDC-like score in author-labelled activated DC versus cDC2; it was designated for this module but not prospectively preregistered. All other comparisons, score variants, leave-one-gene-out analyses and minimum-cell-per-label sensitivities are exploratory and receive Benjamini-Hochberg FDR adjustment across Supplementary Table S31. Median-difference 95% confidence intervals use 10,000 bootstrap resamples. `counts` and `X` had the same sparse zero pattern, allowing detection to be reported against raw UMI counts.

The dataset contains author-labelled activated DCs but no author-defined Tfh label. Consequently, it can support only higher fixed mregDC-like scores in activated DCs and document the Tfh-like program's non-specificity boundary. It cannot validate the spatial score-field association, Tfh identity, cellular contact, TLS morphology, mechanism, or causality. S30 records the author-label crosswalk; S31 gives patient-level primary and sensitivity comparisons; S32 gives program-gene detection and scoreability by author label; S33 gives activated-DC versus cDC2 sensitivities requiring at least 3 or 5 cells in each label per patient.

## Supplementary figure legends

### Supplementary Figure S1. Marker-gated program detectability and gating-feasibility QC in GSE217845.

This marker-gated analysis documents primary-program gene detectability and feasibility of the stated gates in an independent public PDAC count matrix. Because gate definitions partly use the scored genes, score enrichment is circular and the figure is not evidence of cellular identity, localization, cellular interaction, tissue organization, or an independent single-cell validation.

### Supplementary Figure S2. Representative GSE277116 spatial maps.

Representative FFPE and frozen tumour packages display the tissue image and the corresponding local residual fields under the specified focal-spot-excluded six-neighbour analysis. These are sample-level technical-replication displays because public patient identifiers and package-to-patient linkage are unavailable. They do not provide pathology-defined tissue regions or cell annotations.

### Supplementary Figure S3. Complementary Moran calibration and cross-scale audits.

For each cohort and each 1,000-draw graph sensitivity, the calibrated target has near-zero per-draw section-median absolute Moran's I mismatch, whereas the uncalibrated alternative scale is an audit. Global-targeted batches match whole-section outcome Moran's I; mask-targeted batches match Moran's I on the induced DC-core-mask subgraph. Neither procedure is asserted to match both scales simultaneously.

### Supplementary Figure S4. GSE277116 full-package spatial QC.

All tumour-designated packages passing coordinate, barcode and expression-matrix audit are summarized by tissue-spot count, library size, detected-gene count and per-program marker detection. The main technical-replication set comprises 18 complete scoreable tumour packages; this figure documents platform and package-level constraints rather than biological cell abundance.

### Supplementary Figure S5. Post-hoc primary-program scoreability threshold sensitivity.

The specified local correlation is summarized after retaining units with joint mregDC-like and Tfh-like nonzero-score proportions at or above 0.20, 0.30, 0.40, or 0.50. Points are cohort medians with percentile bootstrap 95% confidence intervals; labels give retained unit counts. Declining counts, especially in GSE277116, mean these are descriptive robustness summaries and not replacement primary analyses.

### Supplementary Figure S6. GSE278687 all-patient local-field atlas.

One section nearest each patient's patient-level primary correlation is shown for all 18 patients, ordered by that patient-level effect. The display is the DC-core-mask joint field `min[z(local mregDC-like residual), z(local Tfh-like residual)]`; it is a visual audit of score-defined fields, not a cell map, tissue-organization map, or an additional inferential analysis.

### Supplementary Figure S7. Post-hoc DC-core conditioning and mask-threshold sensitivity.

Whole-tissue and DC-core masks above the within-section median, mean, or 75th percentile are shown with and without DC-core residualization of the mregDC-like score. Points are patient- or sample-level medians with percentile bootstrap 95% confidence intervals. This post-hoc analysis tests dependence on the conditional estimand definition, including its double DC-core conditioning; it does not replace the primary endpoint.

### Supplementary Figure S8. Neighbourhood-size sensitivity and joint-model collinearity audit.

The primary adjusted mean-DC-core-mask correlation is summarized for k=4, 6, and 12 nearest neighbours. The companion bars give patient- or sample-level median VIF values for the local-field joint model. Coordinate units differ between cohorts and VIF does not address unmodeled spatial residual correlation; the joint model is descriptive.

### Supplementary Figure S9. Fixed-program attribution in an independent author-annotated PDAC single-cell reference (GSE202051).

Author-provided cell labels were used without relabelling. Each point is a within-patient median score for an author label, not an individual cell. The mregDC-like score is preferentially higher in activated DCs than cDC2, cDC1, pDC, and macrophage states, and this direction persists after `CCL19` deletion. The operational Tfh-like score is distributed across broad T-cell states; the reference contains no author-defined Tfh label. This is program-attribution and non-specificity QC, not validation of spatial proximity, a cell-cell interaction, TLS morphology, or mechanism.
