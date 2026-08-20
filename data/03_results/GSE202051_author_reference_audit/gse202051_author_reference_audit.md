# GSE202051 author-reference audit

- Remote expected bytes: `2589715190`.
- Local integrated object bytes: `2589715190`.
- Exact-length complete: `True`.
- Decision: Eligible for fixed-program attribution/specificity QC: the complete integrated object carries final author labels and provenance fields.

## Individual local objects

### GSE202051_adata_010nuc_10x.h5ad
- Cells x genes: `2607 x 33538`.
- Final author label fields: `none`.
- Provenance fields: `pid`; PID categories: `010nuc_10x`.
- Eligibility: ineligible: no final author cell-type/state field in this object

### GSE202051_adata_010orgCRT_10x.h5ad
- Cells x genes: `341 x 33538`.
- Final author label fields: `none`.
- Provenance fields: `pid`; PID categories: `010orgCRT_10x`.
- Eligibility: ineligible: no final author cell-type/state field in this object

### GSE202051_totaldata-final-toshare.h5ad
- Cells x genes: `224988 x 22164`.
- Final author label fields: `Level 1 Annotation, Level 2 Annotation, Level 3 Annotation, celltypes, new_celltypes`.
- Provenance fields: `pid, sampleid, treatment_status, new_treatment`; PID categories: `none`.
- Eligibility: eligible in principle: requires barcode and provenance audit before scoring

The audit is a data-eligibility check only; no biological inference has been calculated.
