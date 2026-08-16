# PDAC spatial program association

Code, derived results, figures, and manuscript sources for a public-data study
of a specified local association between mregDC-like and T-cell/lymphoid
transcriptional programs in pancreatic ductal adenocarcinoma (PDAC).

## Scope and interpretation

The primary estimand is the correlation between composition-adjusted,
focal-spot-excluded six-nearest-neighbour score fields among DC-core-above-
within-section-mean spots. GSE278687 provides the primary patient-level
analysis (18 patients; 21 sections). GSE277116 is an external sample-level
technical replication because public patient identifiers are unavailable.
GSE202051 is used only as an independent, author-annotated fixed-program
attribution and specificity QC.

This repository does **not** establish a cellular niche, direct cell-cell
contact, TLS morphology or maturity, a mechanism, causality, prognosis, or
treatment prediction. The spatial claim is conditional on the stated score,
residualization, mask, and smoothing definitions. Details, including the
analysis timeline and all limitations, are in
`manuscript/drafts/Manuscript_full_v7_single_cell_attribution.md`.

## What is versioned

The repository contains the computational scripts, frozen gene programs,
derived tables, figures, manuscript source, and supplementary material. Raw
GEO archives are deliberately excluded. `config/public_inputs.json` records
the accession-specific download URLs, intended analytical role, and the
published byte count/SHA-256 for the GSE202051 archive.

| Input | Role |
| --- | --- |
| GSE278687 | primary patient-level PDAC spatial analysis |
| GSE277116 | sample-level technical replication |
| GSE217847 (GSE217845 PDAC subseries) | limited marker-gated program-detectability QC |
| GSE202051 | author-annotated fixed-program attribution and specificity QC |

## Reproduce the full v7 release

This workflow requires Windows PowerShell, Conda (or Mamba), Python 3.12, an
internet connection to NCBI GEO, and enough free disk space for the raw
archives plus the decompressed GSE202051 H5AD. The compressed GSE202051 input
alone is 2,589,715,190 bytes. Runtime is hardware-dependent; the 999-draw
spatial null analyses are intentionally computationally substantial.

Create and activate the pinned environment:

```powershell
conda env create -f environment.yml
conda activate pdac-spatial-program-association
```

Then run the complete pipeline from the repository root. The first invocation
downloads all registered public inputs, verifies the GSE202051 checksum,
decompresses its H5AD, prepares the required GSE278687 files, recalculates the
primary and sensitivity analyses, rebuilds figures and manuscript files, and
verifies the release manifest.

```powershell
.\RUN_ALL.ps1 -Python ((Get-Command python).Source) -DownloadInputs
```

If inputs have already been downloaded and checked, omit `-DownloadInputs`:

```powershell
.\RUN_ALL.ps1 -Python ((Get-Command python).Source)
```

To download or audit a single accession without analysis:

```powershell
python scripts/fetch_public_inputs.py --download --only GSE202051
python scripts/fetch_public_inputs.py --only GSE202051
```

`RUN_ALL.ps1` is fail-fast. A successful full run ends with
`release_verification.json` reporting `"status": "PASS"`. The final manuscript
is rebuilt at `manuscript/drafts/Manuscript_full_v7_single_cell_attribution.docx`;
main and supplementary figures are under `figures/`.

## Reproducibility records

`environment.yml` and `requirements.txt` define the canonical public Python
environment. `environment-v6.yml` / `requirements-v2.txt` and
`environment-v7-single-cell.yml` / `requirements-v7-single-cell.txt` are
retained as historical observed-runtime records. The historical records are
not a claim that one earlier environment was used for every intermediate
analysis.

`scripts/build_release_manifest.py` hashes the public source-and-results
release. `scripts/verify_release.py --require-inputs` checks those hashes and
the presence of the four local GEO inputs; it also verifies the published
GSE202051 SHA-256. The GitHub Action performs syntax and input-registry checks
only; it intentionally does not download public data or rerun analyses.

## Public synchronized release

Release v1.0.5 synchronizes the public manuscript, supplementary materials, metadata, manifests, verification report, and build scripts at their documented repository paths. Cite the stable Zenodo concept DOI [10.5281/zenodo.21951887](https://doi.org/10.5281/zenodo.21951887), which resolves to the latest published archive. The v1.0.5 package contains the submission manuscript, final supplementary methods/data dictionary, supplementary-pack manifest, public reproducibility documentation, and the SHA-256 release manifest.

## Publication metadata

`LICENSE`, `CITATION.cff`, and `zenodo.json` hold the release metadata. The
tagged GitHub release is archived by Zenodo; the stable concept DOI above is
used consistently in the manuscript and supplementary materials.

## Data and citation

Please cite the original GEO studies and the spatial-statistics methods cited
in the manuscript. This repository redistributes no raw GEO data; accessions,
selection rules, platforms, and source-study citations are audited in
`manuscript/tables/Supplementary_Table_S21_GEO_accession_audit_v1.tsv`.
