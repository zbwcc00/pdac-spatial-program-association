# GSE277116 full-package spatial QC

## Scope and design

This QC reads all 28 nested Visium packages directly from the GEO outer archive. Two packages lacking expression matrices/barcodes/features are excluded before quantitative QC. Among the 26 complete packages, GEO metadata identifies 21 tumour samples (8 FFPE and 13 frozen) and 5 lymph-node samples (4 FFPE and 1 frozen). No patient identifier is available; samples remain sample-level units.

## Technical summary

- **package_status**: complete=26; incomplete=2
- **complete_packages_tissue_class**: lymph_node=5; tumor=21
- **complete_packages_sample_type_norm**: ffpe=12; frozen=14
- **complete_packages_technical_candidate**: 1=26
- **complete_packages_program_scoreability_candidate**: 0=4; 1=22
- **complete_packages_tumor_replication_candidate**: 0=5; 1=21
- **complete_packages_tumor_program_scoreability_candidate**: 0=8; 1=18

- Complete packages quantitatively assessed: **26**.
- Candidate tumour replication packages under the prespecified technical screen: **21**.
- Tumour packages meeting the additional marker-detection scoreability screen: **18**.
- Across complete packages, median of per-package median tissue UMI: **1457.0**; median of per-package median detected genes: **957.8**.

## Candidate screen

The technical-candidate flag requires coordinate/barcode overlap >=0.95, >=100 tissue spots, median tissue UMI >=10, median detected genes >=10, largest 8-neighbour tissue component >=50% of tissue spots, and at least 50% of Tfh-like and strict mregDC-like marker genes detected in tissue. The separate program-scoreability flag requires at least 3/6 strict mregDC-like and 5/9 Tfh-like markers to be detected in >=1% of tissue spots. These are permissive operational gates for deciding which packages merit program-level replication QC; they are not biological truth thresholds and must be reviewed by FFPE/frozen stratum.

Only packages labeled Tumor/tumor are eligible for the replication manifest. Lymph-node packages are retained as anatomical controls and are excluded from the PDAC tumour replication endpoint. The CODA annotation resource covers only the J1568 FFPE subset and has no positive TLS labels, so annotation cannot be used as a TLS-positive endpoint.

## Files

- `GSE277116_package_qc.tsv`: one row per nested package with expression, coordinate, complexity, marker and spatial-connectivity metrics.
- `GSE277116_marker_detection_by_sample.tsv`: per-sample tissue detection fractions for every frozen-program marker.
- `GSE277116_tumor_replication_manifest.tsv`: tumour-only sample-level candidate manifest.
- `GSE277116_full_package_qc_report.md`: this audit and interpretation.
