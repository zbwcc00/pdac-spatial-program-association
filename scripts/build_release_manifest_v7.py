"""Create the local, content-addressed v7 pre-publication release manifest."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "local_release_manifest_v7.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(relative: str) -> dict[str, object]:
    path = PROJECT / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    return {"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)}


def main() -> None:
    files = [
        "README.md",
        "RUN_V7_REPRODUCTION.md",
        "REPRODUCE_ALL_V7.ps1",
        "RELEASE_CHECKLIST_v7.md",
        "environment-v6.yml",
        "requirements-v2.txt",
        "environment-v7-single-cell.yml",
        "requirements-v7-single-cell.txt",
        "scripts/run_unified_primary_pipeline_v2.py",
        "scripts/run_gse278687_block_null_v3_999.py",
        "scripts/run_spatial_nulls_v3_999.py",
        "scripts/combine_spatial_null_batches_v3.py",
        "scripts/run_mask_moran_sensitivity_v4.py",
        "scripts/combine_mask_moran_batches_v4.py",
        "scripts/run_scoreability_threshold_sensitivity_v4.py",
        "scripts/run_conditional_scale_joint_sensitivities_v7.py",
        "scripts/run_spatial_holdout_validation.py",
        "scripts/run_gse217845_marker_gated_single_cell_validation.py",
        "scripts/audit_gse202051_author_reference.py",
        "scripts/run_gse202051_author_annotated_attribution.py",
        "scripts/generate_main_figures_v4.py",
        "scripts/generate_spatial_and_supplementary_figures_v6.py",
        "scripts/generate_conditional_scale_figures_v7.py",
        "scripts/generate_holdout_figure_v1.py",
        "scripts/build_manuscript_docx_v7.py",
        "scripts/build_manuscript_docx_v5.py",
        "scripts/assemble_submission_supplements_v7.py",
        "scripts/build_release_manifest_v7.py",
        "scripts/verify_release_v7.py",
        "data/03_results/GSE154778_program_freeze/frozen_programs.json",
        "data/03_results/unified_primary_pipeline_v2/pipeline_manifest.json",
        "data/03_results/unified_primary_pipeline_v2/GSE278687_v2_per_patient.tsv",
        "data/03_results/unified_primary_pipeline_v2/GSE277116_v2_per_package.tsv",
        "data/03_results/spatial_nulls_v3_999/GSE278687_v3_patient_block_null_999.tsv",
        "data/03_results/spatial_nulls_v3_999/v3_combined_1000_graph_nulls.tsv",
        "data/03_results/spatial_nulls_v4_mask_moran/GSE278687_v4_patient_mask_moran_combined_1000.tsv",
        "data/03_results/spatial_nulls_v4_mask_moran/GSE277116_v4_sample_mask_moran_combined_1000.tsv",
        "data/03_results/conditional_scale_joint_sensitivities_v7/cohort_conditional_scale_summary.tsv",
        "data/03_results/GSE202051_author_reference_audit/gse202051_author_reference_audit.json",
        "data/03_results/GSE202051_author_annotated_attribution_v1/analysis_protocol.json",
        "data/03_results/GSE202051_author_annotated_attribution_v1/patient_level_primary_and_sensitivity_tests.tsv",
        "data/03_results/GSE202051_author_annotated_attribution_v1/mregDC_activated_dc_minimum_cell_sensitivity.tsv",
        "data/03_results/GSE202051_author_annotated_attribution_v1/program_scoreability_by_author_label.tsv",
        "manuscript/analysis_timeline_v1.md",
        "manuscript/external_single_cell_attribution_audit_v3.md",
        "manuscript/reviews/V7_FINAL_CONSISTENCY_AUDIT_20260814.md",
        "manuscript/reviews/REFERENCE_EVIDENCE_AUDIT_v1.md",
        "manuscript/tables/Supplementary_Table_S21_GEO_accession_audit_v1.tsv",
        "manuscript/drafts/Manuscript_full_v7_single_cell_attribution.md",
        "manuscript/drafts/Manuscript_full_v7_single_cell_attribution.docx",
        "figures/main_v4/Figure1_v4_primary_and_spatial_nulls.png",
        "figures/main_v4/Figure1_v4_primary_and_spatial_nulls.pdf",
        "figures/submission_v6/Figure2_v6_GSE278687_representative_spatial_maps.png",
        "figures/submission_v6/Figure2_v6_GSE278687_representative_spatial_maps.pdf",
        "figures/main_v4/Figure2_v4_robustness_and_competition.png",
        "figures/main_v4/Figure2_v4_robustness_and_competition.pdf",
        "figures/main_v4/Figure3_v4_external_technical_replication.png",
        "figures/main_v4/Figure3_v4_external_technical_replication.pdf",
        "figures/main_v2/Figure5_v1_spatial_holdout_validation.png",
        "figures/main_v2/Figure5_v1_spatial_holdout_validation.pdf",
    ]
    supplementary = PROJECT / "manuscript/supplementary_v7_single_cell_attribution"
    supplementary_files = [
        path.relative_to(PROJECT).as_posix()
        for path in sorted(supplementary.rglob("*"))
        if path.is_file()
    ]
    raw_inputs = [
        {
            "accession": "GSE202051",
            "path": "data/00_raw/single_cell/GSE202051/GSE202051_totaldata-final-toshare.h5ad.gz",
            "bytes": 2589715190,
            "sha256": "e9b43be8b5bf8d7a606b9cb3c972b1bee93826d2c48f7ac565aea0fe57bb1a43",
            "role": "author-annotated program-attribution and specificity QC only",
        }
    ]
    payload = {
        "manifest_schema": "local-content-addressed-pdac-spatial-release-v7",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "release_status": {
            "state": "local pre-publication record; not externally published",
            "github_repository": None,
            "github_release_tag": None,
            "zenodo_doi": None,
            "publication_authorization": "not granted",
        },
        "primary_endpoint": {
            "definition": "composition-adjusted focal-spot-excluded six-nearest-neighbour local mregDC-like--Tfh-like Pearson correlation among DC-core-above-within-section-mean spots",
            "GSE278687_aggregation": "section effects -> patient medians -> cohort median (18 patients, 21 sections)",
            "GSE277116_aggregation": "sample-level cohort median (18 complete tumour packages; no public patient identifiers)",
            "interpretive_boundary": "score-defined local association only; not a cell interaction, TLS morphology, mechanism, causality, prognosis, or treatment prediction",
        },
        "single_cell_module": {
            "dataset": "GSE202051",
            "role": "independent author-annotated fixed-program attribution and specificity QC",
            "inferential_unit": "within-patient median score per author label",
            "sole_designated_primary_test": "matched-control mregDC-like score: author-labelled activated DC versus cDC2; two-sided paired Wilcoxon signed-rank test",
            "exploratory_family": "all other comparisons, scoring variants, leave-one-gene-out tests and >=3/>=5 cells-per-label sensitivities; Benjamini-Hochberg FDR adjustment",
            "boundary": "no author-defined Tfh label; not spatial validation or Tfh identity assignment",
        },
        "analysis_status": {
            "primary_operational_analysis": "locked before unified v2 rerun, not prospectively preregistered",
            "later_sensitivities": "post-hoc and documented in manuscript/analysis_timeline_v1.md",
        },
        "environments": {
            "spatial_historical_record": "environment-v6.yml and requirements-v2.txt",
            "single_cell_observed_runtime": "environment-v7-single-cell.yml and requirements-v7-single-cell.txt",
            "manifest_builder_runtime": {"python": sys.version, "platform": platform.platform()},
        },
        "raw_input_policy": "Raw public archives remain local and are excluded from any future source repository. Hashes are retained for integrity verification.",
        "raw_inputs": raw_inputs,
        "content_addressed_files": [record(path) for path in dict.fromkeys(files + supplementary_files)],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
