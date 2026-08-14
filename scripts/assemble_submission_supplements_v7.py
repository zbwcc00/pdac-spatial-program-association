"""Synchronize v7 single-cell outputs into the versioned supplementary package."""
from __future__ import annotations

import json
import shutil
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SUPPLEMENT = PROJECT / "manuscript/supplementary_v7_single_cell_attribution"


def copy(relative_source: str, relative_destination: str) -> None:
    source = PROJECT / relative_source
    destination = SUPPLEMENT / relative_destination
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> None:
    if not SUPPLEMENT.is_dir():
        raise FileNotFoundError(SUPPLEMENT)
    # Remove generated indexes from superseded assemblies to retain one v7 index.
    for obsolete in (
        "Supplementary_Data_S1_release_manifest_v1.json",
        "supplementary_pack_manifest_v3.json",
    ):
        (SUPPLEMENT / obsolete).unlink(missing_ok=True)
    copy(
        "figures/GSE202051_author_annotated_attribution_v1/Figure_S9_author_annotated_program_attribution.png",
        "Supplementary_Figure_S9_GSE202051_author_annotated_program_attribution.png",
    )
    copy(
        "figures/GSE202051_author_annotated_attribution_v1/Figure_S9_author_annotated_program_attribution.pdf",
        "Supplementary_Figure_S9_GSE202051_author_annotated_program_attribution.pdf",
    )
    copy(
        "data/03_results/GSE202051_author_annotated_attribution_v1/author_label_crosswalk.tsv",
        "Supplementary_Table_S30_GSE202051_author_annotation_crosswalk.tsv",
    )
    copy(
        "data/03_results/GSE202051_author_annotated_attribution_v1/patient_level_primary_and_sensitivity_tests.tsv",
        "Supplementary_Table_S31_GSE202051_patient_level_program_tests.tsv",
    )
    copy(
        "data/03_results/GSE202051_author_annotated_attribution_v1/program_scoreability_by_author_label.tsv",
        "Supplementary_Table_S32_GSE202051_program_scoreability_by_author_label.tsv",
    )
    copy(
        "data/03_results/GSE202051_author_annotated_attribution_v1/mregDC_activated_dc_minimum_cell_sensitivity.tsv",
        "Supplementary_Table_S33_GSE202051_activated_DC_minimum_cell_sensitivity.tsv",
    )
    copy(
        "manuscript/tables/Supplementary_Table_S21_GEO_accession_audit_v1.tsv",
        "Supplementary_Table_S21_GEO_accession_audit.tsv",
    )
    payload = {
        "version": "v7_single_cell_attribution",
        "principal_randomization_evidence": "999-draw patient-level array-block null",
        "global_moran_graph_sensitivity": "1000 draws; target global outcome Moran I",
        "mask_moran_graph_sensitivity": "1000 draws; target DC-core-mask induced-subgraph outcome Moran I",
        "single_cell_program_attribution": "GSE202051 author-annotated reference; higher fixed mregDC-like scores in activated DCs and Tfh-like non-specificity boundary only",
        "release_manifest": "local only; no GitHub/Zenodo DOI claimed",
    }
    (SUPPLEMENT / "supplementary_pack_manifest_v7.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(SUPPLEMENT)


if __name__ == "__main__":
    main()
