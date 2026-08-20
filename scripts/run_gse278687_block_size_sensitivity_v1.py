"""Evaluate the GSE278687 patient-level block null across fixed block sizes."""
from __future__ import annotations

import csv
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts" / "run_unified_primary_pipeline_v2.py"
DATA_ROOT = PROJECT / "data"
OUT = PROJECT / "data" / "03_results" / "block_size_sensitivity_v1"
BLOCK_SIZES = (4, 8, 12)
N_DRAWS = 999
SEED = 20260813


def load_locked_pipeline():
    specification = importlib.util.spec_from_file_location("locked_block_size", PIPELINE)
    if specification is None or specification.loader is None:
        raise RuntimeError("Could not load the locked pipeline")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    module.G278 = DATA_ROOT / "01_unpacked" / "spatial" / "GSE278687"
    module.G278_COORDS = DATA_ROOT / "01_unpacked" / "spatial" / "GSE278687_spatial"
    return module


def load_section_states(locked, block_size):
    sections = []
    states = []
    locked.BLOCK = block_size
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, rows, cols, patient = locked.read_g278(path)
        scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
        mreg, tfh, mask, neighbours, effects = locked.prepared_fields(scores, coords)
        sections.append(
            {
                "sample": sample,
                "patient": patient,
                "primary_mregDC_Tfh_local_r": effects["primary_mregDC_Tfh_local_r"],
            }
        )
        states.append(
            {
                "mreg": mreg,
                "tfh": tfh,
                "mask": mask,
                "neighbours": neighbours,
                "block_ids": locked.block_ids(rows, cols),
            }
        )
    return sections, states


def patient_block_null_with_draws(locked, section_states, sections, rng):
    observed_by_patient = defaultdict(list)
    for section in sections:
        observed_by_patient[section["patient"]].append(section["primary_mregDC_Tfh_local_r"])
    observed = float(np.median([np.median(values) for values in observed_by_patient.values()]))
    draws = []
    for _ in range(N_DRAWS):
        null_by_patient = defaultdict(list)
        for section, state in zip(sections, section_states):
            surrogate = locked.block_surrogate(state["tfh"], state["block_ids"], rng)
            null_by_patient[section["patient"]].append(
                locked.local_r(state["mreg"], surrogate, state["neighbours"], state["mask"])
            )
        draws.append(float(np.median([np.median(values) for values in null_by_patient.values()])))
    draws = np.asarray(draws)
    extreme = int(np.sum(np.abs(draws) >= abs(observed)))
    return {
        "observed": observed,
        "n_patients": len(observed_by_patient),
        "n_sections": len(sections),
        "n_null_draws": N_DRAWS,
        "null_median": float(np.median(draws)),
        "null_ci_low": float(np.quantile(draws, 0.025)),
        "null_ci_high": float(np.quantile(draws, 0.975)),
        "n_as_or_more_extreme": extreme,
        "two_sided_add_one_monte_carlo_p": float((1 + extreme) / (N_DRAWS + 1)),
    }, draws


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    locked = load_locked_pipeline()
    summary_rows = []
    draw_rows = []
    for block_size in BLOCK_SIZES:
        sections, states = load_section_states(locked, block_size)
        summary, draws = patient_block_null_with_draws(locked, states, sections, np.random.default_rng(SEED))
        summary_rows.append({"array_coordinate_block_size": block_size, "seed": SEED, **summary})
        draw_rows.extend(
            {
                "array_coordinate_block_size": block_size,
                "seed": SEED,
                "draw": draw_index + 1,
                "cohort_median_patient_local_r": float(draw),
            }
            for draw_index, draw in enumerate(draws)
        )
    write_tsv(OUT / "GSE278687_patient_block_size_sensitivity_summary.tsv", summary_rows)
    write_tsv(OUT / "GSE278687_patient_block_size_sensitivity_null_draws.tsv", draw_rows)
    metadata = {
        "block_sizes": list(BLOCK_SIZES),
        "n_draws_per_block_size": N_DRAWS,
        "seed": SEED,
        "locked_pipeline": str(PIPELINE.relative_to(PROJECT)),
        "description": "Patient-level block-null sensitivity using the locked GSE278687 primary score, mask, local-field, and patient-aggregation definitions.",
    }
    (OUT / "GSE278687_patient_block_size_sensitivity_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary_rows, indent=2))


if __name__ == "__main__":
    main()
