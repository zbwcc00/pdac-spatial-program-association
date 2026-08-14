"""Extend only the corrected GSE278687 patient-level block null to 999 draws."""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts/run_unified_primary_pipeline_v2.py"
OUT = PROJECT / "data/03_results/spatial_nulls_v3_999"
OUT.mkdir(parents=True, exist_ok=True)
N_DRAWS = 999
SEED = 20260813


def load():
    module = importlib.util.module_from_spec(importlib.util.spec_from_loader("locked_block", loader=None))
    source = PIPELINE.read_text(encoding="utf-8")
    source, count = re.subn(r"^PROJECT = Path\([^\n]+\)$", f"PROJECT = Path({str(PROJECT)!r})", source, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError("Could not substitute path")
    exec(compile(source, str(PIPELINE), "exec"), module.__dict__)
    return module


def main():
    locked = load(); rng = np.random.default_rng(SEED)
    sections, states = [], []
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, rows, cols, patient = locked.read_g278(path)
        scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
        mreg, tfh, mask, neighbours, effects = locked.prepared_fields(scores, coords)
        sections.append({"unit": patient, "primary_mregDC_Tfh_local_r": effects["primary_mregDC_Tfh_local_r"]})
        states.append({"mreg": mreg, "tfh": tfh, "mask": mask, "neighbours": neighbours, "block_ids": locked.block_ids(rows, cols)})
    result = locked.patient_block_null(states, sections, rng, n=N_DRAWS)
    result["n_as_or_more_extreme"] = int(round(result["two_sided_p"] * (N_DRAWS + 1) - 1))
    result["two_sided_monte_carlo_p"] = result.pop("two_sided_p")
    with (OUT / "GSE278687_v3_patient_block_null_999.tsv").open("w", encoding="utf-8") as handle:
        keys = list(result); handle.write("\t".join(keys) + "\n"); handle.write("\t".join(str(result[key]) for key in keys) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
