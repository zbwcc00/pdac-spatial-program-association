"""Combine four separately seeded 250-draw v4 mask-Moran batches.

The combined test uses all retained draws once.  It recomputes the add-one
two-sided Monte Carlo P value from the concatenated distribution rather than
combining batch-level P values.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "data" / "03_results" / "spatial_nulls_v4_mask_moran"
BATCHES = ("1", "2", "3", "4")


def read_one_row(path):
    with path.open(encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def combine(cohort, label):
    prefix = f"{cohort}_v4_{label}_mask_moran"
    rows = [read_one_row(OUT / f"{prefix}_null_batch{batch}.tsv") for batch in BATCHES]
    draws = np.concatenate([np.loadtxt(OUT / f"{prefix}_draws_batch{batch}.tsv", skiprows=1) for batch in BATCHES])
    audits = [
        np.genfromtxt(OUT / f"{prefix}_audit_batch{batch}.tsv", delimiter="\t", names=True)
        for batch in BATCHES
    ]
    mask_deltas = np.concatenate([audit["mask_moran_section_median_abs_delta"] for audit in audits])
    global_deltas = np.concatenate([audit["global_moran_section_median_abs_delta"] for audit in audits])
    observed = float(rows[0]["observed"])
    extreme = int(np.sum(np.abs(draws) >= abs(observed)))
    result = {
        "cohort": cohort,
        "statistic": rows[0]["statistic"],
        "n_units": int(rows[0]["n_units"]),
        "n_null_draws": int(len(draws)),
        "observed": observed,
        "null_median": float(np.median(draws)),
        "null_ci_low": float(np.quantile(draws, 0.025)),
        "null_ci_high": float(np.quantile(draws, 0.975)),
        "n_as_or_more_extreme": extreme,
        "two_sided_add_one_monte_carlo_p": float((1 + extreme) / (len(draws) + 1)),
        "mask_moran_abs_delta_draw_section_median": float(np.median(mask_deltas)),
        "mask_moran_abs_delta_draw_section_median_max": float(np.max(mask_deltas)),
        "global_moran_abs_delta_draw_section_median": float(np.median(global_deltas)),
        "global_moran_abs_delta_draw_section_median_max": float(np.max(global_deltas)),
        "calibration_target": rows[0]["calibration_target"],
        "batching": "four separately seeded deterministic batches of 250 draws; concatenated before inference",
    }
    with (OUT / f"{prefix}_combined_1000.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter="\t")
        writer.writeheader(); writer.writerow(result)
    np.savetxt(OUT / f"{prefix}_combined_1000_draws.tsv", draws, fmt="%.12g", header="cohort_median_unit_r", comments="")
    return result


def main():
    results = [combine("GSE278687", "patient"), combine("GSE277116", "sample")]
    (OUT / "v4_mask_moran_combined_1000.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
