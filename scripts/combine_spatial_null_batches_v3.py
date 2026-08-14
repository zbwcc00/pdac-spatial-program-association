"""Combine independent graph-null batches and report 1,000-draw inference."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "data/03_results/spatial_nulls_v3_999"


def read_values(pattern):
    paths = sorted(OUT.glob(pattern))
    if not paths:
        raise RuntimeError(f"No files match {pattern}")
    values = [np.loadtxt(path, skiprows=1) for path in paths]
    return np.concatenate(values), paths


def audit_values(pattern):
    paths = sorted(OUT.glob(pattern))
    records = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            records.extend(csv.DictReader(handle, delimiter="\t"))
    return records, paths


def cohort(name, draw_pattern, audit_pattern, observed):
    draws, draw_paths = read_values(draw_pattern)
    audit, audit_paths = audit_values(audit_pattern)
    n = len(draws); extreme = int(np.sum(np.abs(draws) >= abs(observed)))
    result = {
        "cohort": name, "n_independent_graph_null_draws": n, "observed": observed,
        "n_as_or_more_extreme": extreme,
        "two_sided_monte_carlo_p": float((1 + extreme) / (n + 1)),
        "null_median": float(np.median(draws)), "null_ci_low": float(np.quantile(draws, .025)), "null_ci_high": float(np.quantile(draws, .975)),
        "global_moran_abs_delta_per_draw_section_median": float(np.median([float(row["global_moran_section_median_abs_delta"]) for row in audit])),
        "global_moran_abs_delta_per_draw_section_median_max": float(np.max([float(row["global_moran_section_median_abs_delta"]) for row in audit])),
        "mask_moran_abs_delta_per_draw_section_median": float(np.median([float(row["mask_moran_section_median_abs_delta"]) for row in audit])),
        "mask_moran_abs_delta_per_draw_section_median_max": float(np.max([float(row["mask_moran_section_median_abs_delta"]) for row in audit])),
        "input_draw_files": [path.name for path in draw_paths], "input_audit_files": [path.name for path in audit_paths],
        "graph_null_target": "global outcome Moran I; mask Moran I reported as an audit, not separately calibrated",
    }
    return result


def write_tsv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    results = [
        cohort("GSE278687", "GSE278687_v3_patient_graph_null_draws_batch*.tsv", "GSE278687_v3_patient_graph_null_audit_batch*.tsv", 0.4183273),
        cohort("GSE277116", "GSE277116_v3_main18_graph_null_draws_batch*.tsv", "GSE277116_v3_main18_graph_null_audit_batch*.tsv", 0.4136335),
    ]
    write_tsv(OUT / "v3_combined_1000_graph_nulls.tsv", results)
    (OUT / "v3_combined_1000_graph_nulls.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
