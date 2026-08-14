"""Predefined scoreability-threshold sensitivity for the locked primary endpoint.

Thresholds apply to the fraction of tissue spots with a nonzero score for both
fixed primary programs.  GSE278687 uses the patient median across all sections
for each scoreability measure and then retains the locked patient-level effect;
GSE277116 uses the individual package.  This avoids choosing individual
sections after examining their association estimate.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCOREABILITY = PROJECT / "data" / "03_results" / "program_scoreability_v3" / "program_scoreability_by_sample.tsv"
LOCKED = PROJECT / "data" / "03_results" / "unified_primary_pipeline_v2"
OUT = PROJECT / "data" / "03_results" / "scoreability_threshold_sensitivity_v4"
OUT.mkdir(parents=True, exist_ok=True)
THRESHOLDS = (0.20, 0.30, 0.40, 0.50)
PRIMARY = {"mregDC_like", "Tfh_like"}
BOOTSTRAP_SEED = 20260813


def read_tsv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def quantile_ci(values, rng, n=10000):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return np.nan, np.nan
    replicates = np.median(rng.choice(values, size=(n, len(values)), replace=True), axis=1)
    return float(np.quantile(replicates, .025)), float(np.quantile(replicates, .975))


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def g278(score_rows):
    fractions = defaultdict(lambda: defaultdict(list))
    for row in score_rows:
        if row["cohort"] == "GSE278687" and row["program"] in PRIMARY:
            fractions[row["sample"]][row["program"]].append(float(row["spot_nonzero_fraction"]))
    section_rows = read_tsv(LOCKED / "GSE278687_v2_per_section.tsv")
    effects = defaultdict(list)
    for row in section_rows:
        effects[row["unit"]].append(float(row["primary_mregDC_Tfh_local_r"]))
    patient_scoreability = defaultdict(lambda: defaultdict(list))
    section_patient = {row["sample"]: row["unit"] for row in section_rows}
    for sample, program_values in fractions.items():
        for program, value in program_values.items():
            patient_scoreability[section_patient[sample]][program].append(float(np.median(value)))
    units = []
    for patient in sorted(effects):
        units.append({
            "cohort": "GSE278687", "unit": patient,
            "effect": float(np.median(effects[patient])),
            "mregDC_nonzero_fraction": float(np.median(patient_scoreability[patient]["mregDC_like"])),
            "Tfh_nonzero_fraction": float(np.median(patient_scoreability[patient]["Tfh_like"])),
            "aggregation": "patient median across sections for scoreability and effect",
        })
    return units


def g277(score_rows):
    fractions = defaultdict(dict)
    for row in score_rows:
        if row["cohort"] == "GSE277116" and row["program"] in PRIMARY:
            fractions[row["sample"]][row["program"]] = float(row["spot_nonzero_fraction"])
    package_rows = read_tsv(LOCKED / "GSE277116_v2_per_package.tsv")
    units = []
    for row in package_rows:
        if row["analysis_set"] != "main_18":
            continue
        # Scoreability was recorded using the study's human-readable package
        # identifier (e.g. pdac2), whereas the locked results retain GSM IDs.
        sample = row["sample"]
        units.append({
            "cohort": "GSE277116", "unit": row["gsm"],
            "effect": float(row["primary_mregDC_Tfh_local_r"]),
            "mregDC_nonzero_fraction": fractions[sample]["mregDC_like"],
            "Tfh_nonzero_fraction": fractions[sample]["Tfh_like"],
            "aggregation": "individual technical-replication package",
        })
    return units


def summarize(units, cohort, rng):
    rows = []
    for threshold in THRESHOLDS:
        kept = [unit for unit in units if min(unit["mregDC_nonzero_fraction"], unit["Tfh_nonzero_fraction"]) >= threshold]
        values = [unit["effect"] for unit in kept]
        low, high = quantile_ci(values, rng)
        rows.append({
            "cohort": cohort,
            "threshold_rule": "both fixed primary program nonzero fractions >= threshold",
            "threshold": threshold,
            "n_eligible_units": len(kept),
            "n_total_units": len(units),
            "median_primary_local_r": float(np.median(values)) if values else np.nan,
            "bootstrap_ci_low": low,
            "bootstrap_ci_high": high,
            "n_positive": int(np.sum(np.asarray(values) > 0)) if values else 0,
            "selection_note": "Thresholds were predefined; no competitive-program scoreability was used for selection.",
        })
    return rows


def main():
    score_rows = read_tsv(SCOREABILITY)
    g278_units = g278(score_rows)
    g277_units = g277(score_rows)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    summary = summarize(g278_units, "GSE278687", rng) + summarize(g277_units, "GSE277116", rng)
    write_tsv(OUT / "scoreability_units_v4.tsv", g278_units + g277_units)
    write_tsv(OUT / "scoreability_threshold_summary_v4.tsv", summary)
    print("\n".join("\t".join(str(value) for value in row.values()) for row in summary))


if __name__ == "__main__":
    main()
