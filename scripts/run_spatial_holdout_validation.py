"""Leakage-aware spatial holdout and leave-one-unit sensitivity analyses.

This script is deliberately downstream of the locked v2 primary pipeline.  It
does not replace its primary estimand or generate a second P value.  Instead,
it asks whether the sign and magnitude remain when a contiguous array block is
held out.  Residualization coefficients and DC-core threshold are fitted only
in the remaining spots, and every held-out local field is constructed only
from held-out spots.  Thus neither expression nor neighbours from the fitting
set enter a held-out local correlation.
"""
from __future__ import annotations

import csv
import importlib.util
import io
import json
import re
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts/run_unified_primary_pipeline_v2.py"
OUT = PROJECT / "data/03_results/spatial_holdout_validation_v1"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260813
BLOCK = 8
MIN_TEST_SPOTS = 25
MIN_MASKED_SPOTS = 20

# The locked source predates a Windows code-page change and contains a literal
# project path that can be decoded differently across shells.  Execute its
# unchanged analytical functions after substituting only that path assignment;
# the original locked file itself is never edited by this validation script.
locked = importlib.util.module_from_spec(importlib.util.spec_from_loader("locked_v2", loader=None))
source = PIPELINE.read_text(encoding="utf-8")
source, replacements = re.subn(
    r'^PROJECT = Path\([^\n]+\)$',
    f"PROJECT = Path({str(PROJECT)!r})",
    source,
    count=1,
    flags=re.MULTILINE,
)
if replacements != 1:
    raise RuntimeError("Could not replace the locked pipeline project-path assignment")
exec(compile(source, str(PIPELINE), "exec"), locked.__dict__)


def fit_residual(train_y, train_covariates, test_y, test_covariates):
    """Apply a training-only standardized least-squares residualization."""
    columns = []
    test_columns = []
    for train, test in zip(train_covariates, test_covariates):
        mean = float(np.mean(train))
        sd = float(np.std(train))
        if sd == 0:
            sd = 1.0
        columns.append((np.asarray(train, dtype=float) - mean) / sd)
        test_columns.append((np.asarray(test, dtype=float) - mean) / sd)
    train_design = np.column_stack([np.ones(len(train_y)), *columns])
    test_design = np.column_stack([np.ones(len(test_y)), *test_columns])
    beta, *_ = np.linalg.lstsq(train_design, train_y, rcond=None)
    return np.asarray(test_y, dtype=float) - test_design @ beta


def holdout_rows(cohort, sample, unit, matrix, genes, coords, array_rows, array_cols):
    scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
    ids = locked.block_ids(array_rows, array_cols)
    result = []
    for block in sorted(set(ids.tolist())):
        test = ids == block
        train = ~test
        if int(test.sum()) < MIN_TEST_SPOTS or int(train.sum()) < MIN_TEST_SPOTS:
            continue
        mreg_test = fit_residual(
            scores["mregDC_strict"][train],
            [scores[name][train] for name in ("library", "PTPRC", "epithelial", "stromal", "DC_core")],
            scores["mregDC_strict"][test],
            [scores[name][test] for name in ("library", "PTPRC", "epithelial", "stromal", "DC_core")],
        )
        tfh_test = fit_residual(
            scores["Tfh_like"][train],
            [scores[name][train] for name in ("library", "PTPRC", "epithelial", "stromal")],
            scores["Tfh_like"][test],
            [scores[name][test] for name in ("library", "PTPRC", "epithelial", "stromal")],
        )
        # The threshold is fitted in training spots; no held-out DC-core value
        # affects that decision other than its own comparison to the fixed cut.
        mask = scores["DC_core"][test] > float(np.mean(scores["DC_core"][train]))
        neighbours = locked.neighbour_index(coords[test])
        if int(mask.sum()) < MIN_MASKED_SPOTS or neighbours.shape[1] < 6:
            continue
        r = locked.local_r(mreg_test, tfh_test, neighbours, mask)
        if not np.isfinite(r):
            continue
        result.append({
            "cohort": cohort,
            "sample": sample,
            "unit": unit,
            "held_out_block": block,
            "test_spots": int(test.sum()),
            "test_dc_enriched_spots": int(mask.sum()),
            "training_spots": int(train.sum()),
            "held_out_local_r": float(r),
        })
    return result


def write_tsv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def bootstrap_median(values, rng, draws=2000):
    values = np.asarray(values, dtype=float)
    boot = np.array([np.median(rng.choice(values, len(values), replace=True)) for _ in range(draws)])
    return float(np.median(values)), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def patient_aggregate(rows):
    per_unit = defaultdict(list)
    for row in rows:
        per_unit[row["unit"]].append(row["held_out_local_r"])
    return [{"cohort": rows[0]["cohort"], "unit": unit, "holdout_block_median_r": float(np.median(values)), "n_eligible_blocks": len(values)} for unit, values in sorted(per_unit.items())]


def leave_one_out(rows, unit_label):
    values = {row[unit_label]: float(row["primary_mregDC_Tfh_local_r"]) for row in rows}
    output = []
    for omitted in sorted(values):
        retained = [value for key, value in values.items() if key != omitted]
        output.append({"omitted_unit": omitted, "n_retained": len(retained), "cohort_median_r_without_unit": float(np.median(retained))})
    return output


def run_g278():
    rows = []
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, array_rows, array_cols, patient = locked.read_g278(path)
        rows.extend(holdout_rows("GSE278687", sample, patient, matrix, genes, coords, array_rows, array_cols))
    patient = patient_aggregate(rows)
    primary = []
    with (locked.OUT / "GSE278687_v2_per_patient.tsv").open(encoding="utf-8") as handle:
        primary = list(csv.DictReader(handle, delimiter="\t"))
    loo = leave_one_out(primary, "patient")
    return rows, patient, loo


def run_g277():
    manifest = locked.g277_manifest()
    rows = []
    with tarfile.open(locked.G277_RAW, "r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith(".tar.gz"):
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in manifest or manifest[gsm]["tumor_program_scoreability_candidate"] != "1":
                continue
            loaded = locked.read_g277(gsm, outer.extractfile(member).read(), manifest[gsm])
            if loaded is None:
                continue
            sample, matrix, genes, coords, array_rows, array_cols, unit = loaded
            rows.extend(holdout_rows("GSE277116", sample, sample, matrix, genes, coords, array_rows, array_cols))
    sample = patient_aggregate(rows)
    primary = []
    with (locked.OUT / "GSE277116_v2_per_package.tsv").open(encoding="utf-8") as handle:
        primary = [row for row in csv.DictReader(handle, delimiter="\t") if row["analysis_set"] == "main_18"]
    loo = leave_one_out(primary, "sample")
    return rows, sample, loo


def main():
    rng = np.random.default_rng(SEED)
    g278_blocks, g278_patient, g278_loo = run_g278()
    g277_blocks, g277_sample, g277_loo = run_g277()
    for name, rows in {
        "GSE278687_block_holdout.tsv": g278_blocks,
        "GSE278687_patient_holdout_summary.tsv": g278_patient,
        "GSE278687_leave_one_patient_out.tsv": g278_loo,
        "GSE277116_block_holdout.tsv": g277_blocks,
        "GSE277116_sample_holdout_summary.tsv": g277_sample,
        "GSE277116_leave_one_sample_out.tsv": g277_loo,
    }.items():
        write_tsv(OUT / name, rows)
    g278_med, g278_lo, g278_hi = bootstrap_median([row["holdout_block_median_r"] for row in g278_patient], rng)
    g277_med, g277_lo, g277_hi = bootstrap_median([row["holdout_block_median_r"] for row in g277_sample], rng)
    summary = {
        "definition": "within-section contiguous 8x8 array-coordinate held-out block; training-only residualization and DC-core threshold; held-out-only six-neighbour fields",
        "not_a_new_primary_test": True,
        "GSE278687": {
            "sections_with_eligible_blocks": len(set(row["sample"] for row in g278_blocks)),
            "eligible_blocks": len(g278_blocks),
            "patients": len(g278_patient),
            "patient_median_of_holdout_block_medians_r": g278_med,
            "bootstrap_ci_low": g278_lo,
            "bootstrap_ci_high": g278_hi,
            "positive_patients": int(sum(row["holdout_block_median_r"] > 0 for row in g278_patient)),
            "leave_one_patient_out_range": [float(min(row["cohort_median_r_without_unit"] for row in g278_loo)), float(max(row["cohort_median_r_without_unit"] for row in g278_loo))],
        },
        "GSE277116": {
            "samples_with_eligible_blocks": len(set(row["sample"] for row in g277_blocks)),
            "eligible_blocks": len(g277_blocks),
            "samples": len(g277_sample),
            "sample_median_of_holdout_block_medians_r": g277_med,
            "bootstrap_ci_low": g277_lo,
            "bootstrap_ci_high": g277_hi,
            "positive_samples": int(sum(row["holdout_block_median_r"] > 0 for row in g277_sample)),
            "leave_one_sample_out_range": [float(min(row["cohort_median_r_without_unit"] for row in g277_loo)), float(max(row["cohort_median_r_without_unit"] for row in g277_loo))],
        },
    }
    (OUT / "holdout_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
