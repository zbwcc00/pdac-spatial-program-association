"""Mask-targeted Moran spatial-null sensitivity for the locked PDAC endpoint.

This v4 sensitivity keeps the v2 scores, covariates, focal-spot-excluded
six-neighbour fields, DC-core mask, and unit-level aggregation unchanged.  In
contrast to v3, each rank-preserving graph surrogate is calibrated to Moran's
I *within the exact DC-core mask where the local-field statistic is evaluated*.
Whole-section Moran's I is retained as an audit only.  The two targets cannot
generally be matched simultaneously by this one-parameter graph filter, so
v4 must never be interpreted as confirmation of the v3 global-Moran null.
"""
from __future__ import annotations

import csv
import importlib.util
import io
import gzip
import json
import os
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts" / "run_unified_primary_pipeline_v2.py"
LOCKED_OUT = PROJECT / "data" / "03_results" / "unified_primary_pipeline_v2"
OUT = PROJECT / "data" / "03_results" / "spatial_nulls_v4_mask_moran"
OUT.mkdir(parents=True, exist_ok=True)

BASE_SEED = 20260813
BATCH = os.environ.get("MASK_MORAN_BATCH", "full")
N_DRAWS = int(os.environ.get("MASK_MORAN_DRAWS", "999"))


def load_locked_module():
    specification = importlib.util.spec_from_file_location("locked_v2_for_v4", PIPELINE)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not load pipeline: {PIPELINE}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


locked = load_locked_module()


def numericize(row):
    result = {}
    for key, value in row.items():
        try:
            result[key] = float(value) if value not in {"", None} else value
        except ValueError:
            result[key] = value
    return result


def write_tsv(path, rows):
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def mask_weights(neighbours, mask):
    """Build a symmetric row-standardized kNN graph restricted to masked spots."""
    index = np.flatnonzero(mask)
    if len(index) < 20:
        return index, None
    remap = np.full(len(mask), -1, dtype=int)
    remap[index] = np.arange(len(index))
    rows, cols = [], []
    for old_index in index:
        retained = remap[neighbours[old_index]]
        retained = retained[retained >= 0]
        rows.extend([remap[old_index]] * len(retained))
        cols.extend(retained.tolist())
    if len(rows) < 2:
        return index, None
    from scipy import sparse
    adjacency = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index), len(index)))
    adjacency = ((adjacency + adjacency.T) > 0).astype(float).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    degree[degree == 0] = 1
    return index, sparse.diags(1 / degree) @ adjacency


def graph_surrogate_mask_target(values, global_weights, mask_index, mask_graph, target, rng):
    """Generate one rank-preserved field calibrated to masked-subgraph Moran I.

    The candidate field is generated over the whole section, then its masked
    values are evaluated on the induced mask graph.  This preserves the local
    field definition, including any mask-to-unmasked neighbour links, while
    making mask Moran the sole calibration target.
    """
    raw = rng.normal(size=len(values))
    states = [raw]
    for _ in range(locked.GRAPH_MAX_STEPS):
        raw = 0.5 * raw + 0.5 * np.asarray(global_weights @ raw).ravel()
        states.append(raw)

    def evaluate(state):
        candidate = locked.rank_map(state, values)
        return candidate, locked.moran(candidate[mask_index], mask_graph)

    candidates = [evaluate(state) for state in states]
    morans = [item[1] for item in candidates]
    best = int(np.nanargmin(np.abs(np.asarray(morans) - target)))
    candidate_best, achieved = candidates[best]
    error_best = abs(achieved - target)
    intervals = [
        (i, i + 1) for i in range(len(morans) - 1)
        if np.isfinite(morans[i]) and np.isfinite(morans[i + 1])
        and (morans[i] - target) * (morans[i + 1] - target) <= 0
    ]
    if intervals:
        start, end = min(intervals, key=lambda pair: abs(morans[pair[0]] - target) + abs(morans[pair[1]] - target))
        low, high = 0.0, 1.0
        low_moran, high_moran = morans[start], morans[end]
        for _ in range(locked.GRAPH_MATCH_ITERATIONS):
            middle = (low + high) / 2
            candidate, candidate_moran = evaluate((1 - middle) * states[start] + middle * states[end])
            candidate_error = abs(candidate_moran - target)
            if candidate_error < error_best:
                candidate_best, achieved, error_best = candidate, candidate_moran, candidate_error
            if (low_moran - target) * (candidate_moran - target) <= 0:
                high, high_moran = middle, candidate_moran
            else:
                low, low_moran = middle, candidate_moran
    return candidate_best, achieved


def aggregate_observed(rows):
    per_unit = defaultdict(list)
    for row in rows:
        per_unit[row["unit"]].append(row["primary_mregDC_Tfh_local_r"])
    return float(np.median([np.median(values) for values in per_unit.values()])), per_unit


def prepare_state(matrix, genes, coords):
    scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
    mreg, tfh, mask, neighbours, _ = locked.prepared_fields(scores, coords)
    global_graph = locked.graph_weights(neighbours)
    index, induced_graph = mask_weights(neighbours, mask)
    if induced_graph is None:
        raise RuntimeError("DC-core mask cannot support a Moran sensitivity graph")
    return {
        "mreg": mreg, "tfh": tfh, "mask": mask, "neighbours": neighbours,
        "global_graph": global_graph, "mask_index": index, "mask_graph": induced_graph,
        "mask_target": locked.moran(tfh[index], induced_graph),
        "global_target": locked.moran(tfh, global_graph),
    }


def cohort_null(rows, states, cohort, label, rng):
    observed, per_unit = aggregate_observed(rows)
    draws, mask_deltas, global_deltas = [], [], []
    for _ in range(N_DRAWS):
        values_by_unit = defaultdict(list)
        draw_mask_delta, draw_global_delta = [], []
        for row, state in zip(rows, states):
            surrogate, mask_achieved = graph_surrogate_mask_target(
                state["tfh"], state["global_graph"], state["mask_index"], state["mask_graph"], state["mask_target"], rng,
            )
            values_by_unit[row["unit"]].append(locked.local_r(state["mreg"], surrogate, state["neighbours"], state["mask"]))
            draw_mask_delta.append(abs(mask_achieved - state["mask_target"]))
            global_achieved = locked.moran(surrogate, state["global_graph"])
            draw_global_delta.append(abs(global_achieved - state["global_target"]))
        draws.append(float(np.median([np.median(values) for values in values_by_unit.values()])))
        mask_deltas.append(float(np.median(draw_mask_delta)))
        global_deltas.append(float(np.median(draw_global_delta)))
    draws = np.asarray(draws)
    extreme = int(np.sum(np.abs(draws) >= abs(observed)))
    result = {
        "cohort": cohort,
        "statistic": f"cohort median of {label}-median local r",
        "n_units": len(per_unit),
        "n_null_draws": N_DRAWS,
        "observed": observed,
        "null_median": float(np.median(draws)),
        "null_ci_low": float(np.quantile(draws, 0.025)),
        "null_ci_high": float(np.quantile(draws, 0.975)),
        "n_as_or_more_extreme": extreme,
        "two_sided_add_one_monte_carlo_p": float((1 + extreme) / (N_DRAWS + 1)),
        "mask_moran_abs_delta_draw_section_median": float(np.median(mask_deltas)),
        "mask_moran_abs_delta_draw_section_median_max": float(np.max(mask_deltas)),
        "global_moran_abs_delta_draw_section_median": float(np.median(global_deltas)),
        "global_moran_abs_delta_draw_section_median_max": float(np.max(global_deltas)),
        "calibration_target": "DC-core-mask induced-subgraph outcome Moran I; whole-section Moran I is audit only",
    }
    return result, draws, mask_deltas, global_deltas


def gse278687():
    rows_by_sample = {}
    with (LOCKED_OUT / "GSE278687_v2_per_section.tsv").open(encoding="utf-8") as handle:
        rows_by_sample = {row["sample"]: numericize(row) for row in csv.DictReader(handle, delimiter="\t")}
    rows, states = [], []
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, _, _, _ = locked.read_g278(path)
        rows.append(rows_by_sample[sample])
        states.append(prepare_state(matrix, genes, coords))
    return rows, states


def gse277116():
    manifest = locked.g277_manifest()
    with (LOCKED_OUT / "GSE277116_v2_per_package.tsv").open(encoding="utf-8") as handle:
        locked_rows = {row["gsm"]: numericize(row) for row in csv.DictReader(handle, delimiter="\t")}
    pairs = []
    with tarfile.open(locked.G277_RAW, "r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith(".tar.gz"):
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in manifest or locked_rows[gsm]["analysis_set"] != "main_18":
                continue
            loaded = locked.read_g277(gsm, outer.extractfile(member).read(), manifest[gsm])
            if loaded is None:
                raise RuntimeError(f"Locked main package became unreadable: {gsm}")
            _, matrix, genes, coords, _, _, _ = loaded
            pairs.append((locked_rows[gsm], prepare_state(matrix, genes, coords)))
    pairs.sort(key=lambda pair: pair[0]["gsm"])
    rows = [pair[0] for pair in pairs]
    states = [pair[1] for pair in pairs]
    return rows, states


def main():
    if len(sys.argv) != 2 or sys.argv[1].lower() not in {"gse278687", "gse277116"}:
        raise SystemExit("Usage: run_mask_moran_sensitivity_v4.py [gse278687|gse277116]")
    target = sys.argv[1].lower()
    batch_offset = int(BATCH) * 100_000 if BATCH.isdigit() else 0
    seed = BASE_SEED + batch_offset
    rng = np.random.default_rng(seed)
    rows, states = gse278687() if target == "gse278687" else gse277116()
    cohort = "GSE278687" if target == "gse278687" else "GSE277116"
    label = "patient" if target == "gse278687" else "sample"
    result, draws, mask_deltas, global_deltas = cohort_null(rows, states, cohort, label, rng)
    suffix = f"_batch{BATCH}" if BATCH != "full" else ""
    prefix = f"{cohort}_v4_{label}_mask_moran"
    write_tsv(OUT / f"{prefix}_null{suffix}.tsv", [result])
    np.savetxt(OUT / f"{prefix}_draws{suffix}.tsv", draws, fmt="%.12g", header="cohort_median_unit_r", comments="")
    write_tsv(OUT / f"{prefix}_audit{suffix}.tsv", [
        {"draw": i + 1, "mask_moran_section_median_abs_delta": mask_delta, "global_moran_section_median_abs_delta": global_delta}
        for i, (mask_delta, global_delta) in enumerate(zip(mask_deltas, global_deltas))
    ])
    manifest = {
        "version": "v4_mask_moran_sensitivity", "cohort": cohort, "batch": BATCH, "seed": seed,
        "draws": N_DRAWS, "endpoint": "locked v2 primary endpoint", "result": result,
        "caveat": "Only mask-subgraph Moran I is calibrated. Global Moran mismatch is an audit. This is a sensitivity null, not a simultaneous global-and-mask matching method.",
    }
    (OUT / f"{prefix}_manifest{suffix}.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
