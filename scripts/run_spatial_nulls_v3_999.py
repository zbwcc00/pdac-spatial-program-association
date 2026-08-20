"""Extended spatial-null inference with 999 draws and mask-Moran audit.

This is a versioned follow-up to the locked v2 primary analysis. It preserves
the primary endpoint, all scoring, covariates, masking and aggregation rules.
Only the Monte Carlo resolution is increased and the graph-null audit is
expanded to report outcome Moran matching both globally and within the exact
DC-core mask used for the local correlation.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts/run_unified_primary_pipeline_v2.py"
OUT = PROJECT / "data/03_results/spatial_nulls_v3_999"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260813
N_DRAWS = int(os.environ.get("SPATIAL_NULL_DRAWS", "999"))
BATCH = os.environ.get("SPATIAL_NULL_BATCH", "full")


def load_locked_module():
    specification = importlib.util.spec_from_file_location("locked_v2_for_v3", PIPELINE)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not load pipeline: {PIPELINE}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    module.OUT = OUT
    return module


locked = load_locked_module()
LOCKED_OUT = PROJECT / "data/03_results/unified_primary_pipeline_v2"


def numericize(row):
    """Convert TSV numeric fields while leaving dataset identifiers intact."""
    output = {}
    for key, value in row.items():
        try:
            output[key] = float(value) if value not in {"", None} else value
        except ValueError:
            output[key] = value
    return output


def make_mask_weights(neighbours, mask):
    """Build the masked subgraph once per section, never once per surrogate."""
    index = np.flatnonzero(mask)
    if len(index) < 20:
        return index, None
    lookup = np.full(len(mask), -1, dtype=int)
    lookup[index] = np.arange(len(index))
    rows, cols = [], []
    for original in index:
        retained = lookup[neighbours[original]]
        retained = retained[retained >= 0]
        rows.extend([lookup[original]] * len(retained)); cols.extend(retained.tolist())
    if len(rows) < 2:
        return index, None
    from scipy import sparse
    adjacency = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index), len(index)))
    adjacency = ((adjacency + adjacency.T) > 0).astype(float).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel(); degree[degree == 0] = 1
    return index, sparse.diags(1 / degree) @ adjacency


def add_audit_graphs(states):
    for state in states:
        state["global_weights"] = locked.graph_weights(state["neighbours"])
        state["mask_index"], state["mask_weights"] = make_mask_weights(state["neighbours"], state["mask"])
        state["global_target"] = locked.moran(state["tfh"], state["global_weights"])
        state["mask_target"] = (
            locked.moran(state["tfh"][state["mask_index"]], state["mask_weights"])
            if state["mask_weights"] is not None else np.nan
        )


def graph_null_with_mask_audit(section_state, section_rows, rng, cohort, unit_label):
    """Generate one matched graph field per section, then aggregate per draw."""
    observed_by_unit = defaultdict(list)
    for row in section_rows:
        observed_by_unit[row["unit"]].append(row["primary_mregDC_Tfh_local_r"])
    observed = float(np.median([np.median(values) for values in observed_by_unit.values()]))
    draws = []
    global_deltas, mask_deltas = [], []
    for _ in range(N_DRAWS):
        by_unit = defaultdict(list)
        section_global, section_mask = [], []
        for row, state in zip(section_rows, section_state):
            surrogate, global_achieved, _ = locked.graph_surrogate(state["tfh"], state["global_weights"], state["global_target"], rng)
            mask_target = state["mask_target"]
            mask_achieved = (
                locked.moran(surrogate[state["mask_index"]], state["mask_weights"])
                if state["mask_weights"] is not None else np.nan
            )
            by_unit[row["unit"]].append(locked.local_r(state["mreg"], surrogate, state["neighbours"], state["mask"]))
            section_global.append(abs(global_achieved - state["global_target"]))
            if np.isfinite(mask_target) and np.isfinite(mask_achieved):
                section_mask.append(abs(mask_achieved - mask_target))
        draws.append(float(np.median([np.median(values) for values in by_unit.values()])))
        global_deltas.append(float(np.median(section_global)))
        mask_deltas.append(float(np.median(section_mask)) if section_mask else np.nan)
    draws = np.asarray(draws)
    extreme = int(np.sum(np.abs(draws) >= abs(observed)))
    return {
        "cohort": cohort,
        "statistic": f"cohort median of {unit_label}-median local r",
        "n_units": len(observed_by_unit),
        "n_null_draws": N_DRAWS,
        "observed": observed,
        "null_median": float(np.median(draws)),
        "null_ci_low": float(np.quantile(draws, .025)),
        "null_ci_high": float(np.quantile(draws, .975)),
        "n_as_or_more_extreme": extreme,
        "two_sided_monte_carlo_p": float((1 + extreme) / (N_DRAWS + 1)),
        "global_moran_abs_delta_per_draw_section_median": float(np.nanmedian(global_deltas)),
        "global_moran_abs_delta_per_draw_section_median_max": float(np.nanmax(global_deltas)),
        "mask_moran_abs_delta_per_draw_section_median": float(np.nanmedian(mask_deltas)),
        "mask_moran_abs_delta_per_draw_section_median_max": float(np.nanmax(mask_deltas)),
        "graph_null_target": "global outcome Moran I; mask Moran I is an audit, not separately calibrated",
    }, draws, global_deltas, mask_deltas


def write_tsv(path, rows):
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def run_g278(rng):
    sections, states = [], []
    with (LOCKED_OUT / "GSE278687_v2_per_section.tsv").open(encoding="utf-8") as handle:
        locked_rows = {row["sample"]: numericize(row) for row in csv.DictReader(handle, delimiter="\t")}
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, rows, cols, patient = locked.read_g278(path)
        scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
        mreg, tfh, mask, neighbours, _ = locked.prepared_fields(scores, coords)
        state = {"mreg": mreg, "tfh": tfh, "mask": mask, "neighbours": neighbours, "block_ids": locked.block_ids(rows, cols)}
        sections.append(locked_rows[sample]); states.append(state)
    patients = locked.g278_patient_rows(sections)
    add_audit_graphs(states)
    block = locked.patient_block_null(states, sections, rng, n=N_DRAWS)
    block["n_as_or_more_extreme"] = int(round(block["two_sided_p"] * (N_DRAWS + 1) - 1))
    block["two_sided_monte_carlo_p"] = block.pop("two_sided_p")
    graph, graph_draws, global_deltas, mask_deltas = graph_null_with_mask_audit(states, sections, rng, "GSE278687", "patient")
    write_tsv(OUT / "GSE278687_v3_per_section.tsv", sections)
    write_tsv(OUT / "GSE278687_v3_per_patient.tsv", patients)
    write_tsv(OUT / "GSE278687_v3_patient_block_null.tsv", [block])
    write_tsv(OUT / "GSE278687_v3_patient_graph_null.tsv", [graph])
    suffix = f"_batch{BATCH}" if BATCH != "full" else ""
    np.savetxt(OUT / f"GSE278687_v3_patient_graph_null_draws{suffix}.tsv", graph_draws, fmt="%.12g", header="cohort_median_patient_r", comments="")
    write_tsv(OUT / f"GSE278687_v3_patient_graph_null_audit{suffix}.tsv", [
        {"draw": index + 1, "global_moran_section_median_abs_delta": global_delta, "mask_moran_section_median_abs_delta": mask_delta}
        for index, (global_delta, mask_delta) in enumerate(zip(global_deltas, mask_deltas))
    ])
    return sections, patients, block, graph


def run_g277(rng):
    manifest = locked.g277_manifest()
    rows, main_rows, states = [], [], []
    with (LOCKED_OUT / "GSE277116_v2_per_package.tsv").open(encoding="utf-8") as handle:
        locked_rows = {row["gsm"]: numericize(row) for row in csv.DictReader(handle, delimiter="\t")}
    with tarfile.open(locked.G277_RAW, "r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith(".tar.gz"):
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in manifest:
                continue
            loaded = locked.read_g277(gsm, outer.extractfile(member).read(), manifest[gsm])
            if loaded is None:
                continue
            sample, matrix, genes, coords, array_rows, array_cols, unit = loaded
            result = locked_rows[gsm]
            scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
            mreg, tfh, mask, neighbours, _ = locked.prepared_fields(scores, coords)
            state = {"mreg": mreg, "tfh": tfh, "mask": mask, "neighbours": neighbours, "block_ids": locked.block_ids(array_rows, array_cols)}
            rows.append(result)
            if result["analysis_set"] == "main_18":
                main_rows.append(result); states.append(state)
    rows.sort(key=lambda row: row["gsm"])
    add_audit_graphs(states)
    graph, graph_draws, global_deltas, mask_deltas = graph_null_with_mask_audit(states, main_rows, rng, "GSE277116", "sample")
    write_tsv(OUT / "GSE277116_v3_per_package.tsv", rows)
    write_tsv(OUT / "GSE277116_v3_main18_graph_null.tsv", [graph])
    suffix = f"_batch{BATCH}" if BATCH != "full" else ""
    np.savetxt(OUT / f"GSE277116_v3_main18_graph_null_draws{suffix}.tsv", graph_draws, fmt="%.12g", header="cohort_median_sample_r", comments="")
    write_tsv(OUT / f"GSE277116_v3_main18_graph_null_audit{suffix}.tsv", [
        {"draw": index + 1, "global_moran_section_median_abs_delta": global_delta, "mask_moran_section_median_abs_delta": mask_delta}
        for index, (global_delta, mask_delta) in enumerate(zip(global_deltas, mask_deltas))
    ])
    return rows, main_rows, graph


def main():
    # Full runs retain the locked seed. Batched runs use disjoint deterministic
    # streams so that concatenated batches are independent Monte Carlo draws.
    batch_offset = int(BATCH) * 100_000 if BATCH.isdigit() else 0
    rng = np.random.default_rng(SEED + batch_offset)
    target = sys.argv[1].lower() if len(sys.argv) == 2 else "all"
    if target not in {"all", "gse278687", "gse277116"}:
        raise SystemExit("Usage: run_spatial_nulls_v3_999.py [all|gse278687|gse277116]")
    block = g278_graph = g277_graph = None
    if target in {"all", "gse278687"}:
        print("Starting GSE278687 extended nulls", flush=True)
        _, _, block, g278_graph = run_g278(rng)
        print("Completed GSE278687 extended nulls", flush=True)
    if target in {"all", "gse277116"}:
        print("Starting GSE277116 extended nulls", flush=True)
        _, _, g277_graph = run_g277(rng)
        print("Completed GSE277116 extended nulls", flush=True)
    summary = {
        "version": "v3_999", "seed": SEED + batch_offset, "batch": BATCH, "n_null_draws": N_DRAWS,
        "endpoint": "unchanged v2 composition-adjusted focal-spot-excluded six-neighbour local mregDC-like--Tfh-like correlation",
        "GSE278687_patient_block_null": block,
        "GSE278687_patient_graph_null": g278_graph,
        "GSE277116_sample_graph_null": g277_graph,
        "mask_moran_note": "Graph surrogates are calibrated to global outcome Moran I; within-mask Moran matching is reported as an audit rather than claimed as calibrated.",
    }
    (OUT / "v3_999_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
