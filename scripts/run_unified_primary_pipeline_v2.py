"""Versioned, unified spatial inference for PDAC immune-program organization.

The primary estimand is the patient- (GSE278687) or sample-level (GSE277116)
correlation of six-neighbour local fields.  The focal spot is excluded from
each field.  This version repairs two inferential problems in v1:

1. A block-constrained null is aggregated at the patient level before a cohort
   P value is calculated; section-level P values are never pooled.
2. A graph-diffusion spatial-surrogate null is added.  It rank-preserves the
   observed outcome values and calibrates a graph-smoothed surrogate to the
   observed global Moran statistic.  It is a sensitivity analysis, not a
   proof of cell-cell interaction.

The GSE277116 analysis uses the same scores, covariates, mask, neighbour rule,
and null procedures.  It remains an external sample-level technical
replication because GEO does not supply patient identifiers for that series.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import tarfile
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import scipy.io
from scipy import sparse, stats
from scipy.spatial import cKDTree

PROJECT = Path(__file__).resolve().parents[1]
FROZEN = json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8"))
G278 = PROJECT / "data/01_unpacked/spatial/GSE278687"
G278_COORDS = PROJECT / "data/01_unpacked/spatial/GSE278687_spatial"
G277_RAW = PROJECT / "data/00_raw/spatial/GSE277116_RAW.tar"
G277_QC = PROJECT / "data/03_results/GSE277116_full_package_qc/GSE277116_tumor_replication_manifest.tsv"
OUT = PROJECT / "data/03_results/unified_primary_pipeline_v2"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260813
K = 6
BLOCK = 8
N_BLOCK_NULL = 199
N_GRAPH_NULL = 199
N_COHORT_GRAPH_NULL = 199
GRAPH_MAX_STEPS = 12
GRAPH_MATCH_ITERATIONS = 12

EPITHELIAL = ("KRT8", "KRT18", "KRT19", "EPCAM", "KRT7", "MUC1", "CEACAM6", "KRT17")
STROMAL = ("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2", "SPARC", "TAGLN")

# None of these reference programs shares a gene with the frozen Tfh-like
# program. They are adjustment/competition references, not asserted cell calls.
COMPETITORS = {
    "broad_T_nonoverlap": ("TRAC", "TRBC1", "TRBC2", "CD247", "LCK", "CD2", "CD7"),
    "non_Tfh_CD4_nonoverlap": ("LTB", "LEF1", "TCF7", "SELL", "MAL"),
    "Treg_nonoverlap": ("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF4"),
    "exhausted_CD8_nonoverlap": ("CD8A", "CD8B", "LAG3", "HAVCR2", "ENTPD1", "LAYN"),
    "cDC_nonoverlap": ("CLEC10A", "CD1C", "FCER1A", "CD1E"),
    "macrophage_nonoverlap": ("C1QA", "C1QB", "C1QC", "APOC1", "GPNMB", "SPP1"),
}


def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / sd if sd > 0 else np.zeros_like(x)


def residual(y, *covariates):
    design = np.column_stack([np.ones(len(y)), *[z(v) for v in covariates]])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return np.asarray(y, dtype=float) - design @ beta


def correlation(left, right):
    keep = np.isfinite(left) & np.isfinite(right)
    if keep.sum() < 20 or np.nanstd(np.asarray(left)[keep]) == 0 or np.nanstd(np.asarray(right)[keep]) == 0:
        return np.nan
    return float(np.corrcoef(np.asarray(left)[keep], np.asarray(right)[keep])[0, 1])


def neighbour_index(coords, k=K):
    q = min(k + 1, len(coords))
    _, idx = cKDTree(coords).query(coords, k=q)
    return idx[:, 1:] if idx.ndim == 2 and idx.shape[1] > 1 else np.empty((len(coords), 0), dtype=int)


def local_field(values, neighbours):
    return np.asarray(values)[neighbours].mean(axis=1) if neighbours.shape[1] else np.zeros(len(values))


def graph_weights(neighbours):
    n = len(neighbours)
    if neighbours.shape[1] == 0:
        return sparse.eye(n, format="csr")
    rows = np.repeat(np.arange(n), neighbours.shape[1])
    cols = neighbours.reshape(-1)
    adjacency = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    adjacency = ((adjacency + adjacency.T) > 0).astype(float).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    degree[degree == 0] = 1
    return sparse.diags(1 / degree) @ adjacency


def moran(values, weights):
    centered = z(values)
    denominator = float(centered @ centered)
    return float(centered @ (weights @ centered) / denominator) if denominator > 0 else np.nan


def score(matrix, gene_index, genes):
    idx = [gene_index[g] for g in genes if g in gene_index]
    return np.log1p(np.asarray(matrix[idx, :].sum(axis=0)).ravel()) if idx else np.zeros(matrix.shape[1])


def make_scores(matrix, genes):
    gene_index = {g.upper(): i for i, g in enumerate(genes)}
    programs = {name: score(matrix, gene_index, meta["genes"]) for name, meta in FROZEN.items() if name in {"DC_core", "mregDC_strict", "Tfh_like"}}
    programs["library"] = np.log1p(np.asarray(matrix.sum(axis=0)).ravel())
    programs["PTPRC"] = score(matrix, gene_index, ("PTPRC",))
    programs["epithelial"] = score(matrix, gene_index, EPITHELIAL)
    programs["stromal"] = score(matrix, gene_index, STROMAL)
    for name, geneset in COMPETITORS.items():
        programs[name] = score(matrix, gene_index, geneset)
    return programs


def prepared_fields(scores, coords):
    neighbours = neighbour_index(coords)
    mreg = residual(scores["mregDC_strict"], scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"], scores["DC_core"])
    tfh = residual(scores["Tfh_like"], scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"])
    broad = residual(scores["broad_T_nonoverlap"], scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"])
    tfh_after_broad = residual(tfh, broad)
    mreg_no_ccl19 = residual(
        score_from_existing(scores, "mregDC_no_CCL19"), scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"], scores["DC_core"]
    )
    mask = z(scores["DC_core"]) > 0
    local_mreg = local_field(z(mreg), neighbours)
    local_tfh = local_field(z(tfh), neighbours)
    fields = {
        "primary_mregDC_Tfh_local_r": correlation(local_mreg[mask], local_tfh[mask]),
        "mregDC_Tfh_after_nonoverlap_broad_T_adjustment_r": correlation(local_mreg[mask], local_field(z(tfh_after_broad), neighbours)[mask]),
        "mregDC_no_CCL19_Tfh_local_r": correlation(local_field(z(mreg_no_ccl19), neighbours)[mask], local_tfh[mask]),
    }
    local_competitors = {}
    for name in COMPETITORS:
        value = residual(scores[name], scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"])
        local_competitors[name] = local_field(z(value), neighbours)
        fields[f"mregDC_{name}_competitive_r"] = correlation(local_mreg[mask], local_competitors[name][mask])
    # Joint model: the Tfh coefficient is interpreted conditionally on the
    # non-overlapping reference fields; it is not a unique cell-state estimate.
    valid = mask & np.isfinite(local_mreg) & np.isfinite(local_tfh)
    design = np.column_stack([np.ones(valid.sum()), z(local_tfh[valid]), *[z(local_competitors[name][valid]) for name in COMPETITORS]])
    beta, *_ = np.linalg.lstsq(design, z(local_mreg[valid]), rcond=None)
    fields["joint_model_Tfh_standardized_beta"] = float(beta[1])
    return mreg, tfh, mask, neighbours, fields


def score_from_existing(scores, name):
    # Stored directly when scores are created; this helper is intentionally
    # small to make the CCL19 deletion definition visible at the call site.
    return scores[name]


def add_mreg_no_ccl19(scores, matrix, genes):
    gene_index = {g.upper(): i for i, g in enumerate(genes)}
    scores["mregDC_no_CCL19"] = score(matrix, gene_index, ("LAMP3", "FSCN1", "CCR7", "CD80", "CD86"))
    return scores


def local_r(left, right, neighbours, mask):
    return correlation(local_field(z(left), neighbours)[mask], local_field(z(right), neighbours)[mask])


def block_ids(array_rows, array_cols):
    return np.array([f"{int(row // BLOCK)}:{int(col // BLOCK)}" for row, col in zip(array_rows, array_cols)], dtype=str)


def block_surrogate(values, ids, rng):
    output = np.asarray(values, dtype=float).copy()
    for block in np.unique(ids):
        index = np.where(ids == block)[0]
        if len(index) >= 2:
            output[index] = rng.permutation(output[index])
    return output


def rank_map(values, template):
    ranks = np.argsort(np.argsort(values, kind="mergesort"), kind="mergesort")
    return np.sort(np.asarray(template, dtype=float))[ranks]


def graph_surrogate(values, weights, target_moran, rng):
    """Rank-preserving graph-filtered null matched to the observed Moran I.

    For every draw, an independent Gaussian field is smoothed on the kNN graph.
    Adjacent graph-filter states are interpolated continuously; bisection picks
    the interpolation weight whose rank-mapped field is closest to the target
    Moran I.  This avoids using a single discrete smoothing step, which can
    fail badly when its autocorrelation jumps over the observed value.
    """
    raw = rng.normal(size=weights.shape[0])
    states = [raw]
    for _ in range(GRAPH_MAX_STEPS):
        raw = 0.5 * raw + 0.5 * np.asarray(weights @ raw).ravel()
        states.append(raw)
    candidates = [rank_map(state, values) for state in states]
    values_moran = [moran(candidate, weights) for candidate in candidates]
    best = int(np.argmin(np.abs(np.asarray(values_moran) - target_moran)))
    candidate_best = candidates[best]
    error_best = abs(values_moran[best] - target_moran)
    # Find a neighbouring pair that brackets the target, if one exists.
    intervals = [(i, i + 1) for i in range(len(values_moran) - 1) if (values_moran[i] - target_moran) * (values_moran[i + 1] - target_moran) <= 0]
    if intervals:
        start, end = min(intervals, key=lambda pair: abs(values_moran[pair[0]] - target_moran) + abs(values_moran[pair[1]] - target_moran))
        low, high = 0.0, 1.0
        low_moran, high_moran = values_moran[start], values_moran[end]
        for _ in range(GRAPH_MATCH_ITERATIONS):
            middle = (low + high) / 2
            candidate = rank_map((1 - middle) * states[start] + middle * states[end], values)
            candidate_moran = moran(candidate, weights)
            candidate_error = abs(candidate_moran - target_moran)
            if candidate_error < error_best:
                candidate_best, error_best = candidate, candidate_error
            if (low_moran - target_moran) * (candidate_moran - target_moran) <= 0:
                high, high_moran = middle, candidate_moran
            else:
                low, low_moran = middle, candidate_moran
    return candidate_best, float(moran(candidate_best, weights)), float(error_best)


def graph_null(left, right, neighbours, mask, rng, n=N_GRAPH_NULL):
    weights = graph_weights(neighbours)
    target = moran(right, weights)
    observed = local_r(left, right, neighbours, mask)
    null, achieved, deltas = [], [], []
    for _ in range(n):
        surrogate, achieved_moran, delta = graph_surrogate(right, weights, target, rng)
        null.append(local_r(left, surrogate, neighbours, mask))
        achieved.append(achieved_moran)
        deltas.append(delta)
    null = np.asarray(null, dtype=float)
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (1 + len(null))
    return {"graph_null_p": float(p), "graph_null_median_r": float(np.nanmedian(null)), "outcome_moran": float(target), "surrogate_moran_median": float(np.nanmedian(achieved)), "surrogate_moran_abs_delta_median": float(np.nanmedian(deltas))}


def read_positions(path):
    result = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("in_tissue") == "1":
                result[row["barcode"]] = (float(row["pxl_col_in_fullres"]), float(row["pxl_row_in_fullres"]), int(row.get("array_row", 0)), int(row.get("array_col", 0)))
    return result


def read_g278(path):
    sample = path.name.replace("_filtered_feature_bc_matrix.h5", "")
    positions = read_positions(G278_COORDS / sample / "spatial/tissue_positions.csv")
    with h5py.File(path, "r") as handle:
        matrix_group = handle["matrix"]
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in matrix_group["barcodes"][:]]
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in matrix_group["features"]["name"][:]]
        matrix = sparse.csc_matrix((matrix_group["data"][:], matrix_group["indices"][:], matrix_group["indptr"][:]), shape=tuple(matrix_group["shape"][:])).tocsr()
    keep = np.array([i for i, barcode in enumerate(barcodes) if barcode in positions], dtype=int)
    meta = np.array([positions[barcodes[i]] for i in keep], dtype=float)
    return sample, matrix[:, keep], genes, meta[:, :2], meta[:, 2].astype(int), meta[:, 3].astype(int), sample.rsplit("_", 1)[-1].split("-")[0]


def parse_coords(raw, filename):
    text = gzip.decompress(raw).decode("utf-8", errors="replace") if filename.endswith(".gz") else raw.decode("utf-8", errors="replace")
    output = {}
    for line in text.splitlines()[1:]:
        fields = line.replace("\t", ",").split(",")
        if len(fields) >= 6:
            try:
                output[fields[0].strip()] = (int(float(fields[1])), int(float(fields[2])), int(float(fields[3])))
            except ValueError:
                pass
    return output


def member_by_basename(members, name):
    return next((member for member in members if Path(member.name).name == name), None)


def g277_manifest():
    with G277_QC.open(encoding="utf-8") as handle:
        return {row["gsm"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_g277(gsm, raw_nested, meta):
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw_nested)), mode="r") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        matrix_member = member_by_basename(members, "matrix.mtx.gz")
        barcode_member = member_by_basename(members, "barcodes.tsv.gz")
        feature_member = member_by_basename(members, "features.tsv.gz")
        coordinate_member = member_by_basename(members, "tissue_positions.csv") or member_by_basename(members, "tissue_positions_list.csv")
        if not all((matrix_member, barcode_member, feature_member, coordinate_member)):
            return None
        barcodes = gzip.decompress(archive.extractfile(barcode_member).read()).decode("utf-8", errors="replace").splitlines()
        feature_lines = gzip.decompress(archive.extractfile(feature_member).read()).decode("utf-8", errors="replace").splitlines()
        genes = [(line.split("\t")[1] if len(line.split("\t")) > 1 and line.split("\t")[1] else line.split("\t")[0]) for line in feature_lines]
        coordinates = parse_coords(archive.extractfile(coordinate_member).read(), coordinate_member.name)
        keep = np.array([i for i, barcode in enumerate(barcodes) if barcode in coordinates and coordinates[barcode][0] == 1], dtype=int)
        if len(keep) < 20:
            return None
        matrix = scipy.io.mmread(io.BytesIO(gzip.decompress(archive.extractfile(matrix_member).read()))).tocsr()[:, keep]
        array = np.array([coordinates[barcodes[i]] for i in keep], dtype=float)
    return meta["title"], matrix, genes, array[:, 1:3], array[:, 1].astype(int), array[:, 2].astype(int), gsm


def section_analysis(cohort, sample, matrix, genes, coords, rows, cols, unit, rng):
    scores = add_mreg_no_ccl19(make_scores(matrix, genes), matrix, genes)
    mreg, tfh, mask, neighbours, effects = prepared_fields(scores, coords)
    ids = block_ids(rows, cols)
    graph = graph_null(mreg, tfh, neighbours, mask, rng)
    result = {"cohort": cohort, "sample": sample, "unit": unit, "spots": len(coords), "dc_enriched_spots": int(mask.sum()), **effects, **graph}
    result["block_null_eligible"] = int(len(set(ids.tolist())) >= 3)
    return result, {"mreg": mreg, "tfh": tfh, "mask": mask, "neighbours": neighbours, "block_ids": ids}


def write_tsv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def median_ci(values, rng):
    values = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if not len(values):
        return np.nan, np.nan, np.nan, 0, np.nan
    boot = [np.median(rng.choice(values, len(values), replace=True)) for _ in range(2000)]
    return float(np.median(values)), float(np.quantile(boot, .025)), float(np.quantile(boot, .975)), int((values > 0).sum()), float(stats.binomtest(int((values > 0).sum()), len(values), .5).pvalue)


def summarize_sections(rows, unit_label, rng):
    endpoints = [key for key in rows[0] if (key.endswith("_r") and key != "graph_null_median_r") or key.endswith("_beta")]
    output = []
    for endpoint in endpoints:
        median, lo, hi, positive, sign_p = median_ci([row[endpoint] for row in rows], rng)
        output.append({"cohort": rows[0]["cohort"], "analysis_unit": unit_label, "endpoint": endpoint, "n_units": len(rows), "median": median, "bootstrap_ci_low": lo, "bootstrap_ci_high": hi, "positive": positive, "sign_test_p": sign_p})
    return output


def g278_patient_rows(section_rows):
    patients = sorted({row["unit"] for row in section_rows})
    endpoints = [key for key in section_rows[0] if (key.endswith("_r") and key != "graph_null_median_r") or key.endswith("_beta")]
    output = []
    for patient in patients:
        sections = [row for row in section_rows if row["unit"] == patient]
        output.append({"cohort": "GSE278687", "patient": patient, **{endpoint: float(np.median([row[endpoint] for row in sections])) for endpoint in endpoints}})
    return output


def patient_block_null(section_state, section_rows, rng, n=N_BLOCK_NULL):
    """Synchronously apply a section-level block null and aggregate per draw.

    This returns one cohort null distribution. It deliberately does not take a
    median of section P values, which has no patient-level null interpretation.
    """
    observed_patient = {row["unit"]: [] for row in section_rows}
    for row in section_rows:
        observed_patient[row["unit"]].append(row["primary_mregDC_Tfh_local_r"])
    observed = float(np.median([np.median(values) for values in observed_patient.values()]))
    draws = []
    for _ in range(n):
        per_patient = defaultdict(list)
        for row, state in zip(section_rows, section_state):
            right = block_surrogate(state["tfh"], state["block_ids"], rng)
            per_patient[row["unit"]].append(local_r(state["mreg"], right, state["neighbours"], state["mask"]))
        draws.append(float(np.median([np.median(values) for values in per_patient.values()])))
    draws = np.asarray(draws)
    p = float((1 + np.sum(np.abs(draws) >= abs(observed))) / (len(draws) + 1))
    return {"cohort": "GSE278687", "statistic": "cohort median of patient-median local r", "n_patients": len(observed_patient), "n_null_draws": n, "observed": observed, "null_median": float(np.median(draws)), "null_ci_low": float(np.quantile(draws, .025)), "null_ci_high": float(np.quantile(draws, .975)), "two_sided_p": p}


def aggregate_graph_null(section_state, section_rows, rng, n=N_COHORT_GRAPH_NULL, cohort="GSE278687", unit_label="patient"):
    """Global null from independently Moran-matched fields, aggregated by unit."""
    observed_patient = {row["unit"]: [] for row in section_rows}
    for row in section_rows:
        observed_patient[row["unit"]].append(row["primary_mregDC_Tfh_local_r"])
    observed = float(np.median([np.median(values) for values in observed_patient.values()]))
    draws, matched_errors = [], []
    for _ in range(n):
        per_patient = defaultdict(list)
        errors = []
        for row, state in zip(section_rows, section_state):
            weights = graph_weights(state["neighbours"])
            surrogate, _, error = graph_surrogate(state["tfh"], weights, moran(state["tfh"], weights), rng)
            per_patient[row["unit"]].append(local_r(state["mreg"], surrogate, state["neighbours"], state["mask"]))
            errors.append(error)
        draws.append(float(np.median([np.median(values) for values in per_patient.values()])))
        matched_errors.append(float(np.median(errors)))
    draws = np.asarray(draws)
    p = float((1 + np.sum(np.abs(draws) >= abs(observed))) / (len(draws) + 1))
    return {"cohort": cohort, "statistic": f"cohort median of {unit_label}-median local r", "n_units": len(observed_patient), "n_null_draws": n, "observed": observed, "null_median": float(np.median(draws)), "null_ci_low": float(np.quantile(draws, .025)), "null_ci_high": float(np.quantile(draws, .975)), "two_sided_p": p, "per_draw_section_median_abs_moran_delta_median": float(np.median(matched_errors)), "per_draw_section_median_abs_moran_delta_max": float(np.max(matched_errors))}


def run_g278(rng):
    sections, states = [], []
    for path in sorted(G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, rows, cols, patient = read_g278(path)
        result, state = section_analysis("GSE278687", sample, matrix, genes, coords, rows, cols, patient, rng)
        sections.append(result)
        states.append(state)
    if len(sections) != 21:
        raise RuntimeError(f"Expected 21 GSE278687 sections, found {len(sections)}")
    patients = g278_patient_rows(sections)
    block = patient_block_null(states, sections, rng)
    graph = aggregate_graph_null(states, sections, rng, cohort="GSE278687", unit_label="patient")
    write_tsv(OUT / "GSE278687_v2_per_section.tsv", sections)
    write_tsv(OUT / "GSE278687_v2_per_patient.tsv", patients)
    write_tsv(OUT / "GSE278687_v2_summary.tsv", summarize_sections(patients, "patient", rng))
    write_tsv(OUT / "GSE278687_v2_patient_level_block_null.tsv", [block])
    write_tsv(OUT / "GSE278687_v2_patient_level_graph_null.tsv", [graph])
    return sections, patients, block, graph


def run_g277(rng):
    manifest = g277_manifest()
    rows, main_states, main_rows = [], [], []
    with tarfile.open(G277_RAW, "r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith(".tar.gz"):
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in manifest:
                continue
            loaded = read_g277(gsm, outer.extractfile(member).read(), manifest[gsm])
            if loaded is None:
                continue
            sample, matrix, genes, coords, array_rows, array_cols, unit = loaded
            result, state = section_analysis("GSE277116", sample, matrix, genes, coords, array_rows, array_cols, unit, rng)
            result["gsm"] = gsm
            result["stratum"] = manifest[gsm]["sample_type_norm"]
            result["analysis_set"] = "main_18" if manifest[gsm]["tumor_program_scoreability_candidate"] == "1" else "sensitivity_low_detection_frozen"
            rows.append(result)
            if result["analysis_set"] == "main_18":
                main_states.append(state)
                main_rows.append(result)
    rows.sort(key=lambda row: row["gsm"])
    if len(rows) != 21:
        raise RuntimeError(f"Expected 21 GSE277116 tumour packages, found {len(rows)}")
    main = [row for row in rows if row["analysis_set"] == "main_18"]
    graph = aggregate_graph_null(main_states, main_rows, rng, cohort="GSE277116", unit_label="sample")
    summary = summarize_sections(main, "sample (main 18)", rng)
    for stratum in ("ffpe", "frozen"):
        summary.extend(summarize_sections([row for row in main if row["stratum"] == stratum], f"sample ({stratum})", rng))
    summary.extend(summarize_sections(rows, "sample (all 21 sensitivity)", rng))
    write_tsv(OUT / "GSE277116_v2_per_package.tsv", rows)
    write_tsv(OUT / "GSE277116_v2_summary.tsv", summary)
    write_tsv(OUT / "GSE277116_v2_main18_graph_null.tsv", [graph])
    return rows, main, graph


def main():
    rng = np.random.default_rng(SEED)
    g278_sections, g278_patients, block, graph = run_g278(rng)
    g277_rows, g277_main, g277_graph = run_g277(rng)
    manifest = {
        "version": "v2", "seed": SEED,
        "primary_endpoint": "composition-adjusted six-neighbour local mregDC-like--Tfh-like correlation",
        "covariates": {"both": ["log_library", "PTPRC", "epithelial_proxy", "stromal_proxy"], "mregDC_extra": ["DC_core"]},
        "local_field": "six nearest neighbours; focal spot excluded", "mask": "DC-core above within-sample mean",
        "GSE278687_aggregation": "patient median across sections", "GSE277116_aggregation": "sample-level only; no patient IDs available",
        "block_null": {"status": "completed", "method": "within 8x8 array-block outcome permutation; one synchronized cohort-level null draw aggregates section effects to patient medians before cohort median", "draws": N_BLOCK_NULL},
        "graph_null": {"status": "completed", "method": "rank-preserving Moran-matched graph-filtered null on a symmetric six-nearest-neighbour graph; each surrogate draw continuously calibrates graph-filter strength to the observed outcome Moran I", "draws_per_section": N_GRAPH_NULL, "max_graph_filter_steps": GRAPH_MAX_STEPS, "bisection_iterations": GRAPH_MATCH_ITERATIONS},
        "competition": {"status": "completed", "reference_programs": list(COMPETITORS), "constraint": "no genes overlap frozen Tfh-like signature", "additional": ["joint local-field model", "mregDC CCL19-deletion sensitivity"]},
    }
    (OUT / "pipeline_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (OUT / "run_summary.json").write_text(json.dumps({"GSE278687_sections": len(g278_sections), "GSE278687_patients": len(g278_patients), "GSE278687_patient_block_null": block, "GSE278687_patient_graph_null": graph, "GSE277116_packages": len(g277_rows), "GSE277116_main18": len(g277_main), "GSE277116_main18_sample_graph_null": g277_graph}, indent=2), encoding="utf-8")
    print(json.dumps(json.loads((OUT / "run_summary.json").read_text(encoding="utf-8")), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
