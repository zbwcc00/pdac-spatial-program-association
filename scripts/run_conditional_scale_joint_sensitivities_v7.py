"""Post-hoc sensitivity analyses for the conditional local-program estimand.

This script intentionally does not replace the locked v2 primary endpoint.
It quantifies how the observed local correlation varies with (i) DC-core
conditioning, (ii) the DC-core mask threshold, and (iii) neighbourhood size.
It also reports coordinate-scale neighbour distances and descriptive joint
model VIF diagnostics. GSE278687 is aggregated to patients; GSE277116 remains
sample-level because patient identifiers are unavailable.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "data/03_results/conditional_scale_joint_sensitivities_v7"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260814


def load_primary():
    path = PROJECT / "scripts/run_unified_primary_pipeline_v2.py"
    spec = importlib.util.spec_from_file_location("unified_v2_sensitivity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.PROJECT = PROJECT
    module.FROZEN = json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8"))
    module.G278 = PROJECT / "data/01_unpacked/spatial/GSE278687"
    module.G278_COORDS = PROJECT / "data/01_unpacked/spatial/GSE278687_spatial"
    module.G277_RAW = PROJECT / "data/00_raw/spatial/GSE277116_RAW.tar"
    module.G277_QC = PROJECT / "data/03_results/GSE277116_full_package_qc/GSE277116_tumor_replication_manifest.tsv"
    return module


def write_tsv(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def median_ci(values, rng):
    values = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if not len(values): return (np.nan, np.nan, np.nan, 0)
    boot = np.asarray([np.median(rng.choice(values, len(values), replace=True)) for _ in range(2000)])
    return float(np.median(values)), float(np.quantile(boot, .025)), float(np.quantile(boot, .975)), int((values > 0).sum())


def vif(values):
    """Columnwise variance inflation factors, descriptive only."""
    output = []
    for j in range(values.shape[1]):
        y = values[:, j]
        x = np.delete(values, j, axis=1)
        design = np.column_stack([np.ones(len(y)), x])
        fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
        total = np.sum((y - y.mean()) ** 2)
        r2 = 1 - np.sum((y - fitted) ** 2) / total if total > 0 else np.nan
        output.append(float(1 / max(1 - r2, 1e-12)) if np.isfinite(r2) else np.nan)
    return output


def analyse_unit(module, cohort, sample, matrix, genes, coords, rows, cols, unit):
    scores = module.add_mreg_no_ccl19(module.make_scores(matrix, genes), matrix, genes)
    base_covariates = (scores["library"], scores["PTPRC"], scores["epithelial"], scores["stromal"])
    mreg_unadjusted = module.residual(scores["mregDC_strict"], *base_covariates)
    mreg_adjusted = module.residual(scores["mregDC_strict"], *base_covariates, scores["DC_core"])
    tfh = module.residual(scores["Tfh_like"], *base_covariates)
    dc = scores["DC_core"]
    threshold_defs = {
        "whole_tissue": np.ones(len(dc), dtype=bool),
        "DC_core_above_median": dc > np.median(dc),
        "DC_core_above_mean": module.z(dc) > 0,
        "DC_core_above_q75": dc > np.quantile(dc, .75),
    }
    fields = []
    distances = []
    vif_rows = []
    for k in (4, 6, 12):
        neighbours = module.neighbour_index(coords, k=k)
        # kNN distances are raw coordinate units, not assumed to be micrometres.
        delta = coords[:, None, :] - coords[neighbours]
        d = np.sqrt(np.sum(delta ** 2, axis=2)).ravel()
        distances.append({"cohort": cohort, "sample": sample, "unit": unit, "k": k, "coordinate_unit": "image_pixels" if cohort == "GSE278687" else "archive_array_coordinates", "n_edges": len(d), "distance_q05": float(np.quantile(d,.05)), "distance_q25": float(np.quantile(d,.25)), "distance_median": float(np.median(d)), "distance_q75": float(np.quantile(d,.75)), "distance_q95": float(np.quantile(d,.95)), "distance_mean": float(np.mean(d))})
        for adjustment, mreg in (("without_DC_core_adjustment", mreg_unadjusted), ("with_DC_core_adjustment", mreg_adjusted)):
            left = module.local_field(module.z(mreg), neighbours)
            right = module.local_field(module.z(tfh), neighbours)
            for mask_name, mask in threshold_defs.items():
                fields.append({"cohort": cohort, "sample": sample, "unit": unit, "k": k, "mregDC_residualization": adjustment, "mask_definition": mask_name, "mask_spots": int(mask.sum()), "mask_fraction": float(mask.mean()), "local_r": module.correlation(left[mask], right[mask])})
        # Joint fields and VIF are only descriptive. Match the locked primary
        # k=6/DC-core mask definition so they can be interpreted beside Figure 3.
        if k == 6:
            mask = threshold_defs["DC_core_above_mean"]
            local_mreg = module.local_field(module.z(mreg_adjusted), neighbours)
            predictors = [module.local_field(module.z(tfh), neighbours)]
            labels = ["Tfh_like"]
            for name in module.COMPETITORS:
                residual = module.residual(scores[name], *base_covariates)
                predictors.append(module.local_field(module.z(residual), neighbours))
                labels.append(name)
            valid = mask & np.isfinite(local_mreg) & np.all(np.isfinite(np.column_stack(predictors)), axis=1)
            x = np.column_stack([module.z(field[valid]) for field in predictors])
            beta = np.linalg.lstsq(np.column_stack([np.ones(valid.sum()), x]), module.z(local_mreg[valid]), rcond=None)[0][1:]
            for label, this_vif, this_beta in zip(labels, vif(x), beta):
                vif_rows.append({"cohort": cohort, "sample": sample, "unit": unit, "predictor": label, "spots": int(valid.sum()), "VIF": this_vif, "standardized_beta_descriptive": float(this_beta)})
    return fields, distances, vif_rows


def main():
    module = load_primary(); rng = np.random.default_rng(SEED)
    all_fields, all_distances, all_vif = [], [], []
    section_to_patient = {}
    for path in sorted(module.G278.glob("*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, coords, rows, cols, patient = module.read_g278(path)
        fields, distances, vifs = analyse_unit(module, "GSE278687", sample, matrix, genes, coords, rows, cols, patient)
        all_fields.extend(fields); all_distances.extend(distances); all_vif.extend(vifs); section_to_patient[sample] = patient
    manifest = module.g277_manifest()
    # The full-package QC manifest predates the final main-18 designation.
    # Read the locked v2 per-package table rather than assuming that every
    # technically viable tumour package entered the external main analysis.
    with (PROJECT / "data/03_results/unified_primary_pipeline_v2/GSE277116_v2_per_package.tsv").open(encoding="utf-8") as handle:
        main18 = {row["gsm"] for row in csv.DictReader(handle, delimiter="\t") if row.get("analysis_set") == "main_18"}
    import tarfile
    with tarfile.open(module.G277_RAW, "r") as outer:
        members = {Path(m.name).name: m for m in outer.getmembers() if m.isfile()}
        for gsm, meta in manifest.items():
            if gsm not in main18: continue
            member = next((m for name, m in members.items() if name.startswith(f"{gsm}_") and name.endswith(".tar.gz")), None)
            if member is None: continue
            loaded = module.read_g277(gsm, outer.extractfile(member).read(), meta)
            if loaded is None: continue
            sample, matrix, genes, coords, rows, cols, unit = loaded
            fields, distances, vifs = analyse_unit(module, "GSE277116", sample, matrix, genes, coords, rows, cols, unit)
            all_fields.extend(fields); all_distances.extend(distances); all_vif.extend(vifs)
    write_tsv(OUT / "per_section_conditional_scale_effects.tsv", all_fields)
    write_tsv(OUT / "per_section_neighbour_distance_distributions.tsv", all_distances)
    write_tsv(OUT / "per_section_joint_model_vif.tsv", all_vif)
    # Aggregate GSE278687 at patient level for every sensitivity combination.
    combined = []
    for cohort in ("GSE278687", "GSE277116"):
        subset = [r for r in all_fields if r["cohort"] == cohort]
        groups = sorted({(r["k"], r["mregDC_residualization"], r["mask_definition"]) for r in subset})
        for key in groups:
            rows = [r for r in subset if (r["k"],r["mregDC_residualization"],r["mask_definition"]) == key]
            if cohort == "GSE278687":
                by_unit = {u: [] for u in sorted({r["unit"] for r in rows})}
                for row in rows: by_unit[row["unit"]].append(row["local_r"])
                values = [float(np.nanmedian(v)) for v in by_unit.values()]
            else:
                values = [r["local_r"] for r in rows]
            median, lo, hi, positive = median_ci(values, rng)
            combined.append({"cohort": cohort, "analysis_unit": "patient" if cohort == "GSE278687" else "sample", "k":key[0], "mregDC_residualization":key[1], "mask_definition":key[2], "n_units":len(values), "median_local_r":median, "bootstrap_ci_low":lo, "bootstrap_ci_high":hi, "positive_units":positive})
    write_tsv(OUT / "cohort_conditional_scale_summary.tsv", combined)
    # Summarise VIF/betas at the same hierarchy, descriptive only.
    vif_summary=[]
    for cohort in ("GSE278687","GSE277116"):
        for predictor in sorted({r["predictor"] for r in all_vif if r["cohort"]==cohort}):
            rows=[r for r in all_vif if r["cohort"]==cohort and r["predictor"]==predictor]
            if cohort == "GSE278687":
                unit_vif={u:[] for u in {r["unit"] for r in rows}}; unit_beta={u:[] for u in {r["unit"] for r in rows}}
                for r in rows: unit_vif[r["unit"]].append(r["VIF"]); unit_beta[r["unit"]].append(r["standardized_beta_descriptive"])
                vifs=[np.nanmedian(x) for x in unit_vif.values()]; betas=[np.nanmedian(x) for x in unit_beta.values()]
            else:
                vifs=[r["VIF"] for r in rows]; betas=[r["standardized_beta_descriptive"] for r in rows]
            vif_summary.append({"cohort":cohort,"analysis_unit":"patient" if cohort=="GSE278687" else "sample","predictor":predictor,"n_units":len(vifs),"median_VIF":float(np.nanmedian(vifs)),"max_VIF":float(np.nanmax(vifs)),"median_standardized_beta_descriptive":float(np.nanmedian(betas))})
    write_tsv(OUT / "joint_model_vif_summary.tsv", vif_summary)
    (OUT / "README.md").write_text("# Post-hoc conditional, scale, and joint-model sensitivity analyses\n\nThese analyses were added after review of the locked v2 endpoint. They do not alter the primary estimand or create confirmatory claims. Coordinates are measured in per-archive raw units; cross-cohort physical-distance equivalence is not inferred. Joint-model betas and VIF values are descriptive because residual spatial correlation is not modeled.\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__": main()
