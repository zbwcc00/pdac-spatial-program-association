"""Generate the v6 spatial display and formal supplementary Figures S2-S5.

All displays use locked v2/v3/v4 results. Representative GSE278687 sections
are selected deterministically as the patients nearest the empirical 25th and
75th percentiles of the locked patient-level effect distribution, not image
appearance.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from PIL import Image

PROJECT = Path(__file__).resolve().parents[1]
FIGURES = PROJECT / "figures" / "submission_v6"
SUPP = PROJECT / "manuscript" / "supplementary_v3_final_mask_moran"
V2 = PROJECT / "data" / "03_results" / "unified_primary_pipeline_v2"
V3 = PROJECT / "data" / "03_results" / "spatial_nulls_v3_999"
V4 = PROJECT / "data" / "03_results" / "spatial_nulls_v4_mask_moran"
SCORE = PROJECT / "data" / "03_results" / "scoreability_threshold_sensitivity_v4"
G277 = PROJECT / "data" / "03_results" / "GSE277116_spatial_maps_and_robustness"
QC = PROJECT / "data" / "03_results" / "GSE277116_full_package_qc"
SPATIAL = PROJECT / "data" / "01_unpacked" / "spatial" / "GSE278687_spatial"

BLUE, ORANGE, TEAL, PURPLE, GREY = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6B7280"
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42})


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def read_tsv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_primary_module():
    spec = importlib.util.spec_from_file_location("unified_v2_for_figures", PROJECT / "scripts" / "run_unified_primary_pipeline_v2.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.PROJECT = PROJECT
    module.FROZEN = json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8"))
    module.G278 = PROJECT / "data/01_unpacked/spatial/GSE278687"
    module.G278_COORDS = SPATIAL
    return module


def spatial_panel(ax, image, coords, values, title, cmap, norm=None, mask=None):
    ax.imshow(image, origin="upper")
    if values is not None:
        keep = np.ones(len(coords), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
        scatter = ax.scatter(coords[keep, 0], coords[keep, 1], c=np.asarray(values)[keep], s=7, cmap=cmap, norm=norm, linewidths=0, alpha=0.88)
    else:
        scatter = None
    ax.set_title(title, fontsize=8.6, weight="bold", pad=4)
    ax.set_axis_off()
    return scatter


def generate_main_figure2():
    module = load_primary_module()
    selections = [
        ("GSM8552956_PA40", "PA40", "nearest empirical 25th-percentile effect", 0.3098014452941048),
        ("GSM8552957_PA47", "PA47", "nearest empirical 75th-percentile effect", 0.49527827218234743),
    ]
    records = []
    prepared = []
    for sample, patient, rule, effect in selections:
        loaded = module.read_g278(PROJECT / "data/01_unpacked/spatial/GSE278687" / f"{sample}_filtered_feature_bc_matrix.h5")
        _, matrix, genes, coords, _, _, _ = loaded
        scores = module.add_mreg_no_ccl19(module.make_scores(matrix, genes), matrix, genes)
        mreg, tfh, mask, neighbours, _ = module.prepared_fields(scores, coords)
        local_mreg = module.local_field(module.z(mreg), neighbours)
        local_tfh = module.local_field(module.z(tfh), neighbours)
        joint = np.minimum(module.z(local_mreg), module.z(local_tfh))
        spatial_dir = SPATIAL / sample / "spatial"
        scale = json.loads((spatial_dir / "scalefactors_json.json").read_text(encoding="utf-8"))["tissue_hires_scalef"]
        image = np.asarray(Image.open(spatial_dir / "tissue_hires_image.png").convert("RGB"))
        image_coords = coords * scale
        prepared.append((sample, patient, rule, effect, image, image_coords, scores["DC_core"], local_mreg, local_tfh, joint, mask))
        records.append({"sample": sample, "patient": patient, "selection_rule": rule, "patient_level_primary_local_r": effect, "dc_core_mask_spots": int(mask.sum()), "all_tissue_spots": int(len(mask))})

    all_local = np.concatenate([np.concatenate((item[7], item[8])) for item in prepared])
    lim = float(np.quantile(np.abs(all_local), 0.98))
    joint_lim = float(np.quantile(np.abs(np.concatenate([item[9] for item in prepared])), 0.98))
    fig, axes = plt.subplots(2, 5, figsize=(10.3, 4.35))
    for row, item in enumerate(prepared):
        sample, patient, rule, effect, image, coords, dc, local_mreg, local_tfh, joint, mask = item
        spatial_panel(axes[row, 0], image, coords, None, f"{patient}: tissue image\npatient r={effect:.3f}", "gray")
        dc_norm = plt.Normalize(vmin=np.quantile(dc, 0.02), vmax=np.quantile(dc, 0.98))
        s1 = spatial_panel(axes[row, 1], image, coords, dc, "DC-core score", "YlGnBu", dc_norm)
        s2 = spatial_panel(axes[row, 2], image, coords, local_mreg, "Local mregDC-like\nresidual field", "RdBu_r", TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), mask)
        s3 = spatial_panel(axes[row, 3], image, coords, local_tfh, "Local Tfh-like\nresidual field", "RdBu_r", TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), mask)
        s4 = spatial_panel(axes[row, 4], image, coords, joint, "Local joint residual\nmin(z)", "viridis", TwoSlopeNorm(vcenter=0, vmin=-joint_lim, vmax=joint_lim), mask)
        for axis, scatter in zip((axes[row, 1], axes[row, 2], axes[row, 3], axes[row, 4]), (s1, s2, s3, s4)):
            cbar = fig.colorbar(scatter, ax=axis, fraction=.045, pad=.015)
            cbar.ax.tick_params(labelsize=6)
        axes[row, 0].text(.02, .02, rule, transform=axes[row, 0].transAxes, fontsize=6.2, color="black", bbox={"facecolor":"white", "edgecolor":"none", "alpha":.84, "pad":1.2})
    fig.suptitle("GSE278687 spatial displays selected by patient-effect quantile rule", x=.5, y=1.01, fontsize=11, weight="bold")
    save(fig, FIGURES / "Figure2_v6_GSE278687_representative_spatial_maps")
    with (FIGURES / "Figure2_v6_selection_audit.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(records)


def collect_audits(pattern, column):
    values = []
    for audit_path in sorted(candidate for candidate in pattern.parent.glob(pattern.name) if candidate.is_file()):
        for row in read_tsv(audit_path):
            values.append(float(row[column]))
    return np.asarray(values)


def generate_s3_moran_calibration():
    data = []
    settings = [
        ("GSE278687", "Global target", V3 / "GSE278687_v3_patient_graph_null_audit_batch*.tsv", "global_moran_section_median_abs_delta", "matched global Moran I", TEAL),
        ("GSE278687", "Global target", V3 / "GSE278687_v3_patient_graph_null_audit_batch*.tsv", "mask_moran_section_median_abs_delta", "audited mask Moran I", "#9CA3AF"),
        ("GSE278687", "Mask target", V4 / "GSE278687_v4_patient_mask_moran_audit_batch*.tsv", "mask_moran_section_median_abs_delta", "matched mask Moran I", PURPLE),
        ("GSE278687", "Mask target", V4 / "GSE278687_v4_patient_mask_moran_audit_batch*.tsv", "global_moran_section_median_abs_delta", "audited global Moran I", "#9CA3AF"),
        ("GSE277116", "Global target", V3 / "GSE277116_v3_main18_graph_null_audit_batch*.tsv", "global_moran_section_median_abs_delta", "matched global Moran I", TEAL),
        ("GSE277116", "Global target", V3 / "GSE277116_v3_main18_graph_null_audit_batch*.tsv", "mask_moran_section_median_abs_delta", "audited mask Moran I", "#9CA3AF"),
        ("GSE277116", "Mask target", V4 / "GSE277116_v4_sample_mask_moran_audit_batch*.tsv", "mask_moran_section_median_abs_delta", "matched mask Moran I", PURPLE),
        ("GSE277116", "Mask target", V4 / "GSE277116_v4_sample_mask_moran_audit_batch*.tsv", "global_moran_section_median_abs_delta", "audited global Moran I", "#9CA3AF"),
    ]
    for cohort, target, pattern, column, label, color in settings:
        values = collect_audits(pattern, column)
        data.append((cohort, target, label, color, values))
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.15), sharey=True)
    for ax, cohort in zip(axes, ("GSE278687", "GSE277116")):
        subset = [item for item in data if item[0] == cohort]
        positions = np.arange(4)
        box = ax.boxplot([item[4] for item in subset], positions=positions, widths=.58, patch_artist=True, showfliers=False, medianprops={"color":"black", "linewidth":1.1})
        for patch, item in zip(box["boxes"], subset):
            patch.set_facecolor(item[3]); patch.set_alpha(.78)
        for pos, item in zip(positions, subset):
            ax.text(pos, max(np.quantile(item[4], .95), 1e-5) * 1.65, f"median {np.median(item[4]):.2e}", ha="center", fontsize=6.1, rotation=25)
        ax.set_xticks(positions, ["Global\nmatched", "Mask\naudit", "Mask\nmatched", "Global\naudit"])
        ax.set_title(cohort, loc="left", weight="bold")
        ax.set_yscale("log"); ax.grid(axis="y", alpha=.2); ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Per-draw section-median absolute Moran I mismatch")
    fig.suptitle("Complementary Moran calibration and cross-scale audits (1,000 draws each)", y=1.02, fontsize=10.5, weight="bold")
    save(fig, FIGURES / "Supplementary_Figure_S3_Moran_calibration_and_audits")


def generate_s5_scoreability():
    rows = read_tsv(SCORE / "scoreability_threshold_summary_v4.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 2.9), sharey=True)
    for ax, cohort, color in zip(axes, ("GSE278687", "GSE277116"), (BLUE, ORANGE)):
        subset = [row for row in rows if row["cohort"] == cohort]
        x = np.asarray([float(row["threshold"]) for row in subset])
        y = np.asarray([float(row["median_primary_local_r"]) for row in subset])
        lo = np.asarray([float(row["bootstrap_ci_low"]) for row in subset])
        hi = np.asarray([float(row["bootstrap_ci_high"]) for row in subset])
        n = [int(row["n_eligible_units"]) for row in subset]
        ax.errorbar(x, y, yerr=[y-lo, hi-y], marker="o", color=color, lw=1.5, capsize=2)
        for xi, yi, ni in zip(x, y, n): ax.text(xi, yi+.034, f"n={ni}", ha="center", color=color, fontsize=7)
        ax.set_title(cohort, loc="left", weight="bold"); ax.set_xlabel("Joint primary-program nonzero-score threshold")
        ax.set_xticks(x, [f"{value:.1f}" for value in x]); ax.grid(axis="y", alpha=.2); ax.spines[["top","right"]].set_visible(False)
        ax.axhline(.4183273 if cohort == "GSE278687" else .4136335, color=GREY, lw=.8, ls="--", label="Locked main estimate")
        ax.legend(frameon=False, fontsize=6.8, loc="lower left")
    axes[0].set_ylabel("Median local correlation (95% bootstrap CI)")
    fig.suptitle("Predefined scoreability-threshold sensitivity", y=1.02, fontsize=10.5, weight="bold")
    save(fig, FIGURES / "Supplementary_Figure_S5_scoreability_threshold_sensitivity")


def copy_supplementary_existing():
    shutil.copy2(G277 / "GSE277116_representative_spatial_maps.png", FIGURES / "Supplementary_Figure_S2_GSE277116_representative_spatial_maps.png")
    shutil.copy2(G277 / "GSE277116_representative_spatial_maps.pdf", FIGURES / "Supplementary_Figure_S2_GSE277116_representative_spatial_maps.pdf")
    shutil.copy2(QC / "GSE277116_full_package_qc_overview.png", FIGURES / "Supplementary_Figure_S4_GSE277116_full_package_QC.png")
    image = plt.imread(QC / "GSE277116_full_package_qc_overview.png")
    fig, ax = plt.subplots(figsize=(7.2, 7.0)); ax.imshow(image); ax.set_axis_off(); fig.tight_layout(pad=0)
    fig.savefig(FIGURES / "Supplementary_Figure_S4_GSE277116_full_package_QC.pdf", bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def publish_to_supplement():
    SUPP.mkdir(parents=True, exist_ok=True)
    for stem in ("Supplementary_Figure_S2_GSE277116_representative_spatial_maps", "Supplementary_Figure_S3_Moran_calibration_and_audits", "Supplementary_Figure_S4_GSE277116_full_package_QC", "Supplementary_Figure_S5_scoreability_threshold_sensitivity"):
        for extension in (".png", ".pdf"):
            source = FIGURES / f"{stem}{extension}"
            if source.exists(): shutil.copy2(source, SUPP / source.name)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    generate_main_figure2()
    copy_supplementary_existing()
    generate_s3_moran_calibration()
    generate_s5_scoreability()
    publish_to_supplement()
    print(FIGURES)


if __name__ == "__main__":
    main()
