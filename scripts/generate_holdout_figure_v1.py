"""Publication-style Figure 5 for leakage-aware spatial holdout sensitivity."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data/03_results/spatial_holdout_validation_v1"
OUT = PROJECT / "figures/main_v2"
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, TEAL, GREY = "#0072B2", "#D55E00", "#009E73", "#7A7A7A"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.titleweight": "bold", "axes.labelsize": 8, "figure.dpi": 300,
    "savefig.dpi": 300, "savefig.bbox": "tight", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": .17,
})


def read(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(row, name):
    return float(row[name])


def plot_panel(ax, rows, title, unit, color):
    rows = sorted(rows, key=lambda x: f(x, "holdout_block_median_r"))
    values = [f(row, "holdout_block_median_r") for row in rows]
    y = np.arange(len(rows))
    ax.scatter(values, y, s=28, color=color, edgecolor="white", linewidth=.5, zorder=3)
    ax.axvline(0, color=GREY, lw=.8)
    ax.axvline(np.median(values), color=TEAL, lw=1.15, ls="--")
    ax.set_yticks(y, [row["unit"] for row in rows])
    ax.set_xlabel("Median held-out-block local r")
    ax.set_title(title, loc="left", pad=9)
    positive = sum(value > 0 for value in values)
    ax.text(.02, .02, f"Median = {np.median(values):.3f}\n{positive}/{len(values)} positive", transform=ax.transAxes, va="bottom")


def main():
    g278 = read(DATA / "GSE278687_patient_holdout_summary.tsv")
    g277 = read(DATA / "GSE277116_sample_holdout_summary.tsv")
    loo278 = read(DATA / "GSE278687_leave_one_patient_out.tsv")
    loo277 = read(DATA / "GSE277116_leave_one_sample_out.tsv")
    fig, axes = plt.subplots(
        1, 3, figsize=(10.8, 4.6), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.15, 1.15, .85]},
    )
    plot_panel(axes[0], g278, "GSE278687\npatient-level held-out blocks", "patient", BLUE)
    plot_panel(axes[1], g277, "GSE277116\nsample-level held-out blocks", "sample", ORANGE)
    ax = axes[2]
    vals278 = [f(row, "cohort_median_r_without_unit") for row in loo278]
    vals277 = [f(row, "cohort_median_r_without_unit") for row in loo277]
    ax.scatter(vals278, np.zeros(len(vals278)), color=BLUE, edgecolor="white", lw=.5, s=25)
    ax.scatter(vals277, np.ones(len(vals277)), color=ORANGE, edgecolor="white", lw=.5, s=25)
    ax.axvline(0, color=GREY, lw=.8)
    ax.set_yticks([0, 1], ["GSE278687", "GSE277116"])
    ax.set_xlabel("Cohort median primary r after omission")
    ax.set_title("Leave-one-unit sensitivity", loc="left", pad=9)
    ax.text(.02, .03, "GSE278687: 0.368-0.468\nGSE277116: 0.407-0.420", transform=ax.transAxes, va="bottom", fontsize=7.2)
    fig.suptitle("Leakage-aware spatial holdout sensitivity", x=.02, ha="left", fontsize=13, weight="bold")
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"Figure5_v1_spatial_holdout_validation.{suffix}")
    print(OUT / "Figure5_v1_spatial_holdout_validation.png")


if __name__ == "__main__":
    main()
