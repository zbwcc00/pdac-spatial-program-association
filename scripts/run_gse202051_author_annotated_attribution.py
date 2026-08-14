"""Author-annotated PDAC single-cell program-attribution and specificity QC.

The analysis uses the frozen spatial gene lists without reselecting genes or
reannotating cells. It is deliberately not a spatial-validation, interaction,
TLS, causal, prognostic, or treatment-response analysis.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "00_raw" / "single_cell" / "GSE202051"
H5AD = DATA / "GSE202051_totaldata-final-toshare.h5ad"
GZIP = DATA / "GSE202051_totaldata-final-toshare.h5ad.gz"
FROZEN = json.loads(
    (PROJECT / "data" / "03_results" / "GSE154778_program_freeze" / "frozen_programs.json").read_text(encoding="utf-8")
)
OUT = PROJECT / "data" / "03_results" / "GSE202051_author_annotated_attribution_v1"
FIG = PROJECT / "figures" / "GSE202051_author_annotated_attribution_v1"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260814
N_BINS = 25
CONTROLS_PER_TARGET = 50
N_BOOTSTRAP = 10_000
LABEL = "Level 3 Annotation"
PID = "pid"
SAMPLE = "sampleid"
TREATMENT = "new_treatment"

PROGRAMS = {
    "mregDC_like": tuple(FROZEN["mregDC_strict"]["genes"]),
    "Tfh_like": tuple(FROZEN["Tfh_like"]["genes"]),
    "DC_core": tuple(FROZEN["DC_core"]["genes"]),
    "broad_T": ("TRAC", "TRBC1", "TRBC2", "CD247", "LCK", "CD2", "CD7"),
    "non_Tfh_CD4": ("LTB", "LEF1", "TCF7", "SELL", "MAL"),
    "Treg": ("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF4"),
    "exhausted_CD8": ("CD8A", "CD8B", "LAG3", "HAVCR2", "ENTPD1", "LAYN"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def same_sparse_pattern(path: Path) -> bool:
    """Verify that X and the raw-UMI counts layer have the same zero pattern."""
    with h5py.File(path, "r") as handle:
        x, counts = handle["X"], handle["layers"]["counts"]
        if x["indices"].shape != counts["indices"].shape or x["indptr"].shape != counts["indptr"].shape:
            return False
        for key in ("indices", "indptr"):
            left, right = x[key], counts[key]
            for start in range(0, left.shape[0], 5_000_000):
                stop = min(start + 5_000_000, left.shape[0])
                if not np.array_equal(left[start:stop], right[start:stop]):
                    return False
    return True


def expression_bins(matrix: sparse.csr_matrix) -> np.ndarray:
    means = np.asarray(matrix.mean(axis=0)).ravel()
    order = np.argsort(means, kind="mergesort")
    ranks = np.empty_like(order)
    ranks[order] = np.arange(len(order))
    return np.minimum((ranks * N_BINS) // len(order), N_BINS - 1)


def select_controls(
    program: tuple[str, ...], gene_index: dict[str, int], bins: np.ndarray, excluded: set[int], rng: np.random.Generator
) -> tuple[list[int], list[str]]:
    available = [gene for gene in program if gene in gene_index]
    chosen: set[int] = set()
    for gene in available:
        candidates = np.flatnonzero(bins == bins[gene_index[gene]])
        candidates = np.array([index for index in candidates if index not in excluded], dtype=int)
        if len(candidates) == 0:
            continue
        take = min(CONTROLS_PER_TARGET, len(candidates))
        chosen.update(rng.choice(candidates, size=take, replace=False).tolist())
    return sorted(chosen), available


def vector_mean(matrix: sparse.csr_matrix, columns: list[int]) -> np.ndarray:
    if not columns:
        return np.zeros(matrix.shape[0], dtype=float)
    return np.asarray(matrix[:, columns].mean(axis=1)).ravel()


def score_programs(matrix: sparse.csr_matrix, var_names: pd.Index) -> tuple[dict[str, np.ndarray], dict[str, dict]]:
    upper = pd.Index([str(name).upper() for name in var_names])
    if not upper.is_unique:
        raise RuntimeError("Gene symbols are not unique; cannot score frozen lists deterministically.")
    gene_index = {gene: index for index, gene in enumerate(upper)}
    all_frozen = {gene_index[gene] for genes in PROGRAMS.values() for gene in genes if gene in gene_index}
    variants = dict(PROGRAMS)
    variants["mregDC_minus_CCL19"] = tuple(gene for gene in PROGRAMS["mregDC_like"] if gene != "CCL19")
    for gene in PROGRAMS["mregDC_like"]:
        variants[f"mregDC_LOO_{gene}"] = tuple(item for item in PROGRAMS["mregDC_like"] if item != gene)
    for gene in PROGRAMS["Tfh_like"]:
        variants[f"Tfh_LOO_{gene}"] = tuple(item for item in PROGRAMS["Tfh_like"] if item != gene)
    bins = expression_bins(matrix)
    rng = np.random.default_rng(SEED)
    scores, manifest = {}, {}
    for name, genes in variants.items():
        controls, available = select_controls(genes, gene_index, bins, all_frozen, rng)
        module = vector_mean(matrix, [gene_index[gene] for gene in available])
        control = vector_mean(matrix, controls)
        scores[f"{name}_uncontrolled_mean"] = module
        scores[name] = module - control
        manifest[f"{name}_uncontrolled_mean"] = {
            "frozen_genes": list(genes),
            "available_genes": available,
            "missing_genes": [gene for gene in genes if gene not in gene_index],
            "n_controls": 0,
            "control_gene_symbols": [],
            "control_selection": "none; sensitivity score is the mean log-normalized X expression of the frozen genes",
        }
        manifest[name] = {
            "frozen_genes": list(genes),
            "available_genes": available,
            "missing_genes": [gene for gene in genes if gene not in gene_index],
            "n_controls": len(controls),
            "control_gene_symbols": [str(var_names[index]) for index in controls],
            "control_selection": f"expression-bin matched; {N_BINS} bins; up to {CONTROLS_PER_TARGET} controls per target; seed={SEED}",
        }
    return scores, manifest


def patient_label_table(meta: pd.DataFrame, score: np.ndarray, score_name: str) -> pd.DataFrame:
    frame = meta[[PID, SAMPLE, TREATMENT, LABEL]].copy()
    frame[score_name] = score
    return (
        frame.groupby([PID, LABEL], observed=True)
        .agg(
            cells=(score_name, "size"),
            samples=(SAMPLE, "nunique"),
            treatment_states=(TREATMENT, "nunique"),
            median_score=(score_name, "median"),
            mean_score=(score_name, "mean"),
            score_iqr_low=(score_name, lambda value: value.quantile(0.25)),
            score_iqr_high=(score_name, lambda value: value.quantile(0.75)),
        )
        .reset_index()
        .assign(program=score_name)
    )


def bootstrap_median_ci(values: pd.Series, seed_offset: int) -> tuple[float, float]:
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(SEED + seed_offset)
    draws = rng.choice(values.to_numpy(dtype=float), size=(N_BOOTSTRAP, len(values)), replace=True)
    medians = np.median(draws, axis=1)
    return float(np.quantile(medians, 0.025)), float(np.quantile(medians, 0.975))


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna()
    if valid.empty:
        return adjusted
    ordered = valid.sort_values()
    ranks = np.arange(1, len(ordered) + 1, dtype=float)
    values = ordered.to_numpy(dtype=float) * len(ordered) / ranks
    values = np.minimum.accumulate(values[::-1])[::-1]
    adjusted.loc[ordered.index] = np.minimum(values, 1.0)
    return adjusted


def paired_test(
    table: pd.DataFrame,
    program: str,
    left: str,
    right: str,
    purpose: str,
    *,
    min_cells_per_label: int = 1,
    inferential_role: str = "exploratory",
) -> dict:
    subset = table[table["program"] == program]
    score_pivot = subset.pivot(index=PID, columns=LABEL, values="median_score")
    count_pivot = subset.pivot(index=PID, columns=LABEL, values="cells")
    base = {
        "purpose": purpose,
        "program": program,
        "left_label": left,
        "right_label": right,
        "min_cells_per_label": min_cells_per_label,
        "inferential_role": inferential_role,
    }
    if left not in score_pivot or right not in score_pivot:
        return {**base, "paired_patients": 0, "status": "missing author label"}
    retained = (count_pivot[left] >= min_cells_per_label) & (count_pivot[right] >= min_cells_per_label)
    paired = score_pivot.loc[retained, [left, right]].dropna()
    delta = paired[left] - paired[right]
    try:
        result = wilcoxon(delta, alternative="two-sided", method="auto")
        p_value = float(result.pvalue)
    except ValueError:
        p_value = np.nan
    ci_low, ci_high = bootstrap_median_ci(delta, seed_offset=10_000 + min_cells_per_label)
    return {
        **base,
        "paired_patients": int(len(paired)),
        "median_paired_difference": float(delta.median()) if len(delta) else np.nan,
        "median_paired_difference_ci_low": ci_low,
        "median_paired_difference_ci_high": ci_high,
        "median_left": float(paired[left].median()) if len(paired) else np.nan,
        "median_right": float(paired[right].median()) if len(paired) else np.nan,
        "wilcoxon_two_sided_p": p_value,
        "status": "patient-level program-attribution comparison; not a spatial, causal, or cell-interaction test",
    }


def plot_patient_scores(table: pd.DataFrame) -> None:
    panels = [
        ("mregDC_like", ["Dendritic (activated)", "Dendritic (conventional type 2)", "Dendritic (conventional type 1)", "Dendritic (plasmacytoid)", "Macrophage"], "mregDC-like matched-control score"),
        ("mregDC_minus_CCL19", ["Dendritic (activated)", "Dendritic (conventional type 2)", "Dendritic (conventional type 1)", "Dendritic (plasmacytoid)", "Macrophage"], "mregDC-like score without CCL19"),
        ("Tfh_like", ["CD4+ T", "Treg", "Treg (activated)", "CD8+ T", "CD8+ T (dysfunctional)", "B", "Plasma"], "Tfh-like matched-control score"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    for ax, (program, labels, ylabel) in zip(axes, panels):
        data = [table[(table["program"] == program) & (table[LABEL] == label)]["median_score"].dropna().to_numpy() for label in labels]
        plotted = [(label, values) for label, values in zip(labels, data) if len(values)]
        ax.boxplot([values for _, values in plotted], showfliers=False, widths=0.55, patch_artist=True, boxprops={"facecolor": "#a8dadc", "edgecolor": "#264653"}, medianprops={"color": "#e63946", "linewidth": 1.5})
        jitter = np.random.default_rng(SEED).normal(0, 0.045, sum(len(values) for _, values in plotted))
        offset = 0
        for position, (_, values) in enumerate(plotted, start=1):
            ax.scatter(np.full(len(values), position) + jitter[offset:offset + len(values)], values, s=11, color="#264653", alpha=0.55, linewidths=0)
            offset += len(values)
        ax.axhline(0, color="#6c757d", linewidth=0.8)
        ax.set_xticks(range(1, len(plotted) + 1), [label.replace("Dendritic ", "DC ").replace(" (conventional type ", " cDC").replace(")", "") for label, _ in plotted], rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Patient medians; n={sum(len(values) for _, values in plotted)} label-patient units")
    fig.savefig(FIG / "Figure_S9_author_annotated_program_attribution.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "Figure_S9_author_annotated_program_attribution.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not H5AD.exists() or not GZIP.exists():
        raise FileNotFoundError("Both complete H5AD and its compressed source are required.")
    adata = ad.read_h5ad(H5AD, backed="r")
    required = {PID, SAMPLE, TREATMENT, LABEL}
    absent = required.difference(adata.obs.columns)
    if absent:
        raise RuntimeError(f"Required author metadata missing: {sorted(absent)}")
    if not adata.obs_names.is_unique:
        raise RuntimeError("Barcodes/observation identifiers are not unique.")
    meta = adata.obs[[PID, SAMPLE, TREATMENT, LABEL]].copy()
    # X is the authors' non-negative normalized expression matrix. It is loaded
    # once as sparse CSR; raw count detection is audited separately by pattern.
    x = adata.X.to_memory().tocsr()
    if x.shape != adata.shape or x.min() < 0:
        raise RuntimeError("Unexpected expression layer; expected non-negative normalized X.")
    scores, program_manifest = score_programs(x, adata.var_names)
    pattern_matches = same_sparse_pattern(H5AD)
    patient_tables = [patient_label_table(meta, score, name) for name, score in scores.items()]
    patient_table = pd.concat(patient_tables, ignore_index=True)
    label_crosswalk = (
        meta.groupby(LABEL, observed=True)
        .agg(cells=(PID, "size"), patients=(PID, "nunique"), samples=(SAMPLE, "nunique"), treatment_states=(TREATMENT, "nunique"))
        .reset_index()
        .sort_values("cells", ascending=False)
    )
    primary_names = list(PROGRAMS) + ["mregDC_minus_CCL19", "mregDC_like_uncontrolled_mean", "Tfh_like_uncontrolled_mean"]
    label_summaries = []
    label_positions = meta.groupby(LABEL, observed=True).indices
    for name in primary_names:
        for label, index in label_positions.items():
            values = scores[name][np.asarray(index, dtype=int)]
            label_summaries.append({
                "program": name,
                LABEL: label,
                "cells": int(len(values)),
                "patients": int(meta.iloc[np.asarray(index, dtype=int)][PID].nunique()),
                "median_cell_score": float(np.median(values)),
                "cell_score_iqr_low": float(np.quantile(values, 0.25)),
                "cell_score_iqr_high": float(np.quantile(values, 0.75)),
                "positive_fraction": float(np.mean(values > 0)),
            })
    scoreability = []
    upper_var = pd.Index([str(gene).upper() for gene in adata.var_names])
    for name in primary_names:
        genes = program_manifest[name]["available_genes"]
        for label, index in label_positions.items():
            idx = np.asarray(index, dtype=int)
            scoreability.append({
                "program": name,
                LABEL: label,
                "cells": int(len(idx)),
                "gene_coverage": len(genes) / len(program_manifest[name]["frozen_genes"]),
                "score_zero_fraction": float(np.mean(vector_mean(x[idx, :], [upper_var.get_loc(gene) for gene in genes]) == 0)),
                "detection_layer": "raw UMI counts (same sparse zero pattern as X verified globally)" if pattern_matches else "X zero pattern only; counts mismatch",
            })
            for gene in genes:
                expression = x[idx, upper_var.get_loc(gene)]
                scoreability.append({
                    "program": name,
                    LABEL: label,
                    "cells": int(len(idx)),
                    "gene": gene,
                    "gene_detection_fraction": float(expression.count_nonzero() / len(idx)),
                    "detection_layer": "raw UMI counts (same sparse zero pattern as X verified globally)" if pattern_matches else "X zero pattern only; counts mismatch",
                })
    comparisons = [
        paired_test(patient_table, "mregDC_like", "Dendritic (activated)", "Dendritic (conventional type 2)", "mregDC attribution primary", inferential_role="primary"),
        paired_test(patient_table, "mregDC_like", "Dendritic (activated)", "Dendritic (conventional type 1)", "mregDC attribution comparator"),
        paired_test(patient_table, "mregDC_like", "Dendritic (activated)", "Dendritic (plasmacytoid)", "mregDC attribution comparator"),
        paired_test(patient_table, "mregDC_like", "Dendritic (activated)", "Macrophage", "mregDC lineage competitor"),
        paired_test(patient_table, "mregDC_minus_CCL19", "Dendritic (activated)", "Dendritic (conventional type 2)", "CCL19 deletion sensitivity"),
        paired_test(patient_table, "mregDC_like_uncontrolled_mean", "Dendritic (activated)", "Dendritic (conventional type 2)", "matched-control scoring sensitivity"),
        paired_test(patient_table, "Tfh_like", "CD4+ T", "Treg", "Tfh-like non-specificity boundary"),
        paired_test(patient_table, "Tfh_like", "CD4+ T", "CD8+ T", "Tfh-like non-specificity boundary"),
        paired_test(patient_table, "Tfh_like", "CD4+ T", "B", "Tfh-like lymphoid localization"),
        paired_test(patient_table, "Tfh_like_uncontrolled_mean", "CD4+ T", "Treg", "matched-control scoring sensitivity"),
    ]
    for name in [item for item in scores if item.startswith("mregDC_LOO_")]:
        comparisons.append(paired_test(patient_table, name, "Dendritic (activated)", "Dendritic (conventional type 2)", "mregDC leave-one-gene-out sensitivity"))
    for name in [item for item in scores if item.startswith("Tfh_LOO_")]:
        comparisons.append(paired_test(patient_table, name, "CD4+ T", "Treg", "Tfh-like leave-one-gene-out non-specificity sensitivity"))
    for minimum in (3, 5):
        comparisons.append(
            paired_test(
                patient_table,
                "mregDC_like",
                "Dendritic (activated)",
                "Dendritic (conventional type 2)",
                "mregDC attribution minimum-cell sensitivity",
                min_cells_per_label=minimum,
            )
        )
    comparison_table = pd.DataFrame(comparisons)
    exploratory = comparison_table["inferential_role"].eq("exploratory")
    comparison_table["wilcoxon_two_sided_p_fdr_bh_exploratory"] = np.nan
    comparison_table.loc[exploratory, "wilcoxon_two_sided_p_fdr_bh_exploratory"] = benjamini_hochberg(
        comparison_table.loc[exploratory, "wilcoxon_two_sided_p"]
    )
    comparison_table["multiplicity_note"] = np.where(
        comparison_table["inferential_role"].eq("primary"),
        "sole designated primary single-cell attribution comparison; not prospectively preregistered",
        "exploratory comparison; Benjamini-Hochberg FDR adjusted across all non-primary comparisons in this table",
    )
    label_crosswalk.to_csv(OUT / "author_label_crosswalk.tsv", sep="\t", index=False)
    pd.DataFrame(label_summaries).to_csv(OUT / "program_score_by_author_label.tsv", sep="\t", index=False)
    patient_table.to_csv(OUT / "patient_level_label_program_scores.tsv", sep="\t", index=False)
    pd.DataFrame(scoreability).to_csv(OUT / "program_scoreability_by_author_label.tsv", sep="\t", index=False)
    comparison_table.to_csv(OUT / "patient_level_primary_and_sensitivity_tests.tsv", sep="\t", index=False)
    comparison_table.loc[
        comparison_table["purpose"].eq("mregDC attribution minimum-cell sensitivity")
    ].to_csv(OUT / "mregDC_activated_dc_minimum_cell_sensitivity.tsv", sep="\t", index=False)
    protocol = {
        "dataset": "GSE202051; Hwang et al., Nature Genetics 2022; PMID 35902743",
        "input_gzip_bytes": GZIP.stat().st_size,
        "input_gzip_sha256": sha256(GZIP),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "unique_observation_ids": bool(adata.obs_names.is_unique),
        "patients": int(meta[PID].nunique()),
        "samples": int(meta[SAMPLE].nunique()),
        "author_annotation_field": LABEL,
        "author_tfh_label_present": bool(meta[LABEL].astype(str).str.contains("tfh|follicular helper", case=False, regex=True).any()),
        "expression_scoring": "primary: mean log-normalized X expression of frozen genes minus expression-bin matched control-gene mean; sensitivity: unadjusted mean log-normalized X expression",
        "count_detection_audit": "X and layers/counts share an identical sparse zero pattern" if pattern_matches else "X and layers/counts sparse patterns differ",
        "sole_primary_test": "matched-control mregDC-like score: author-labelled activated DC versus cDC2, two-sided paired Wilcoxon test",
        "primary_test_boundary": "designated for this single-cell module; not prospectively preregistered and not a spatial-association test",
        "exploratory_multiplicity": f"all other comparisons, score variants, leave-one-gene-out tests, and >=3/>=5 cell-per-label sensitivities receive Benjamini-Hochberg FDR adjustment; {N_BOOTSTRAP} bootstrap resamples for median-difference 95% confidence intervals",
        "programs": program_manifest,
        "scope_boundary": "program attribution/specificity QC only; no spatial validation, cell interaction, TLS, causality, prognosis, or treatment-response claim",
        "key_limitation": "No author-defined Tfh state is available; Tfh-like results assess broader T-cell/lymphoid localization and non-specificity only.",
    }
    (OUT / "analysis_protocol.json").write_text(json.dumps(protocol, indent=2, ensure_ascii=False), encoding="utf-8")
    plot_patient_scores(patient_table)
    print(json.dumps({"output": str(OUT), "figure": str(FIG), "patients": protocol["patients"], "author_tfh_label_present": protocol["author_tfh_label_present"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
