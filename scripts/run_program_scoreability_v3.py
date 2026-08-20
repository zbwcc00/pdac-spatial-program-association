"""Recompute program scoreability from the locked spatial scoring workflow.

For each primary or competitor program, scoreability is the fraction of
in-tissue spots with a nonzero locked log1p summed-count score. GSE278687
reports all 21 sections. GSE277116 is restricted to the preregistered main
18 tumour packages used by the scoreability sensitivity analysis.
"""
from __future__ import annotations

import csv
import importlib.util
import tarfile
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
LOCKED_SCRIPT = PROJECT / "scripts" / "run_unified_primary_pipeline_v2.py"
OUT = PROJECT / "data" / "03_results" / "program_scoreability_v3"
OUT.mkdir(parents=True, exist_ok=True)
PROGRAM_NAME = {"mregDC_strict": "mregDC_like"}


def load_locked_module():
    specification = importlib.util.spec_from_file_location("locked_primary_pipeline", LOCKED_SCRIPT)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not load locked scoring workflow: {LOCKED_SCRIPT}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def program_definitions(locked):
    primary = (
        ("DC_core", locked.FROZEN["DC_core"]["genes"]),
        ("mregDC_strict", locked.FROZEN["mregDC_strict"]["genes"]),
        ("Tfh_like", locked.FROZEN["Tfh_like"]["genes"]),
    )
    return primary + tuple(locked.COMPETITORS.items())


def scoreability_rows(cohort, sample, matrix, genes, locked):
    gene_index = {gene.upper(): index for index, gene in enumerate(genes)}
    scores = locked.make_scores(matrix, genes)
    rows = []
    for score_name, program_genes in program_definitions(locked):
        values = np.asarray(scores[score_name], dtype=float)
        available = sum(gene.upper() in gene_index for gene in program_genes)
        rows.append({
            "cohort": cohort,
            "sample": sample,
            "program": PROGRAM_NAME.get(score_name, score_name),
            "program_genes": ";".join(program_genes),
            "n_program_genes": len(program_genes),
            "n_genes_available": available,
            "gene_coverage": available / len(program_genes),
            "spot_nonzero_fraction": float(np.mean(values > 0)),
            "score_zero_fraction": float(np.mean(values == 0)),
        })
    return rows


def gse278687_rows(locked):
    rows = []
    for path in sorted(locked.G278.glob("GSM*_filtered_feature_bc_matrix.h5")):
        sample, matrix, genes, *_ = locked.read_g278(path)
        rows.extend(scoreability_rows("GSE278687", sample, matrix, genes, locked))
    samples = {row["sample"] for row in rows}
    if len(samples) != 21:
        raise RuntimeError(f"Expected 21 GSE278687 sections, found {len(samples)}")
    return rows


def gse277116_rows(locked):
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
            sample, matrix, genes, *_ = loaded
            rows.extend(scoreability_rows("GSE277116", sample, matrix, genes, locked))
    samples = {row["sample"] for row in rows}
    if len(samples) != 18:
        raise RuntimeError(f"Expected 18 GSE277116 main packages, found {len(samples)}")
    return rows


def summary_rows(rows, locked):
    summaries = []
    for cohort in sorted({row["cohort"] for row in rows}):
        for program in [name for name, _ in program_definitions(locked)]:
            published_name = PROGRAM_NAME.get(program, program)
            subset = [row for row in rows if row["cohort"] == cohort and row["program"] == published_name]
            summaries.append({
                "cohort": cohort,
                "program": published_name,
                "n_samples": len(subset),
                "median_gene_coverage": float(np.median([float(row["gene_coverage"]) for row in subset])),
                "median_spot_nonzero_fraction": float(np.median([float(row["spot_nonzero_fraction"]) for row in subset])),
                "median_score_zero_fraction": float(np.median([float(row["score_zero_fraction"]) for row in subset])),
            })
    return summaries


def main():
    locked = load_locked_module()
    program_order = {
        PROGRAM_NAME.get(name, name): index
        for index, (name, _) in enumerate(program_definitions(locked))
    }
    rows = gse278687_rows(locked) + gse277116_rows(locked)
    rows.sort(key=lambda row: (row["cohort"], row["sample"], program_order[row["program"]]))
    summaries = summary_rows(rows, locked)
    write_tsv(OUT / "program_scoreability_by_sample.tsv", rows)
    write_tsv(OUT / "program_scoreability_summary.tsv", summaries)
    print(f"Wrote {len(rows)} scoreability rows for {len({row['sample'] for row in rows})} samples.")


if __name__ == "__main__":
    main()
