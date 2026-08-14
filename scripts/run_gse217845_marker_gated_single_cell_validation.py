"""PDAC single-cell program-detectability QC with transparent marker gates.

GSE217845 (a human-PDAC subseries of GSE217847; Caronni et al., Nature 2023) provides raw matrices but no public
cell-type metadata in GEO.  This script therefore does *not* claim author
annotations. It labels only high-confidence marker-gated candidate states and
reports unassigned cells explicitly.  The purpose is to test program
localization, not estimate state prevalence or discover new signatures.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import tarfile
from pathlib import Path

import numpy as np
import scipy.io

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data/00_raw/spatial/GSE217847_RAW.tar"
OUT = PROJECT / "data/03_results/GSE217845_marker_gated_single_cell_validation"
OUT.mkdir(parents=True, exist_ok=True)
FROZEN = json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8"))

MREG = tuple(FROZEN["mregDC_strict"]["genes"])
TFH = tuple(FROZEN["Tfh_like"]["genes"])
GATES = {
    "mregDC_candidate": {"required_any": ("FCER1A", "CD1C", "CLEC10A"), "minimum_detected": 1, "program": MREG, "program_minimum_detected": 2},
    "Tfh_candidate": {"required_any": ("CD3D", "CD3E", "TRAC"), "minimum_detected": 1, "program": TFH, "program_minimum_detected": 3, "required_program_any": ("CXCR5", "PDCD1", "ICOS", "BCL6", "CXCL13")},
    "non_Tfh_CD4_candidate": {"required_all": ("CD3D", "CD3E"), "required_any": ("IL7R", "LTB", "CCR7", "TCF7"), "minimum_detected": 1, "exclude_any": ("CXCR5", "PDCD1", "ICOS", "BCL6", "CXCL13")},
    "Treg_candidate": {"required_any": ("FOXP3", "IL2RA", "CTLA4", "TIGIT"), "minimum_detected": 2},
    "exhausted_CD8_candidate": {"required_any": ("CD8A", "CD8B"), "minimum_detected": 1, "program": ("LAG3", "HAVCR2", "ENTPD1", "LAYN", "TOX"), "program_minimum_detected": 2},
    "macrophage_candidate": {"required_any": ("C1QA", "C1QB", "C1QC", "APOC1", "SPP1"), "minimum_detected": 2},
}


def parse_member(archive, members, suffix):
    return next((member for member in members if member.name.endswith(suffix)), None)


def read_matrix(archive, prefix):
    members = [member for member in archive.getmembers() if member.isfile() and member.name.startswith(prefix)]
    matrix_member = parse_member(archive, members, "_matrix.mtx.gz")
    feature_member = parse_member(archive, members, "_features.tsv.gz")
    barcode_member = parse_member(archive, members, "_barcodes.tsv.gz")
    if not all((matrix_member, feature_member, barcode_member)):
        return None
    matrix = scipy.io.mmread(io.BytesIO(gzip.decompress(archive.extractfile(matrix_member).read()))).tocsr()
    features = gzip.decompress(archive.extractfile(feature_member).read()).decode("utf-8", errors="replace").splitlines()
    genes = [(line.split("\t")[1] if len(line.split("\t")) > 1 and line.split("\t")[1] else line.split("\t")[0]).upper() for line in features]
    barcodes = gzip.decompress(archive.extractfile(barcode_member).read()).decode("utf-8", errors="replace").splitlines()
    return matrix, genes, barcodes


def detect(matrix, gene_index, genes):
    return np.asarray(matrix[[gene_index[g] for g in genes if g in gene_index], :].sum(axis=0)).ravel() > 0


def number_detected(matrix, gene_index, genes):
    index = [gene_index[g] for g in genes if g in gene_index]
    return np.asarray((matrix[index, :] > 0).sum(axis=0)).ravel() if index else np.zeros(matrix.shape[1], dtype=int)


def score(matrix, gene_index, genes):
    index = [gene_index[g] for g in genes if g in gene_index]
    return np.log1p(np.asarray(matrix[index, :].sum(axis=0)).ravel()) if index else np.zeros(matrix.shape[1])


def gate(matrix, gene_index, rule):
    result = np.ones(matrix.shape[1], dtype=bool)
    if "required_all" in rule:
        result &= np.all([(np.asarray(matrix[gene_index[g], :].todense()).ravel() > 0) if g in gene_index else np.zeros(matrix.shape[1], dtype=bool) for g in rule["required_all"]], axis=0)
    if "required_any" in rule:
        result &= number_detected(matrix, gene_index, rule["required_any"]) >= rule.get("minimum_detected", 1)
    if "program" in rule:
        result &= number_detected(matrix, gene_index, rule["program"]) >= rule["program_minimum_detected"]
    if "required_program_any" in rule:
        result &= number_detected(matrix, gene_index, rule["required_program_any"]) >= 1
    if "exclude_any" in rule:
        result &= number_detected(matrix, gene_index, rule["exclude_any"]) == 0
    return result


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    sample_rows, cell_rows = [], []
    with tarfile.open(RAW, "r") as archive:
        human_tumour_gsms = {f"GSM67275{i}" for i in range(42, 52)}
        prefixes = sorted({member.name.split("_barcodes.tsv.gz")[0] for member in archive.getmembers() if member.name.split("_", 1)[0] in human_tumour_gsms and member.name.endswith("_barcodes.tsv.gz") and "_Tumor_" in member.name})
        # GSE217845 is a subseries within the multi-species GSE217847 archive.
        # The accession list in GEO identifies GSM6727542--GSM6727551 as the
        # ten human resected PDAC tumour matrices; this explicit list excludes
        # all peripheral blood, mouse and spatial samples.
        for prefix in prefixes:
            parsed = read_matrix(archive, prefix)
            if parsed is None:
                continue
            matrix, genes, barcodes = parsed
            gene_index = {gene: i for i, gene in enumerate(genes)}
            calls = {name: gate(matrix, gene_index, rule) for name, rule in GATES.items()}
            # Priority avoids dual labels while preserving the constituent gates.
            label = np.full(len(barcodes), "unassigned", dtype=object)
            for name in ("mregDC_candidate", "Tfh_candidate", "Treg_candidate", "exhausted_CD8_candidate", "non_Tfh_CD4_candidate", "macrophage_candidate"):
                label[(label == "unassigned") & calls[name]] = name
            sample = prefix.split("_", 1)[1].rsplit("_", 1)[0]
            row = {"sample": sample, "gsm": prefix.split("_", 1)[0], "cells": len(barcodes)}
            for name, value in calls.items():
                row[name] = int(value.sum())
                row[f"{name}_fraction"] = float(value.mean())
            for program, geneset in {"mregDC_strict": MREG, "Tfh_like": TFH}.items():
                row[f"{program}_gene_coverage"] = len([gene for gene in geneset if gene in gene_index]) / len(geneset)
                row[f"{program}_score_zero_fraction"] = float((score(matrix, gene_index, geneset) == 0).mean())
            sample_rows.append(row)
            for name in ("mregDC_candidate", "Tfh_candidate", "non_Tfh_CD4_candidate", "Treg_candidate", "exhausted_CD8_candidate", "macrophage_candidate", "unassigned"):
                mask = label == name
                if not mask.any():
                    continue
                cell_rows.append({"sample": sample, "label": name, "cells": int(mask.sum()), "median_mregDC_score": float(np.median(score(matrix, gene_index, MREG)[mask])), "median_Tfh_score": float(np.median(score(matrix, gene_index, TFH)[mask])), "mregDC_score_positive_fraction": float((score(matrix, gene_index, MREG)[mask] > 0).mean()), "Tfh_score_positive_fraction": float((score(matrix, gene_index, TFH)[mask] > 0).mean())})
    write_tsv(OUT / "GSE217845_marker_gated_sample_summary.tsv", sample_rows)
    write_tsv(OUT / "GSE217845_marker_gated_state_scores.tsv", cell_rows)
    (OUT / "GSE217845_marker_gating_protocol.json").write_text(json.dumps({"dataset": "GSE217845", "source": "GEO family SOFT: human-PDAC subseries of GSE217847; raw matrices traced to GSM6727542-GSM6727551 in GSE217847_RAW.tar; cell-type metadata not publicly supplied in GEO", "purpose": "marker-gated program-detectability and gating-feasibility QC", "gates": GATES, "priority_order": ["mregDC_candidate", "Tfh_candidate", "Treg_candidate", "exhausted_CD8_candidate", "non_Tfh_CD4_candidate", "macrophage_candidate"], "interpretation_limit": "candidate gate labels are not author annotations and do not establish canonical cell identity"}, indent=2), encoding="utf-8")
    print(json.dumps({"tumour_samples": len(sample_rows), "output": str(OUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
