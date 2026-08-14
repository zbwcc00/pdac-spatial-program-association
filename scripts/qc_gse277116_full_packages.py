"""Cohort-specific QC for the complete GSE277116 Visium packages.

Reads nested sample tarballs directly from the GEO outer archive. No sample
package is permanently unpacked. The output deliberately separates technical
quality from study-design eligibility: lymph-node samples remain QC controls,
but are not candidates for PDAC tumour replication.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import tarfile
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
import scipy.io
import matplotlib.pyplot as plt

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data/00_raw/spatial/GSE277116_RAW.tar"
MANIFEST = PROJECT / "data/03_results/GSE277116_audit/GSE277116_geo_sample_manifest.tsv"
OUT = PROJECT / "data/03_results/GSE277116_full_package_qc"
OUT.mkdir(parents=True, exist_ok=True)
FROZEN = json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8"))

PROGRAMS = {
    "DC_core": tuple(FROZEN["DC_core"]["genes"]),
    "mregDC_strict": tuple(FROZEN["mregDC_strict"]["genes"]),
    "Tfh_like": tuple(FROZEN["Tfh_like"]["genes"]),
    "B_cell_strict": ("MS4A1", "CD79A", "CD79B", "CD22", "CD37"),
    "Plasma_cell": tuple(FROZEN["Plasma_cell"]["genes"]),
}
ALL_MARKERS = sorted({g for genes in PROGRAMS.values() for g in genes})


def read_manifest():
    rows = {}
    with MANIFEST.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["gsm"] = row["gsm"].strip()
            row["tissue_class"] = "lymph_node" if row.get("tissue subtype", "").lower() == "lymph node" else "tumor"
            row["sample_type_norm"] = row.get("sample type", "").lower()
            rows[row["gsm"]] = row
    return rows


def gsm_from_nested(name: str) -> str:
    return name.split("_", 1)[0]


def member_by_basename(members, basename):
    return next((m for m in members if Path(m.name).name == basename), None)


def parse_coordinates(raw: bytes, filename: str):
    text = gzip.decompress(raw).decode("utf-8", errors="replace") if filename.endswith(".gz") else raw.decode("utf-8", errors="replace")
    rows = []
    lines = text.splitlines()
    has_header = bool(lines and lines[0].lower().startswith("barcode"))
    start = 1 if has_header else 0
    for line in lines[start:]:
        if not line.strip():
            continue
        fields = line.split(",")
        if len(fields) < 6:
            fields = line.split("\t")
        if len(fields) < 6:
            continue
        try:
            rows.append({
                "barcode": fields[0].strip(),
                "in_tissue": int(float(fields[1])),
                "array_row": int(float(fields[2])),
                "array_col": int(float(fields[3])),
                "pxl_row": float(fields[4]),
                "pxl_col": float(fields[5]),
            })
        except (TypeError, ValueError):
            continue
    return rows


def parse_features(raw: bytes):
    text = gzip.decompress(raw).decode("utf-8", errors="replace")
    genes = []
    for line in text.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        # 10x features are feature_id, gene_symbol, feature_type. Marker
        # programs use HGNC symbols; retain the symbol and fall back only for
        # non-standard one-column feature files.
        genes.append(fields[1] if len(fields) > 1 and fields[1] else fields[0])
    return genes


def connected_fraction(coords):
    """Largest 8-neighbour array-grid component / tissue spots."""
    points = {(r["array_row"], r["array_col"]) for r in coords if r["in_tissue"] == 1}
    if not points:
        return float("nan"), 0
    unseen = set(points)
    largest = 0
    while unseen:
        start = unseen.pop()
        q = deque([start])
        size = 1
        while q:
            r, c = q.popleft()
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    nb = (r + dr, c + dc)
                    if nb in unseen:
                        unseen.remove(nb)
                        q.append(nb)
                        size += 1
        largest = max(largest, size)
    return largest / len(points), largest


def safe_quantile(values, q):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.quantile(values, q)) if values.size else float("nan")


def qc_one(gsm, nested_name, raw_nested, metadata):
    result = {"gsm": gsm, "nested_sample": nested_name, **metadata}
    if not raw_nested:
        result["status"] = "incomplete"
        return result
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw_nested)), mode="r") as inner:
        members = [m for m in inner.getmembers() if m.isfile()]
        matrix_m = member_by_basename(members, "matrix.mtx.gz")
        barcode_m = member_by_basename(members, "barcodes.tsv.gz")
        feature_m = member_by_basename(members, "features.tsv.gz")
        coord_m = member_by_basename(members, "tissue_positions.csv") or member_by_basename(members, "tissue_positions_list.csv")
        result["matrix_present"] = int(matrix_m is not None)
        result["barcodes_present"] = int(barcode_m is not None)
        result["features_present"] = int(feature_m is not None)
        result["coordinates_present"] = int(coord_m is not None)
        result["images_present"] = int(any(Path(m.name).name in {"tissue_hires_image.png", "tissue_lowres_image.png"} for m in members))
        result["spatial_enrichment_present"] = int(any(Path(m.name).name == "spatial_enrichment.csv" for m in members))
        if not all((matrix_m, barcode_m, feature_m, coord_m)):
            result["status"] = "incomplete"
            return result
        barcodes = gzip.decompress(inner.extractfile(barcode_m).read()).decode("utf-8", errors="replace").splitlines()
        genes = parse_features(inner.extractfile(feature_m).read())
        coords = parse_coordinates(inner.extractfile(coord_m).read(), coord_m.name)
        coord_map = {r["barcode"]: r for r in coords}
        matrix = scipy.io.mmread(io.BytesIO(gzip.decompress(inner.extractfile(matrix_m).read()))).tocsr()
        if matrix.shape != (len(genes), len(barcodes)):
            result["status"] = "matrix_dimension_mismatch"
            result["matrix_rows"], result["matrix_cols"] = matrix.shape
            return result
        n = len(barcodes)
        tissue_idx = np.array([i for i, b in enumerate(barcodes) if coord_map.get(b, {}).get("in_tissue") == 1], dtype=int)
        coord_barcode_overlap = sum(1 for b in barcodes if b in coord_map)
        result.update({
            "status": "complete",
            "n_genes": len(genes),
            "n_barcodes": n,
            "n_coordinate_rows": len(coords),
            "n_coordinate_barcode_overlap": coord_barcode_overlap,
            "coordinate_barcode_overlap_fraction": coord_barcode_overlap / n if n else float("nan"),
            "n_tissue_spots": int(len(tissue_idx)),
            "tissue_fraction": len(tissue_idx) / n if n else float("nan"),
        })
        totals = np.asarray(matrix[:, tissue_idx].sum(axis=0)).ravel() if len(tissue_idx) else np.array([])
        detected = np.asarray((matrix[:, tissue_idx] > 0).sum(axis=0)).ravel() if len(tissue_idx) else np.array([])
        result.update({
            "median_umi_tissue": safe_quantile(totals, .5),
            "q10_umi_tissue": safe_quantile(totals, .1),
            "q90_umi_tissue": safe_quantile(totals, .9),
            "median_genes_tissue": safe_quantile(detected, .5),
            "q10_genes_tissue": safe_quantile(detected, .1),
            "q90_genes_tissue": safe_quantile(detected, .9),
            "zero_umi_tissue_fraction": float(np.mean(totals == 0)) if len(totals) else float("nan"),
        })
        gene_to_idx = {g: i for i, g in enumerate(genes)}
        marker_present = []
        for gene in ALL_MARKERS:
            gi = gene_to_idx.get(gene)
            nz = int((matrix[gi, tissue_idx] > 0).sum()) if gi is not None and len(tissue_idx) else 0
            frac = nz / len(tissue_idx) if len(tissue_idx) else float("nan")
            marker_present.append((gene, gi is not None, nz, frac))
            result[f"marker_{gene}_tissue_detection_fraction"] = frac
        result["n_requested_markers_in_features"] = int(sum(x[1] for x in marker_present))
        result["n_requested_markers_detected_in_tissue"] = int(sum(x[2] > 0 for x in marker_present))
        result["marker_feature_coverage_fraction"] = result["n_requested_markers_in_features"] / len(ALL_MARKERS)
        result["marker_tissue_detection_fraction"] = result["n_requested_markers_detected_in_tissue"] / len(ALL_MARKERS)
        for program, geneset in PROGRAMS.items():
            present = [(g, gene_to_idx.get(g)) for g in geneset]
            detected_program = [g for g, gi in present if gi is not None and len(tissue_idx) and (matrix[gi, tissue_idx] > 0).sum() > 0]
            one_percent = [g for g, gi in present if gi is not None and len(tissue_idx) and (matrix[gi, tissue_idx] > 0).sum() / len(tissue_idx) >= 0.01]
            result[f"{program}_genes_in_features"] = len([1 for _, gi in present if gi is not None])
            result[f"{program}_genes_detected_in_tissue"] = len(detected_program)
            result[f"{program}_genes_detected_at_1pct"] = len(one_percent)
            result[f"{program}_gene_coverage_fraction"] = len(detected_program) / len(geneset)
        frac, largest = connected_fraction(coords)
        result["largest_tissue_component_fraction"] = frac
        result["largest_tissue_component_spots"] = largest
        result["array_rows_tissue"] = len({coord_map[b]["array_row"] for b in barcodes if b in coord_map and coord_map[b]["in_tissue"] == 1})
        result["array_cols_tissue"] = len({coord_map[b]["array_col"] for b in barcodes if b in coord_map and coord_map[b]["in_tissue"] == 1})
        # Candidate flag is intentionally permissive and is not a biological filter.
        result["technical_candidate"] = int(
            result["coordinate_barcode_overlap_fraction"] >= 0.95
            and result["n_tissue_spots"] >= 100
            and result["median_umi_tissue"] >= 10
            and result["median_genes_tissue"] >= 10
            and result["largest_tissue_component_fraction"] >= 0.50
            and result["Tfh_like_gene_coverage_fraction"] >= 0.50
            and result["mregDC_strict_gene_coverage_fraction"] >= 0.50
        )
        result["program_scoreability_candidate"] = int(
            result["mregDC_strict_genes_detected_at_1pct"] >= 3
            and result["Tfh_like_genes_detected_at_1pct"] >= 5
        )
        result["tumor_replication_candidate"] = int(result["technical_candidate"] and result.get("tissue_class") == "tumor")
        result["tumor_program_scoreability_candidate"] = int(result["tumor_replication_candidate"] and result["program_scoreability_candidate"])
    return result


def main():
    metadata = read_manifest()
    rows = []
    with tarfile.open(RAW, "r") as outer:
        for member in outer.getmembers():
            if not member.name.endswith(".tar.gz"):
                continue
            gsm = gsm_from_nested(Path(member.name).name)
            rows.append(qc_one(gsm, Path(member.name).name.replace(".tar.gz", ""), outer.extractfile(member).read(), metadata.get(gsm, {"tissue_class": "unknown", "sample_type_norm": "unknown"})))
    rows.sort(key=lambda r: r["gsm"])
    all_fields = sorted({key for row in rows for key in row})
    with (OUT / "GSE277116_package_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t", extrasaction="ignore"); writer.writeheader(); writer.writerows(rows)
    marker_rows = []
    for row in rows:
        if row.get("status") != "complete":
            continue
        for program, geneset in PROGRAMS.items():
            for gene in geneset:
                marker_rows.append({
                    "gsm": row["gsm"], "title": row.get("title"), "tissue_class": row.get("tissue_class"),
                    "sample_type_norm": row.get("sample_type_norm"), "program": program, "gene": gene,
                    "tissue_detection_fraction": row.get(f"marker_{gene}_tissue_detection_fraction"),
                })
    with (OUT / "GSE277116_marker_detection_by_sample.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = list(marker_rows[0]); writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(marker_rows)
    with (OUT / "GSE277116_tumor_replication_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["gsm", "title", "tissue_class", "sample_type_norm", "tissue subtype", "treatment arm", "technical_candidate", "program_scoreability_candidate", "tumor_replication_candidate", "tumor_program_scoreability_candidate", "n_tissue_spots", "median_umi_tissue", "median_genes_tissue", "Tfh_like_gene_coverage_fraction", "Tfh_like_genes_detected_at_1pct", "mregDC_strict_gene_coverage_fraction", "mregDC_strict_genes_detected_at_1pct", "B_cell_strict_gene_coverage_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"); writer.writeheader(); writer.writerows(r for r in rows if r.get("tissue_class") == "tumor" and r.get("status") == "complete")
    complete = [r for r in rows if r.get("status") == "complete"]
    summary = [
        "- **package_status**: " + "; ".join(f"{k}={v}" for k, v in sorted(Counter(str(r.get("status", "NA")) for r in rows).items())),
    ]
    for key in ["tissue_class", "sample_type_norm", "technical_candidate", "program_scoreability_candidate", "tumor_replication_candidate", "tumor_program_scoreability_candidate"]:
        counts = Counter(str(r.get(key, "NA")) for r in complete)
        summary.append(f"- **complete_packages_{key}**: " + "; ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    candidates = [r for r in rows if r.get("tumor_replication_candidate") == 1]
    scoreable = [r for r in rows if r.get("tumor_program_scoreability_candidate") == 1]
    med_umi = safe_quantile([r["median_umi_tissue"] for r in complete], .5)
    med_genes = safe_quantile([r["median_genes_tissue"] for r in complete], .5)
    report = [
        "# GSE277116 full-package spatial QC", "",
        "## Scope and design", "",
        "This QC reads all 28 nested Visium packages directly from the GEO outer archive. Two packages lacking expression matrices/barcodes/features are excluded before quantitative QC. Among the 26 complete packages, GEO metadata identifies 21 tumour samples (8 FFPE and 13 frozen) and 5 lymph-node samples (4 FFPE and 1 frozen). No patient identifier is available; samples remain sample-level units.", "",
        "## Technical summary", "", *summary, "",
        f"- Complete packages quantitatively assessed: **{len(complete)}**.",
        f"- Candidate tumour replication packages under the prespecified technical screen: **{len(candidates)}**.",
        f"- Tumour packages meeting the additional marker-detection scoreability screen: **{len(scoreable)}**.",
        f"- Across complete packages, median of per-package median tissue UMI: **{med_umi:.1f}**; median of per-package median detected genes: **{med_genes:.1f}**.", "",
        "## Candidate screen", "",
        "The technical-candidate flag requires coordinate/barcode overlap >=0.95, >=100 tissue spots, median tissue UMI >=10, median detected genes >=10, largest 8-neighbour tissue component >=50% of tissue spots, and at least 50% of Tfh-like and strict mregDC-like marker genes detected in tissue. The separate program-scoreability flag requires at least 3/6 strict mregDC-like and 5/9 Tfh-like markers to be detected in >=1% of tissue spots. These are permissive operational gates for deciding which packages merit program-level replication QC; they are not biological truth thresholds and must be reviewed by FFPE/frozen stratum.", "",
        "Only packages labeled Tumor/tumor are eligible for the replication manifest. Lymph-node packages are retained as anatomical controls and are excluded from the PDAC tumour replication endpoint. The CODA annotation resource covers only the J1568 FFPE subset and has no positive TLS labels, so annotation cannot be used as a TLS-positive endpoint.", "",
        "## Files", "",
        "- `GSE277116_package_qc.tsv`: one row per nested package with expression, coordinate, complexity, marker and spatial-connectivity metrics.",
        "- `GSE277116_marker_detection_by_sample.tsv`: per-sample tissue detection fractions for every frozen-program marker.",
        "- `GSE277116_tumor_replication_manifest.tsv`: tumour-only sample-level candidate manifest.",
        "- `GSE277116_full_package_qc_report.md`: this audit and interpretation.",
    ]
    (OUT / "GSE277116_full_package_qc_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    plot_qc(rows)
    print(json.dumps({"packages": len(rows), "complete": len(complete), "tumor_candidates": len(candidates), "tumor_scoreable": len(scoreable), "output": str(OUT)}, ensure_ascii=False, indent=2))


def plot_qc(rows):
    complete = [r for r in rows if r.get("status") == "complete"]
    tissue_order = {"tumor": 0, "lymph_node": 1}
    ordered = sorted(complete, key=lambda r: (tissue_order.get(r.get("tissue_class"), 9), r.get("sample_type_norm"), r.get("title")))
    labels = [r.get("title", r["gsm"]) for r in ordered]
    colors = ["#2A9D8F" if r.get("sample_type_norm") == "ffpe" else "#457B9D" for r in ordered]
    tumour_count = sum(r.get("tissue_class") == "tumor" for r in ordered)
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), constrained_layout=True)
    fig.suptitle("GSE277116 full-package spatial QC: teal = FFPE; blue = frozen; tumour packages precede lymph-node controls", fontsize=13)
    for ax, metric, ylabel, logy in [
        (axes[0, 0], "n_tissue_spots", "Tissue spots", False),
        (axes[0, 1], "median_umi_tissue", "Median UMI per tissue spot", True),
        (axes[1, 0], "median_genes_tissue", "Median detected genes per tissue spot", True),
    ]:
        ax.bar(range(len(ordered)), [float(r[metric]) for r in ordered], color=colors, width=.8)
        ax.set_xticks(range(len(ordered)), labels, rotation=65, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        if logy: ax.set_yscale("log")
        ax.grid(axis="y", alpha=.25)
    matrix = np.array([[float(r.get(f"{p}_genes_detected_at_1pct", 0)) / len(g) for r in ordered] for p, g in PROGRAMS.items()])
    im = axes[1, 1].imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    axes[1, 1].set_yticks(range(len(PROGRAMS)), list(PROGRAMS))
    axes[1, 1].set_xticks(range(len(ordered)), labels, rotation=65, ha="right", fontsize=8)
    axes[1, 1].set_title("Markers detected in >=1% of tissue spots")
    fig.colorbar(im, ax=axes[1, 1], label="Program marker fraction")
    for ax in axes.flat:
        ax.axvline(tumour_count - .5, color="#555555", lw=.8, ls="--", alpha=.7)
    fig.savefig(OUT / "GSE277116_full_package_qc_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
