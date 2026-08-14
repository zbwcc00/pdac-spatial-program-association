"""Forensic eligibility audit for the public GSE202051 processed objects.

This script does not calculate biological program scores.  It verifies whether
the locally available H5AD objects can serve as an independent, author-labelled
PDAC reference under the manuscript's documented eligibility rules.
"""
from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

import h5py


PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "00_raw" / "single_cell" / "GSE202051"
OUT = PROJECT / "data" / "03_results" / "GSE202051_author_reference_audit"
OUT.mkdir(parents=True, exist_ok=True)
REMOTE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202051/suppl/"
    "GSE202051_totaldata-final-toshare.h5ad.gz"
)
TOTAL_NAME = "GSE202051_totaldata-final-toshare.h5ad.gz"
FINAL_LABEL_FIELDS = {
    "celltypes",
    "cell_subsets",
    "detailed_cell_subsets",
    "annot_level_1",
    "annot_level_2",
    "new_celltypes",
    "Level 1 Annotation",
    "Level 2 Annotation",
    "Level 3 Annotation",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def remote_size() -> int | None:
    request = urllib.request.Request(REMOTE_URL, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            value = response.headers.get("Content-Length")
            return int(value) if value else None
    except Exception as error:  # Network availability must not hide local facts.
        print(f"Remote HEAD unavailable: {error}")
        return None


def h5ad_structure(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        var = handle["var"]
        # H5AD 0.7 stores a compound obs dataset; newer files store obs/var as
        # groups. Supporting both prevents a format change being mistaken for
        # absent author annotations.
        obs_fields = list(obs.keys()) if hasattr(obs, "keys") else list(obs.dtype.names or [])
        var_fields = list(var.keys()) if hasattr(var, "keys") else list(var.dtype.names or [])
        x = handle["X"]
        if hasattr(x, "keys"):
            shape = [
                int(x.attrs["shape"][0]) if "shape" in x.attrs else len(handle["obs"]),
                int(x.attrs["shape"][1]) if "shape" in x.attrs else len(handle["var"]),
            ]
            storage = "sparse_group"
        else:
            shape = list(x.shape)
            storage = "dense_dataset"
        label_fields = sorted(FINAL_LABEL_FIELDS.intersection(obs_fields))
        pid_categories = []
        if "pid_categories" in handle.get("uns", {}):
            pid_categories = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["uns"]["pid_categories"][:]]
        provenance_fields = [field for field in ("pid", "sampleid", "treatment_status", "new_treatment") if field in obs_fields]
        return {
            "file": path.name,
            "sha256": sha256(path),
            "bytes": path.stat().st_size,
            "shape_cells_by_genes": shape,
            "expression_storage": storage,
            "obs_fields": obs_fields,
            "var_fields": var_fields,
            "author_final_label_fields": label_fields,
            "provenance_fields": provenance_fields,
            "has_patient_or_sample_field": bool({"pid", "sampleid"}.intersection(obs_fields)),
            "pid_categories": pid_categories,
            "has_leiden_field": "leiden" in obs_fields,
            "eligibility": (
                "ineligible: no final author cell-type/state field in this object"
                if not label_fields
                else "eligible in principle: requires barcode and provenance audit before scoring"
            ),
        }


def main() -> None:
    expected_bytes = remote_size()
    total_path = RAW / TOTAL_NAME
    total_bytes = total_path.stat().st_size if total_path.exists() else 0
    total_complete = bool(expected_bytes and total_bytes == expected_bytes)
    objects = [h5ad_structure(path) for path in sorted(RAW.glob("*.h5ad"))]
    integrated = next((item for item in objects if item["file"] == TOTAL_NAME.removesuffix(".gz")), None)
    integrated_eligible = bool(
        total_complete
        and integrated
        and integrated["author_final_label_fields"]
        and integrated["has_patient_or_sample_field"]
    )
    result = {
        "dataset": "GSE202051",
        "source": "Hwang et al., Nature Genetics 2022, PMID 35902743, DOI 10.1038/s41588-022-01134-8",
        "purpose": "eligibility audit for author-labelled program-attribution/specificity QC only",
        "remote_total_object": {
            "url": REMOTE_URL,
            "expected_bytes_from_head": expected_bytes,
            "local_file": TOTAL_NAME,
            "local_bytes": total_bytes,
            "local_sha256": sha256(total_path) if total_path.exists() else None,
            "complete_by_exact_length": total_complete,
        },
        "local_h5ad_objects": objects,
        "decision": (
            "Eligible for fixed-program attribution/specificity QC: the complete integrated object carries final author labels and provenance fields."
            if integrated_eligible
            else "DO NOT SCORE: the required integrated object is incomplete or does not provide both final author labels and provenance."
        ),
        "interpretation_boundary": (
            "This object can support only program attribution/specificity QC. Do not call its results single-cell spatial validation, "
            "a cell interaction, TLS evidence, or mechanism."
        ),
    }
    (OUT / "gse202051_author_reference_audit.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    lines = [
        "# GSE202051 author-reference audit",
        "",
        f"- Remote expected bytes: `{expected_bytes}`.",
        f"- Local integrated object bytes: `{total_bytes}`.",
        f"- Exact-length complete: `{total_complete}`.",
        f"- Decision: {result['decision']}",
        "",
        "## Individual local objects",
        "",
    ]
    for item in objects:
        lines.extend([
            f"### {item['file']}",
            f"- Cells x genes: `{item['shape_cells_by_genes'][0]} x {item['shape_cells_by_genes'][1]}`.",
            f"- Final author label fields: `{', '.join(item['author_final_label_fields']) or 'none'}`.",
            f"- Provenance fields: `{', '.join(item['provenance_fields']) or 'none'}`; PID categories: `{', '.join(item['pid_categories']) or 'none'}`.",
            f"- Eligibility: {item['eligibility']}",
            "",
        ])
    lines.append("The audit is a data-eligibility check only; no biological inference has been calculated.")
    (OUT / "gse202051_author_reference_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "decision": result["decision"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
