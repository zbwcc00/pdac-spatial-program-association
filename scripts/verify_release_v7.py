"""Verify the v7 local release manifest, raw input, and manuscript boundaries."""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "verification_manifest_v7.json"
REPORT = PROJECT / "release_verification_v7.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if not MANIFEST.is_file():
        raise FileNotFoundError(f"Build the manifest first: {MANIFEST}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures: list[str] = []
    checked: list[str] = []
    for item in manifest["content_addressed_files"]:
        path = PROJECT / item["path"]
        if not path.is_file():
            failures.append(f"missing: {item['path']}")
            continue
        observed = sha256(path)
        checked.append(item["path"])
        if observed != item["sha256"]:
            failures.append(f"hash mismatch: {item['path']}")
    for item in manifest["raw_inputs"]:
        path = PROJECT / item["path"]
        if not path.is_file():
            failures.append(f"missing raw input: {item['path']}")
            continue
        if path.stat().st_size != item["bytes"]:
            failures.append(f"raw byte mismatch: {item['path']}")
        elif sha256(path) != item["sha256"]:
            failures.append(f"raw hash mismatch: {item['path']}")
        else:
            checked.append(item["path"])
    manuscript = (PROJECT / "manuscript/drafts/Manuscript_full_v7_single_cell_attribution.md").read_text(encoding="utf-8")
    required = [
        "GSE202051",
        "two-sided paired Wilcoxon signed-rank tests",
        "no author-defined Tfh label",
        "not the spatial association, cell-cell contact or mechanism",
    "Supplementary Figure S9 and Tables S30-S33",
    ]
    forbidden = [
        "one-sided paired Wilcoxon",
        "single-cell validation of the spatial association",
        "mregDC-Tfh niche",
        "Tfh-specific spatial result",
    ]
    for text in required:
        if text not in manuscript:
            failures.append(f"manuscript missing required boundary: {text}")
    for text in forbidden:
        if text in manuscript:
            failures.append(f"manuscript contains prohibited claim: {text}")
    s21 = (PROJECT / "manuscript/supplementary_v7_single_cell_attribution/Supplementary_Table_S21_GEO_accession_audit.tsv").read_text(encoding="utf-8")
    if "GSE202051" not in s21:
        failures.append("Supplementary Table S21 omits GSE202051")
    report = {
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(MANIFEST.relative_to(PROJECT)),
        "checked_files": checked,
        "failure_count": len(failures),
        "failures": failures,
        "status": "PASS" if not failures else "FAIL",
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
