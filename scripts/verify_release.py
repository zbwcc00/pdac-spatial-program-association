"""Verify the public release manifest and optionally require local GEO inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "release_manifest.json"
INPUTS = PROJECT / "config" / "public_inputs.json"
REPORT = PROJECT / "release_verification.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-inputs", action="store_true", help="also verify that all raw GEO inputs are locally present")
    args = parser.parse_args()
    if not MANIFEST.is_file():
        raise FileNotFoundError("Build release_manifest.json first.")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures, checked = [], []
    for item in manifest["files"]:
        path = PROJECT / item["path"]
        if not path.is_file():
            failures.append(f"missing: {item['path']}")
        elif path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            failures.append(f"hash mismatch: {item['path']}")
        else:
            checked.append(item["path"])
    if args.require_inputs:
        for item in json.loads(INPUTS.read_text(encoding="utf-8"))["inputs"]:
            path = PROJECT / item["destination"]
            if not path.is_file():
                failures.append(f"missing public input: {item['id']}")
            elif "sha256" in item and sha256(path) != item["sha256"]:
                failures.append(f"raw hash mismatch: {item['id']}")
            else:
                checked.append(item["destination"])
    payload = {
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": MANIFEST.name,
        "require_inputs": args.require_inputs,
        "checked_count": len(checked),
        "failure_count": len(failures),
        "failures": failures,
        "status": "PASS" if not failures else "FAIL",
    }
    REPORT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
