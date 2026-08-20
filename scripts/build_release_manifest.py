"""Create a content-addressed manifest for the public source-and-results release."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "release_manifest.json"
EXCLUDED_ROOTS = {".git", "data/00_raw", "data/01_unpacked", "data/02_intermediate"}
EXCLUDED_PATH_PARTS = {"__pycache__"}
EXCLUDED_NAMES = {
    "release_manifest.json",
    "release_verification.json",
    "download_manifest.json",
    "local_release_manifest_v7.json",
    "verification_manifest_v7.json",
    "release_verification_v7.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def include(path: Path) -> bool:
    relative = path.relative_to(PROJECT).as_posix()
    return (
        path.name not in EXCLUDED_NAMES
        and path.suffix.lower() not in {".pyc", ".pyo"}
        and not any(relative == root or relative.startswith(root + "/") for root in EXCLUDED_ROOTS)
        and not (EXCLUDED_PATH_PARTS & set(Path(relative).parts))
        and not any(part.startswith("qa_") for part in relative.split("/"))
    )


def main() -> None:
    records = []
    for path in sorted(PROJECT.rglob("*")):
        if path.is_file() and include(path):
            records.append({"path": path.relative_to(PROJECT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    payload = {
        "schema": "pdac-spatial-program-association-public-release-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "public code, manuscript sources, derived results, and figures; raw GEO data excluded",
        "input_config": "config/public_inputs.json",
        "files": records,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8", newline="\n")
    print(f"Wrote {OUTPUT} with {len(records)} records")


if __name__ == "__main__":
    main()
