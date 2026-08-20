"""Build a public source-and-results ZIP strictly from release_manifest.json."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "release_manifest.json"
VERIFICATION = PROJECT / "release_verification.json"
DEFAULT_OUTPUT = PROJECT.parent / "pdac-spatial-program-association-v1.0.6-source-and-results.zip"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if not MANIFEST.is_file() or not VERIFICATION.is_file():
        raise FileNotFoundError("Build and verify release_manifest.json before creating the ZIP.")
    verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
    if verification.get("status") != "PASS":
        raise RuntimeError("Refusing to archive a release whose verification status is not PASS.")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    members = [item["path"] for item in manifest["files"]] + [MANIFEST.name, VERIFICATION.name]
    if len(members) != len(set(members)):
        raise RuntimeError("Release ZIP member list contains duplicate paths.")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    root = PROJECT.name
    with ZipFile(output, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for relative in members:
            source = PROJECT / relative
            if not source.is_file():
                raise FileNotFoundError(f"Manifest member is missing: {relative}")
            archive.write(source, arcname=f"{root}/{Path(relative).as_posix()}")
    print(json.dumps({"output": str(output), "member_count": len(members), "root": root}, ensure_ascii=False))


if __name__ == "__main__":
    main()
