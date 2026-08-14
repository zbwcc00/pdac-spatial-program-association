"""Download, verify, and prepare the public GEO inputs for this release.

Raw inputs are intentionally local-only.  The script performs no biological
analysis; it records provenance in an ignored download manifest and verifies
the published SHA-256 for the large GSE202051 archive.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT / "config" / "public_inputs.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    temporary = destination.with_name(destination.name + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {destination}")
    request = urllib.request.Request(url, headers={"User-Agent": "pdac-spatial-program-association/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
    temporary.replace(destination)


def verify(item: dict, path: Path) -> dict:
    result = {"id": item["id"], "path": str(path.relative_to(PROJECT)), "exists": path.is_file()}
    if not path.is_file():
        return result
    result["bytes"] = path.stat().st_size
    if "bytes" in item:
        result["expected_bytes"] = item["bytes"]
        result["bytes_match"] = result["bytes"] == item["bytes"]
    if "sha256" in item:
        result["sha256"] = sha256(path)
        result["expected_sha256"] = item["sha256"]
        result["sha256_match"] = result["sha256"] == item["sha256"]
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="download missing selected inputs")
    parser.add_argument("--force", action="store_true", help="re-download selected inputs even if present")
    parser.add_argument("--only", action="append", choices=("GSE278687", "GSE277116", "GSE217847", "GSE202051"), help="restrict to one or more accession IDs")
    parser.add_argument("--no-decompress", action="store_true", help="do not decompress GSE202051 after verification")
    parser.add_argument("--list", action="store_true", help="list registered inputs without checking files")
    args = parser.parse_args()

    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    selected = [item for item in config["inputs"] if not args.only or item["id"] in args.only]
    if args.list:
        for item in selected:
            print(f"{item['id']}: {item['role']}\n  {item['url']}\n  {item['destination']}")
        return

    results = []
    for item in selected:
        destination = PROJECT / item["destination"]
        if args.download and (args.force or not destination.is_file()):
            download(item["url"], destination)
        result = verify(item, destination)
        if result["exists"] and "sha256_match" in result and not result["sha256_match"]:
            raise RuntimeError(f"SHA-256 mismatch for {item['id']}: {destination}")
        if result["exists"] and "bytes_match" in result and not result["bytes_match"]:
            raise RuntimeError(f"byte-count mismatch for {item['id']}: {destination}")
        if result["exists"] and "decompress_to" in item and not args.no_decompress:
            expanded = PROJECT / item["decompress_to"]
            if args.force or not expanded.is_file():
                print(f"Decompressing {destination.name} -> {expanded.name}")
                expanded.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(destination, "rb") as source, expanded.open("wb") as output:
                    shutil.copyfileobj(source, output, length=8 * 1024 * 1024)
            result["decompressed_path"] = str(expanded.relative_to(PROJECT))
            result["decompressed_exists"] = expanded.is_file()
        results.append(result)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(CONFIG.relative_to(PROJECT)),
        "results": results,
    }
    path = PROJECT / "data" / "00_raw" / "download_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    missing = [item["id"] for item in results if not item["exists"]]
    if missing:
        print("Missing required inputs: " + ", ".join(missing), file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
