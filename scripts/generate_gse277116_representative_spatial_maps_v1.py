"""Generate prespecified GSE277116 representative spatial program maps."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import json
import tarfile
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT / "scripts" / "run_unified_primary_pipeline_v2.py"
OUT = PROJECT / "data" / "03_results" / "GSE277116_spatial_maps_and_robustness"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42})


def load_locked_module():
    specification = importlib.util.spec_from_file_location("locked_primary_pipeline", PIPELINE)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not load locked scoring workflow: {PIPELINE}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def read_tsv(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def choose_representatives(locked):
    result_rows = read_tsv(locked.OUT / "GSE277116_v2_per_package.tsv")
    selected = []
    for stratum in ("ffpe", "frozen"):
        subset = [row for row in result_rows if row["analysis_set"] == "main_18" and row["stratum"] == stratum]
        median = float(np.median([float(row["primary_mregDC_Tfh_local_r"]) for row in subset]))
        choice = min(subset, key=lambda row: (abs(float(row["primary_mregDC_Tfh_local_r"]) - median), row["gsm"]))
        selected.append({
            "gsm": choice["gsm"],
            "sample": choice["sample"],
            "stratum": stratum,
            "mregdc_tfh_local_r": float(choice["primary_mregDC_Tfh_local_r"]),
            "selection_rule": "closest to primary-set stratum median local effect",
        })
    return selected


def nested_members(locked, gsm):
    with tarfile.open(locked.G277_RAW, "r") as outer:
        member = next((item for item in outer.getmembers() if Path(item.name).name.startswith(f"{gsm}_") and item.name.endswith(".tar.gz")), None)
        if member is None:
            raise RuntimeError(f"Could not locate {gsm} in {locked.G277_RAW}")
        return outer.extractfile(member).read()


def image_coordinates(raw_nested, locked):
    with tarfile.open(fileobj=io.BytesIO(gzip.decompress(raw_nested)), mode="r") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        barcode_member = locked.member_by_basename(members, "barcodes.tsv.gz")
        position_member = locked.member_by_basename(members, "tissue_positions.csv") or locked.member_by_basename(members, "tissue_positions_list.csv")
        image_member = locked.member_by_basename(members, "tissue_hires_image.png")
        scale_member = locked.member_by_basename(members, "scalefactors_json.json")
        if not all((barcode_member, position_member, image_member, scale_member)):
            raise RuntimeError("Representative package lacks an image, scale factors, barcodes, or positions.")
        barcodes = gzip.decompress(archive.extractfile(barcode_member).read()).decode("utf-8", errors="replace").splitlines()
        raw_positions = archive.extractfile(position_member).read().decode("utf-8", errors="replace").splitlines()
        pixel_positions = {}
        for line in raw_positions[1:]:
            fields = line.replace("\t", ",").split(",")
            if len(fields) >= 6:
                try:
                    if int(float(fields[1])) == 1:
                        pixel_positions[fields[0].strip()] = (float(fields[5]), float(fields[4]))
                except ValueError:
                    continue
        image = np.asarray(Image.open(io.BytesIO(archive.extractfile(image_member).read())).convert("RGB"))
        scale = float(json.loads(archive.extractfile(scale_member).read().decode("utf-8"))["tissue_hires_scalef"])
    keep = [barcode for barcode in barcodes if barcode in pixel_positions]
    return image, np.asarray([pixel_positions[barcode] for barcode in keep], dtype=float) * scale


def fields_for_package(raw_nested, gsm, locked):
    loaded = locked.read_g277(gsm, raw_nested, locked.g277_manifest()[gsm])
    if loaded is None:
        raise RuntimeError(f"Could not parse eligible package {gsm}")
    _, matrix, genes, coords, _, _, _ = loaded
    scores = locked.add_mreg_no_ccl19(locked.make_scores(matrix, genes), matrix, genes)
    mreg, tfh, mask, neighbours, _ = locked.prepared_fields(scores, coords)
    return scores["DC_core"], locked.local_field(locked.z(mreg), neighbours), locked.local_field(locked.z(tfh), neighbours), mask


def draw_panel(axis, image, coordinates, values, title, cmap, norm=None, mask=None):
    axis.imshow(image, origin="upper")
    keep = np.ones(len(coordinates), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    scatter = axis.scatter(coordinates[keep, 0], coordinates[keep, 1], c=np.asarray(values)[keep], s=6, cmap=cmap, norm=norm, linewidths=0, alpha=.88)
    axis.set_title(title, fontsize=8.4, weight="bold", pad=4)
    axis.set_axis_off()
    return scatter


def main():
    locked = load_locked_module()
    selections = choose_representatives(locked)
    prepared = []
    for selection in selections:
        raw_nested = nested_members(locked, selection["gsm"])
        image, pixel_coordinates = image_coordinates(raw_nested, locked)
        dc_core, local_mreg, local_tfh, mask = fields_for_package(raw_nested, selection["gsm"], locked)
        if len(pixel_coordinates) != len(local_mreg):
            raise RuntimeError(f"Coordinate alignment failed for {selection['gsm']}")
        prepared.append((selection, image, pixel_coordinates, dc_core, local_mreg, local_tfh, mask))
    limit = float(np.quantile(np.abs(np.concatenate([np.concatenate((item[4], item[5])) for item in prepared])), .98))
    figure, axes = plt.subplots(2, 4, figsize=(8.3, 4.35))
    for row, item in enumerate(prepared):
        selection, image, coordinates, dc_core, local_mreg, local_tfh, mask = item
        axes[row, 0].imshow(image, origin="upper")
        axes[row, 0].set_title(f"{selection['sample']} ({selection['stratum']})\nr={selection['mregdc_tfh_local_r']:.3f}", fontsize=8.4, weight="bold", pad=4)
        axes[row, 0].set_axis_off()
        dc_norm = plt.Normalize(vmin=np.quantile(dc_core, .02), vmax=np.quantile(dc_core, .98))
        for axis, values, title, cmap, norm in (
            (axes[row, 1], dc_core, "DC-core score", "YlGnBu", dc_norm),
            (axes[row, 2], local_mreg, "Local mregDC-like\nresidual field", "RdBu_r", TwoSlopeNorm(vcenter=0, vmin=-limit, vmax=limit)),
            (axes[row, 3], local_tfh, "Local Tfh-like\nresidual field", "RdBu_r", TwoSlopeNorm(vcenter=0, vmin=-limit, vmax=limit)),
        ):
            scatter = draw_panel(axis, image, coordinates, values, title, cmap, norm, mask if title != "DC-core score" else None)
            colorbar = figure.colorbar(scatter, ax=axis, fraction=.045, pad=.015)
            colorbar.ax.tick_params(labelsize=6)
    figure.suptitle("GSE277116 representative packages selected by within-stratum median effect", y=1.01, fontsize=10.5, weight="bold")
    figure.tight_layout()
    figure.savefig(OUT / "GSE277116_representative_spatial_maps.png", dpi=400, bbox_inches="tight")
    figure.savefig(OUT / "GSE277116_representative_spatial_maps.pdf", bbox_inches="tight")
    plt.close(figure)
    write_tsv(OUT / "GSE277116_representative_map_selection.tsv", selections)
    print(OUT)


if __name__ == "__main__":
    main()
