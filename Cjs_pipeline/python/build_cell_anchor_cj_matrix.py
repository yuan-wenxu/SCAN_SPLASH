#!/usr/bin/env python3
"""Build a cell(barcode) x anchor sparse matrix from SPLASH cjs.tsv.gz.

Input cjs file must contain four columns:
    anchor  sample  barcode  Cj

During matrix construction, barcode values are reverse-mapped from
`mapped_sequence` to `combo_sequence` using barcode_mapping.csv.

Outputs (under --out-dir):
1) cell_anchor_cj.matrix.mtx.gz      (rows=cells/barcodes, cols=anchors)
2) cell_anchor_cj.cells.tsv          (cell metadata)
3) cell_anchor_cj.features.tsv       (anchor metadata)
4) build_cell_anchor_cj_matrix.summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build sparse cell(barcode) x anchor matrix from cjs.tsv.gz"
    )
    parser.add_argument(
        "--cjs-tsv-gz",
        type=Path,
        required=True,
        help="Path to SPLASH result_Cjs/cjs.tsv.gz",
    )
    parser.add_argument(
        "--barcode-mapping-csv",
        type=Path,
        required=True,
        help="Path to barcode_mapping.csv (must contain combo_sequence,mapped_sequence)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--drop-unmapped",
        action="store_true",
        help="Drop barcodes not found in mapping table (default: keep unmapped barcode as-is).",
    )
    return parser.parse_args()


def load_barcode_mapping(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            combo = row.get("combo_sequence", "").strip().upper()
            mapped = row.get("mapped_sequence", "").strip().upper()
            if combo and mapped:
                mapping[mapped] = combo
    if not mapping:
        raise ValueError(f"[ERROR] No valid mapping loaded from: {path}")
    return mapping


def write_cells(path: Path, cells: list[str], cell_to_sample: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["cell_id", "sample_id", "cell_barcode"])
        for i, cell in enumerate(cells):
            sample = cell_to_sample.get(cell, "")
            writer.writerow([f"{i}_{cell}", sample, cell])


def write_features(path: Path, anchors: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature_id", "anchor"])
        for i, anchor in enumerate(anchors):
            writer.writerow([f"anchor_{i:07d}", anchor])


def write_matrix_market_gz(
    path: Path,
    entries: dict[tuple[str, str], float],
    cells: list[str],
    anchors: list[str],
) -> None:
    cell_idx = {c: i + 1 for i, c in enumerate(cells)}
    anchor_idx = {a: i + 1 for i, a in enumerate(anchors)}

    triplets: list[tuple[int, int, float]] = []
    for (cell, anchor), v in entries.items():
        if v != 0.0:
            triplets.append((cell_idx[cell], anchor_idx[anchor], v))
    triplets.sort()

    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write("%%MatrixMarket matrix coordinate real general\n")
        fh.write("%\n")
        fh.write(f"{len(cells)} {len(anchors)} {len(triplets)}\n")
        for r, c, v in triplets:
            fh.write(f"{r} {c} {v:.10g}\n")


def main() -> None:
    args = parse_args()

    if not args.cjs_tsv_gz.is_file():
        raise FileNotFoundError(f"[ERROR] cjs file not found: {args.cjs_tsv_gz}")
    if not args.barcode_mapping_csv.is_file():
        raise FileNotFoundError(f"[ERROR] barcode mapping CSV not found: {args.barcode_mapping_csv}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    out_matrix = args.out_dir / "cell_anchor_cj.matrix.mtx.gz"
    out_cells = args.out_dir / "cell_anchor_cj.cells.tsv"
    out_features = args.out_dir / "cell_anchor_cj.features.tsv"
    out_summary = args.out_dir / "build_cell_anchor_cj_matrix.summary.tsv"

    print(f"[INFO] Loading barcode mapping: {args.barcode_mapping_csv}")
    bc_map = load_barcode_mapping(args.barcode_mapping_csv)
    print(f"[INFO] Loaded {len(bc_map)} mapped_sequence -> combo_sequence entries")

    entries: dict[tuple[str, str], float] = defaultdict(float)
    cells_set: set[str] = set()
    anchors_set: set[str] = set()
    cell_to_sample: dict[str, str] = {}

    total_rows = 0
    mapped_rows = 0
    dropped_unmapped = 0
    bad_rows = 0

    with gzip.open(args.cjs_tsv_gz, "rt", encoding="utf-8") as fh:
        header = fh.readline().strip().split()
        expected = ["anchor", "sample", "barcode", "Cj"]
        if header != expected:
            print(f"[WARNING] Unexpected header: {header} (expected {expected})")

        for line_no, line in enumerate(fh, start=2):
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 4:
                bad_rows += 1
                continue

            anchor, sample, mapped_barcode, cj_str = parts[0], parts[1], parts[2], parts[3]
            total_rows += 1

            try:
                cj = float(cj_str)
            except ValueError:
                bad_rows += 1
                continue

            combo_barcode = bc_map.get(mapped_barcode.upper())
            if combo_barcode is None:
                dropped_unmapped += 1
                if args.drop_unmapped:
                    continue
                combo_barcode = mapped_barcode.upper()

            mapped_rows += 1
            key = (combo_barcode, anchor)
            entries[key] += cj
            cells_set.add(combo_barcode)
            anchors_set.add(anchor)
            if combo_barcode not in cell_to_sample:
                cell_to_sample[combo_barcode] = sample
            elif cell_to_sample[combo_barcode] != sample:
                cell_to_sample[combo_barcode] = "multi"

            if total_rows % 2_000_000 == 0:
                print(f"[INFO] Processed rows: {total_rows}")

    cells = sorted(cells_set)
    anchors = sorted(anchors_set)

    print(f"[INFO] Writing cells: {out_cells}")
    write_cells(out_cells, cells, cell_to_sample)
    print(f"[INFO] Writing features: {out_features}")
    write_features(out_features, anchors)
    print(f"[INFO] Writing matrix: {out_matrix}")
    write_matrix_market_gz(out_matrix, entries, cells, anchors)

    with out_summary.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["total_rows", total_rows])
        writer.writerow(["mapped_rows", mapped_rows])
        writer.writerow(["dropped_unmapped", dropped_unmapped])
        writer.writerow(["bad_rows", bad_rows])
        writer.writerow(["n_cells", len(cells)])
        writer.writerow(["n_anchors", len(anchors)])
        writer.writerow(["nnz", sum(1 for v in entries.values() if v != 0.0)])

    print(f"[OK] total_rows={total_rows}")
    print(f"[OK] mapped_rows={mapped_rows}")
    print(f"[OK] dropped_unmapped={dropped_unmapped}")
    print(f"[OK] n_cells={len(cells)}")
    print(f"[OK] n_anchors={len(anchors)}")
    print(f"[OK] matrix={out_matrix}")


if __name__ == "__main__":
    main()
