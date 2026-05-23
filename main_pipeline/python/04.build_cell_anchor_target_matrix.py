#!/usr/bin/env python3
"""Build a cell x anchor-target matrix from SPLASH significant SATC output.

The input is the binary SATC file produced by SPLASH with
--keep_significant_anchors_satc:

    <out-dir>/<prefix>_satc/after_correction.scores.satc

The script uses satc_dump to decode it, then writes:

1. <prefix>.cell_anchor_target_counts.matrix.mtx.gz
   Sparse Matrix Market matrix with rows=cells and columns=anchor-target pairs.

2. <prefix>.cell_anchor_target_counts.cells.tsv
   Row metadata.

3. <prefix>.cell_anchor_target_counts.features.tsv
   Column metadata.

By default the intermediate SATC text dump is deleted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build sparse cell x anchor-target counts from SPLASH SATC output"
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--satc-dump-bin", type=Path, required=True)
    parser.add_argument("--anchor-len", type=int, default=31)
    parser.add_argument("--target-len", type=int, default=31)
    parser.add_argument(
        "--keep-dump-text",
        action="store_true",
        help="Keep the intermediate SATC text dump",
    )
    parser.add_argument(
        "--cbc-mapping-csv",
        type=Path,
        help="Optional CSV mapping file to reverse-map 16bp CBC back to original combo sequence",
    )
    return parser.parse_args()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def run_satc_dump(satc_dump_bin: Path, satc_path: Path, dump_path: Path) -> None:
    cmd = [str(satc_dump_bin), "--format", "splash", str(satc_path), str(dump_path)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def load_cbc_mapping(mapping_csv: Path) -> Dict[str, str]:
    """Load reverse mapping from mapped_sequence (col 4) to combo_sequence (col 2)."""
    mapping: Dict[str, str] = {}
    with mapping_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"[ERROR] Empty or invalid CSV: {mapping_csv}")
        for row in reader:
            combo_seq = row.get("combo_sequence", "").strip().upper()
            mapped_seq = row.get("mapped_sequence", "").strip().upper()
            if combo_seq and mapped_seq:
                mapping[mapped_seq] = combo_seq
    if not mapping:
        raise ValueError(f"[ERROR] No valid mappings loaded from: {mapping_csv}")
    return mapping


def parse_cell_id(cell_id: str) -> Tuple[str, str]:
    sample_id, cell_barcode = cell_id.rsplit("_", 1)
    return sample_id, cell_barcode


def iter_dump_records(dump_path: Path, anchor_len: int, target_len: int):
    expected_len = anchor_len + target_len
    with dump_path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"[ERROR] Unexpected SATC dump line {line_no}: {line}")

            count_str, anchor_target, cell_id = parts
            if len(anchor_target) != expected_len:
                raise ValueError(
                    f"[ERROR] Unexpected anchor+target length {len(anchor_target)} "
                    f"on line {line_no}; expected {expected_len}"
                )

            anchor = anchor_target[:anchor_len]
            target = anchor_target[anchor_len:]
            yield cell_id, anchor, target, int(count_str)


def load_counts(
    dump_path: Path,
    anchor_len: int,
    target_len: int,
    cbc_mapping: Dict[str, str] | None = None,
):
    counts: Dict[Tuple[str, str], int] = defaultdict(int)
    cells = set()
    features = set()
    unmapped_cbcs = set()

    for cell_id, anchor, target, count in iter_dump_records(dump_path, anchor_len, target_len):
        # Optionally reverse-map CBC from 16bp back to original combo sequence
        if cbc_mapping is not None:
            sample_id, cell_barcode = parse_cell_id(cell_id)
            original_cbc = cbc_mapping.get(cell_barcode.upper())
            if original_cbc is None:
                unmapped_cbcs.add(cell_barcode)
                continue  # Skip unmapped CBCs
            cell_id = f"{sample_id}_{original_cbc}"

        feature_id = f"{anchor}|{target}"
        counts[(cell_id, feature_id)] += count
        cells.add(cell_id)
        features.add(feature_id)

    if unmapped_cbcs:
        print(f"[WARNING] {len(unmapped_cbcs)} CBCs not found in mapping, skipped.")

    return counts, sorted(cells), sorted(features)



def write_cells(path: Path, cells: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cell_id", "sample_id", "cell_barcode"])
        for cell_id in cells:
            sample_id, cell_barcode = parse_cell_id(cell_id)
            writer.writerow([cell_id, sample_id, cell_barcode])


def write_features(path: Path, features: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["anchor_target", "anchor", "target"])
        for feature_id in features:
            anchor, target = feature_id.split("|", 1)
            writer.writerow([feature_id, anchor, target])


def write_matrix_market(
    path: Path,
    counts: Dict[Tuple[str, str], int],
    cells: list[str],
    features: list[str],
) -> None:
    cell_index = {cell_id: i + 1 for i, cell_id in enumerate(cells)}
    feature_index = {feature_id: i + 1 for i, feature_id in enumerate(features)}

    entries = [
        (cell_index[cell_id], feature_index[feature_id], count)
        for (cell_id, feature_id), count in counts.items()
        if count
    ]
    entries.sort()

    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write(f"{len(cells)} {len(features)} {len(entries)}\n")
        for row, col, count in entries:
            handle.write(f"{row} {col} {count}\n")


def main() -> None:
    args = parse_args()

    os.makedirs(args.out_dir / "matrix", exist_ok=True)
    satc_path = args.out_dir / "result_satc" / "after_correction.scores.satc"
    dump_path = args.out_dir / "result_satc" / "after_correction.scores.satc.dump.tsv"
    matrix_path = args.out_dir / "matrix" / "cell_anchor_target_counts.matrix.mtx.gz"
    cells_path = args.out_dir / "matrix" / "cell_anchor_target_counts.cells.tsv"
    features_path = args.out_dir / "matrix" / "cell_anchor_target_counts.features.tsv"

    if not satc_path.is_file():
        raise FileNotFoundError(f"[ERROR] SATC file not found: {satc_path}")
    if not args.satc_dump_bin.is_file():
        raise FileNotFoundError(f"[ERROR] satc_dump binary not found: {args.satc_dump_bin}")

    run_satc_dump(args.satc_dump_bin, satc_path, dump_path)

    cbc_mapping = None
    if args.cbc_mapping_csv:
        if not args.cbc_mapping_csv.is_file():
            raise FileNotFoundError(f"[ERROR] CBC mapping CSV not found: {args.cbc_mapping_csv}")
        print(f"[INFO] Loading CBC mapping from {args.cbc_mapping_csv}")
        cbc_mapping = load_cbc_mapping(args.cbc_mapping_csv)
        print(f"[INFO] Loaded {len(cbc_mapping)} CBC mappings")

    counts, cells, features = load_counts(
        dump_path, args.anchor_len, args.target_len, cbc_mapping
    )

    #write_matrix(counts, cells, features)
    write_cells(cells_path, cells)
    write_features(features_path, features)
    write_matrix_market(matrix_path, counts, cells, features)

    if not args.keep_dump_text:
        dump_path.unlink(missing_ok=True)

    print(f"[INFO] cells={len(cells)}")
    print(f"[INFO] anchor_targets={len(features)}")
    print(f"[INFO] nonzero={sum(1 for count in counts.values() if count)}")
    print(f"[INFO] matrix_market={matrix_path}")
    print(f"[INFO] cells={cells_path}")
    print(f"[INFO] features={features_path}")


if __name__ == "__main__":
    main()
