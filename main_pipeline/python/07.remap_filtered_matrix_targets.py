#!/usr/bin/env python3
"""
Remap targets in a filtered sparse matrix using per-anchor cluster mapping TSV files.

Input matrix directory should contain:
- whitelist_filtered.features.tsv
- whitelist_filtered.cells.tsv
- whitelist_filtered.matrix.mtx.gz
- anchor_align/anchors.fastq.gz

Mapping directory should contain files like:
- anchor_000027.cluster_mapping.tsv
with columns:
- target_sequence
- reference_sequence

This script will:
1. Convert each feature target to mapped reference target (within the same anchor).
2. Merge features that become identical after remapping (same anchor + mapped target).
3. Rebuild MatrixMarket matrix with merged feature rows.
4. Write a new matrix directory with remapped features, copied cells, and remapped matrix.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remap filtered matrix targets to reference targets and merge duplicated features"
    )
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        required=True,
        help="Directory containing whitelist_filtered.{features,cells,matrix} and anchor_align/anchors.fastq.gz",
    )
    parser.add_argument(
        "--mapping-dir",
        type=Path,
        required=True,
        help="Directory containing anchor_*.cluster_mapping.tsv files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for remapped matrix files",
    )
    parser.add_argument(
        "--features-file",
        type=str,
        default="whitelist_filtered.features.tsv",
        help="Features TSV filename inside --matrix-dir (default: whitelist_filtered.features.tsv)",
    )
    parser.add_argument(
        "--cells-file",
        type=str,
        default="whitelist_filtered.cells.tsv",
        help="Cells TSV filename inside --matrix-dir (default: whitelist_filtered.cells.tsv)",
    )
    parser.add_argument(
        "--matrix-file",
        type=str,
        default="whitelist_filtered.matrix.mtx.gz",
        help="MatrixMarket filename inside --matrix-dir (default: whitelist_filtered.matrix.mtx.gz)",
    )
    parser.add_argument(
        "--anchors-fastq",
        type=str,
        default="anchor_align/anchors.fastq.gz",
        help="Anchor FASTQ path relative to --matrix-dir (default: anchor_align/anchors.fastq.gz)",
    )
    return parser.parse_args()


def load_anchor_id_to_sequence(anchor_fastq_gz: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with gzip.open(anchor_fastq_gz, "rt", encoding="utf-8") as fh:
        rec = 0
        header = ""
        for line in fh:
            rec_mod = rec % 4
            if rec_mod == 0:
                header = line.strip()
            elif rec_mod == 1:
                seq = line.strip()
                if header.startswith("@"):
                    anchor_id = header[1:].split("|", 1)[0]
                    mapping[anchor_id] = seq
            rec += 1
    return mapping


def parse_mapping_line(line: str) -> Tuple[str, str] | None:
    s = line.strip()
    if not s:
        return None
    parts_tab = s.split("\t")
    if len(parts_tab) >= 2:
        return parts_tab[0], parts_tab[1]
    parts_ws = s.split()
    if len(parts_ws) >= 2:
        return parts_ws[0], parts_ws[1]
    return None


def load_target_mapping(
    mapping_dir: Path, anchor_id_to_seq: Dict[str, str]
) -> Tuple[Dict[Tuple[str, str], str], set]:
    """Returns (target_map, mapped_anchors) where mapped_anchors is the set of
    anchor sequences that have a cluster_mapping.tsv file."""
    target_map: Dict[Tuple[str, str], str] = {}
    mapped_anchors: set = set()

    for tsv_file in sorted(mapping_dir.glob("*.cluster_mapping.tsv")):
        anchor_id = tsv_file.name.replace(".cluster_mapping.tsv", "")
        anchor_seq = anchor_id_to_seq.get(anchor_id)
        if anchor_seq is None:
            continue

        mapped_anchors.add(anchor_seq)

        with open(tsv_file, "r", encoding="utf-8") as fh:
            first = fh.readline()
            first_parsed = parse_mapping_line(first)
            if first_parsed and first_parsed[0] == "target_sequence":
                pass
            else:
                if first_parsed:
                    target, ref = first_parsed
                    target_map[(anchor_seq, target)] = ref

            for line in fh:
                parsed = parse_mapping_line(line)
                if not parsed:
                    continue
                target, ref = parsed
                target_map[(anchor_seq, target)] = ref

    return target_map, mapped_anchors


def remap_features(
    features_path: Path,
    target_map: Dict[Tuple[str, str], str],
    mapped_anchors: set,
) -> Tuple[List[Tuple[str, str, str]], List[int], int]:
    """
    Returns:
    - new_features rows: (anchor_target, anchor, target)
    - old_to_new_row: 1-based old row index -> 1-based new row index (0 = drop)
    - old_row_count

    Features whose anchor has no cluster_mapping.tsv are dropped (mapped to 0).
    """
    new_features: List[Tuple[str, str, str]] = []
    old_to_new_row: List[int] = [0]
    key_to_new_row: Dict[Tuple[str, str], int] = {}

    with open(features_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row_idx, row in enumerate(reader, start=1):
            anchor = row["anchor"]
            target = row["target"]

            if anchor not in mapped_anchors:
                old_to_new_row.append(0)
                continue

            mapped_target = target_map.get((anchor, target), target)

            key = (anchor, mapped_target)
            new_idx = key_to_new_row.get(key)
            if new_idx is None:
                new_idx = len(new_features) + 1
                key_to_new_row[key] = new_idx
                new_features.append((f"{anchor}|{mapped_target}", anchor, mapped_target))

            old_to_new_row.append(new_idx)

    old_row_count = len(old_to_new_row) - 1
    return new_features, old_to_new_row, old_row_count


def transform_matrix_with_external_sort(
    matrix_path: Path,
    old_to_new_feature: List[int],
    old_feature_count: int,
    new_feature_count: int,
    n_cells_expected: int,
    out_matrix_path: Path,
) -> Tuple[int, int, int, str]:
    """Remap feature indices (on row or column axis) and aggregate duplicates.

    Returns:
    - out_n_rows
    - out_n_cols
    - new_nnz
    - orientation: "feature_on_rows" or "feature_on_cols"
    """
    with tempfile.TemporaryDirectory(prefix="remap_matrix_tmp_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        unsorted_path = tmpdir_path / "triplets.unsorted.tsv"
        sorted_path = tmpdir_path / "triplets.sorted.tsv"
        agg_path = tmpdir_path / "triplets.agg.tsv"

        old_n_rows = -1
        old_n_cols = -1
        orientation = ""

        with gzip.open(matrix_path, "rt", encoding="utf-8") as fin, open(
            unsorted_path, "w", encoding="utf-8"
        ) as fout:
            header_seen = False
            for raw in fin:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("%"):  # MatrixMarket header/comments
                    continue

                if not header_seen:
                    parts = line.split()
                    if len(parts) != 3:
                        raise ValueError(f"[ERROR] Invalid MatrixMarket size line: {line}")
                    old_n_rows = int(parts[0])
                    old_n_cols = int(parts[1])
                    header_seen = True

                    if old_n_rows == old_feature_count and old_n_cols == n_cells_expected:
                        orientation = "feature_on_rows"
                    elif old_n_cols == old_feature_count and old_n_rows == n_cells_expected:
                        orientation = "feature_on_cols"
                    else:
                        raise ValueError(
                            "[ERROR] Matrix dimensions do not match features/cells. "
                            f"matrix=({old_n_rows},{old_n_cols}), features={old_feature_count}, cells={n_cells_expected}"
                        )
                    continue

                parts = line.split()
                if len(parts) < 3:
                    continue

                row = int(parts[0])
                col = int(parts[1])
                val = int(parts[2])

                if orientation == "feature_on_rows":
                    if row >= len(old_to_new_feature):
                        raise IndexError(
                            f"[ERROR] Matrix row index {row} exceeds feature mapping size {len(old_to_new_feature) - 1}"
                        )
                    new_row = old_to_new_feature[row]
                    new_col = col
                else:
                    if col >= len(old_to_new_feature):
                        raise IndexError(
                            f"[ERROR] Matrix col index {col} exceeds feature mapping size {len(old_to_new_feature) - 1}"
                        )
                    new_row = row
                    new_col = old_to_new_feature[col]

                # 0 means the feature was dropped (no mapping file for that anchor)
                if (orientation == "feature_on_rows" and new_row == 0) or (
                    orientation == "feature_on_cols" and new_col == 0
                ):
                    continue

                fout.write(f"{new_row}\t{new_col}\t{val}\n")

        if old_n_rows <= 0 or old_n_cols <= 0:
            raise ValueError("[ERROR] Failed to read MatrixMarket dimensions")

        sort_cmd = [
            "sort",
            "-T",
            str(tmpdir_path),
            "-k1,1n",
            "-k2,2n",
            str(unsorted_path),
        ]
        with open(sorted_path, "w", encoding="utf-8") as sorted_out:
            subprocess.run(sort_cmd, check=True, stdout=sorted_out, env={**os.environ, "LC_ALL": "C"})

        new_nnz = 0
        prev_row = None
        prev_col = None
        acc_val = 0

        with open(sorted_path, "r", encoding="utf-8") as fin, open(
            agg_path, "w", encoding="utf-8"
        ) as fout:
            for line in fin:
                row_s, col_s, val_s = line.rstrip("\n").split("\t")
                row = int(row_s)
                col = int(col_s)
                val = int(val_s)

                if prev_row is None:
                    prev_row, prev_col, acc_val = row, col, val
                    continue

                if row == prev_row and col == prev_col:
                    acc_val += val
                else:
                    if acc_val != 0:
                        fout.write(f"{prev_row} {prev_col} {acc_val}\n")
                        new_nnz += 1
                    prev_row, prev_col, acc_val = row, col, val

            if prev_row is not None and acc_val != 0:
                fout.write(f"{prev_row} {prev_col} {acc_val}\n")
                new_nnz += 1

        if orientation == "feature_on_rows":
            out_n_rows = new_feature_count
            out_n_cols = old_n_cols
        else:
            out_n_rows = old_n_rows
            out_n_cols = new_feature_count

        with gzip.open(out_matrix_path, "wt", encoding="utf-8") as out:
            out.write("%%MatrixMarket matrix coordinate integer general\n")
            out.write("%\n")
            out.write(f"{out_n_rows} {out_n_cols} {new_nnz}\n")
            with open(agg_path, "r", encoding="utf-8") as agg_in:
                shutil.copyfileobj(agg_in, out)

    return out_n_rows, out_n_cols, new_nnz, orientation


def write_features(features_out: Path, features_rows: List[Tuple[str, str, str]]) -> None:
    with open(features_out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["anchor_target", "anchor", "target"])
        writer.writerows(features_rows)


def count_cells(cells_path: Path) -> int:
    count = 0
    with open(cells_path, "r", encoding="utf-8") as fh:
        for _ in fh:
            count += 1
    return max(count - 1, 0)


def main() -> None:
    args = parse_args()

    matrix_dir = args.matrix_dir
    mapping_dir = args.mapping_dir
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    features_path = matrix_dir / args.features_file
    cells_path = matrix_dir / args.cells_file
    matrix_path = matrix_dir / args.matrix_file
    anchors_fastq_path = matrix_dir / args.anchors_fastq

    if not features_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing features file: {features_path}")
    if not cells_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing cells file: {cells_path}")
    if not matrix_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing matrix file: {matrix_path}")
    if not anchors_fastq_path.exists():
        raise FileNotFoundError(f"[ERROR] Missing anchors FASTQ: {anchors_fastq_path}")
    if not mapping_dir.exists():
        raise FileNotFoundError(f"[ERROR] Missing mapping directory: {mapping_dir}")

    print("[INFO] Loading anchor ID to sequence mapping...")
    anchor_id_to_seq = load_anchor_id_to_sequence(anchors_fastq_path)
    print(f"[INFO] Loaded {len(anchor_id_to_seq)} anchor IDs from {anchors_fastq_path.name}")

    print("[INFO] Loading target remapping table from cluster_mapping TSV files...")
    target_map, mapped_anchors = load_target_mapping(mapping_dir, anchor_id_to_seq)
    print(
        f"[INFO] Loaded {len(target_map)} (anchor, target) -> ref_target mappings "
        f"across {len(mapped_anchors)} anchors"
    )

    print("[INFO] Remapping features and building row index mapping...")
    new_features, old_to_new_row, old_feature_count = remap_features(
        features_path, target_map, mapped_anchors
    )
    print(
        f"[INFO] Features: old={old_feature_count}, new={len(new_features)}, "
        f"[INFO] Merged={old_feature_count - len(new_features)}"
    )

    out_features = out_dir / f"remapped.features.tsv"
    out_cells = out_dir / f"remapped.cells.tsv"
    out_matrix = out_dir / f"remapped.matrix.mtx.gz"

    write_features(out_features, new_features)
    shutil.copy2(cells_path, out_cells)

    n_cells = count_cells(cells_path)
    print("[INFO] Transforming and aggregating matrix entries (external sort)...")
    out_n_rows, out_n_cols, new_nnz, orientation = transform_matrix_with_external_sort(
        matrix_path=matrix_path,
        old_to_new_feature=old_to_new_row,
        old_feature_count=old_feature_count,
        new_feature_count=len(new_features),
        n_cells_expected=n_cells,
        out_matrix_path=out_matrix,
    )

    print(f"[INFO] Matrix orientation: {orientation}")
    print(f"[INFO] Matrix remapped: rows={out_n_rows}, cols={out_n_cols}, nnz={new_nnz}")
    print(f"[INFO] Output features: {out_features}")
    print(f"[INFO] Output cells:    {out_cells}")
    print(f"[INFO] Output matrix:   {out_matrix}")


if __name__ == "__main__":
    main()
