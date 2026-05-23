#!/usr/bin/env python3
"""Extract anchor sequences from a features.tsv file and save them as FASTQ.GZ.

Input:
  - features.tsv with columns: anchor_target, anchor, target

Output:
    - anchors.fastq.gz (all anchors in one FASTQ)
    - anchor_targets_fastq/anchor_XXXXXX.targets.fastq.gz (one FASTQ per anchor)
    - anchor_target_fastq_mapping.tsv

Anchors are deduplicated by default. For per-anchor target FASTQ output,
targets are deduplicated by default unless --keep-duplicate-targets is set.
FASTQ quality is filled with a single repeated character.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract anchors from features.tsv and write FASTQ.GZ")
    parser.add_argument("--features", type=Path, help="Input features.tsv")
    parser.add_argument("--out-dir", type=Path, help="Output directory")
    parser.add_argument(
        "--qual-char",
        default="I",
        help="Single quality character used to build pseudo FASTQ qualities",
    )
    parser.add_argument(
        "--keep-duplicate-targets",
        action="store_true",
        help="Keep duplicate targets within each anchor instead of deduplicating",
    )
    return parser.parse_args()


def load_anchor_targets(
    features_path: Path,
    keep_duplicate_targets: bool,
) -> tuple[list[str], dict[str, list[str]]]:
    anchors: list[str] = []
    anchor_seen: set[str] = set()
    anchor_to_targets: dict[str, list[str]] = {}
    target_seen_by_anchor: dict[str, set[str]] = {}

    with features_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "anchor" not in reader.fieldnames or "target" not in reader.fieldnames:
            raise ValueError(f"[ERROR] Missing anchor column in features file: {features_path}")
        for row in reader:
            anchor = (row.get("anchor") or "").strip().upper()
            if not anchor:
                continue
            if anchor not in anchor_seen:
                anchor_seen.add(anchor)
                anchors.append(anchor)
                anchor_to_targets[anchor] = []
                target_seen_by_anchor[anchor] = set()

            target = (row.get("target") or "").strip().upper()
            if not target:
                continue

            if keep_duplicate_targets:
                anchor_to_targets[anchor].append(target)
            else:
                if target not in target_seen_by_anchor[anchor]:
                    target_seen_by_anchor[anchor].add(target)
                    anchor_to_targets[anchor].append(target)

    if not anchors:
        raise ValueError(f"[ERROR] No anchor sequences loaded from: {features_path}")
    return anchors, anchor_to_targets


def write_fastq_gz(path: Path, sequences: list[str], qual_char: str, record_prefix: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for i, seq in enumerate(sequences, start=1):
            header = f"@{record_prefix}_{i:06d}|len={len(seq)}"
            qual = qual_char * len(seq)
            handle.write(f"{header}\n")
            handle.write(f"{seq}\n")
            handle.write("+\n")
            handle.write(f"{qual}\n")


def main() -> None:
    args = parse_args()

    if len(args.qual_char) != 1:
        raise ValueError("[ERROR] --qual-char must be a single character")
    if not args.features.is_file():
        raise FileNotFoundError(f"[ERROR] Features file not found: {args.features}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    anchors_fastq = args.out_dir / "anchors.fastq.gz"
    targets_dir = args.out_dir / "anchor_targets_fastq"
    mapping_tsv = args.out_dir / "anchor_target_fastq_mapping.tsv"
    targets_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Loading anchors from {args.features}")
    anchors, anchor_to_targets = load_anchor_targets(args.features, args.keep_duplicate_targets)
    print(f"[INFO] Loaded {len(anchors)} anchor sequences")

    print(f"[INFO] Writing all anchors FASTQ.GZ to {anchors_fastq}")
    write_fastq_gz(anchors_fastq, anchors, args.qual_char, "anchor")

    print(f"[INFO] Writing per-anchor targets FASTQ.GZ files to {targets_dir}")
    mapping_rows: list[dict[str, str | int]] = []

    for idx, anchor in enumerate(anchors, start=1):
        targets = anchor_to_targets.get(anchor, [])
        out_name = f"anchor_{idx:06d}.targets.fastq.gz"
        out_path = targets_dir / out_name
        write_fastq_gz(out_path, targets, args.qual_char, f"anchor_{idx:06d}_target")
        mapping_rows.append(
            {
                "anchor_index": idx,
                "anchor": anchor,
                "n_targets": len(targets),
                "targets_fastq": out_name,
            }
        )

    with mapping_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["anchor_index", "anchor", "n_targets", "targets_fastq"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    print(f"[OK] anchors_fastq={anchors_fastq}")
    print(f"[OK] per_anchor_targets_dir={targets_dir}")
    print(f"[OK] mapping_tsv={mapping_tsv}")


if __name__ == "__main__":
    main()
