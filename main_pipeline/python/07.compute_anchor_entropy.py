#!/usr/bin/env python3
"""
Compute per-anchor Shannon entropy from a remapped sparse matrix.

The remapped matrix directory should contain:
- remapped.features.tsv  (columns: anchor_target, anchor, target)
- remapped.cells.tsv
- remapped.matrix.mtx.gz  (cells x features, MatrixMarket format)

The anchor_align directory (sibling of the matrix directory, or specified via
--anchor-align-dir) should contain:
- anchors.fastq.gz    (anchor_id -> anchor sequence)
- assigned_reads.tsv  (SAM format; column 1 = anchor_id, XT:Z: tag = gene_id)

For each anchor, the script:
1. Sums counts across all cells for each target (column).
2. Computes Shannon entropy H = -sum(p_i * log2(p_i)) over the target count
   distribution.

Output CSV columns (sorted by entropy descending):
  anchor, gene_id, entropy, n_targets, total_count

Also writes a histogram PNG of entropy values.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute per-anchor Shannon entropy from a remapped sparse matrix"
    )
    parser.add_argument(
        "--matrix-dir",
        required=True,
        help="Directory containing remapped.features.tsv and remapped.matrix.mtx.gz",
    )
    parser.add_argument(
        "--anchor-align-dir",
        default=None,
        help=(
            "Directory containing anchors.fastq.gz and assigned_reads.tsv. "
            "Defaults to <matrix-dir>/../anchor_align"
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (default: <matrix-dir>/anchor_entropy.csv)",
    )
    parser.add_argument(
        "--histogram",
        default=None,
        help="Output histogram PNG path (default: <matrix-dir>/anchor_entropy_histogram.png)",
    )
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=1,
        help="Minimum total count for an anchor to be reported (default: 1)",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Feature loading
# ---------------------------------------------------------------------------

def load_features(features_path: Path) -> tuple[list[str], list[int]]:
    """
    Read remapped.features.tsv (anchor_target, anchor, target).
    Returns (anchor_list, feature_anchor_ids) where feature_anchor_ids[i]
    is the index into anchor_list for feature column i.
    """
    anchor_list: list[str] = []
    anchor_to_id: dict[str, int] = {}
    feature_anchor_ids: list[int] = []

    with open(features_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != ["anchor_target", "anchor", "target"]:
            print(f"[WARNING] unexpected features header: {header}", file=sys.stderr)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            anchor_seq = parts[1] if len(parts) > 1 else parts[0]
            if anchor_seq not in anchor_to_id:
                anchor_to_id[anchor_seq] = len(anchor_list)
                anchor_list.append(anchor_seq)
            feature_anchor_ids.append(anchor_to_id[anchor_seq])

    return anchor_list, feature_anchor_ids


# ---------------------------------------------------------------------------
# Anchor ID -> sequence mapping from anchors.fastq.gz
# ---------------------------------------------------------------------------

def load_anchor_id_to_sequence(fastq_path: Path) -> dict[str, str]:
    """
    Parse anchors.fastq.gz.  Each record header is:
      @anchor_XXXXXX|len=31|count=N
    Returns {anchor_id: sequence}, e.g. {"anchor_000001": "AAAA..."}.
    """
    mapping: dict[str, str] = {}
    open_fn = gzip.open if str(fastq_path).endswith(".gz") else open
    with open_fn(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # +
            fh.readline()  # quality
            anchor_id = header.lstrip("@").split("|")[0]
            mapping[anchor_id] = seq
    return mapping


# ---------------------------------------------------------------------------
# Anchor ID -> gene_id from assigned_reads.tsv (SAM format)
# ---------------------------------------------------------------------------

def load_anchor_gene_mapping(sam_path: Path) -> dict[str, str]:
    """
    Parse assigned_reads.tsv (SAM format).
    Column 0: query name = anchor_XXXXXX|len=31|count=N
    XT:Z:<gene_id> tag carries the gene assignment.
    Returns {anchor_id: gene_id}.  Anchors with no XT:Z: tag get "".
    """
    mapping: dict[str, str] = {}
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            anchor_id = fields[0].split("|")[0]
            gene_id = ""
            for field in fields[11:]:
                if field.startswith("XT:Z:"):
                    gene_id = field[5:]
                    break
            mapping[anchor_id] = gene_id
    return mapping


# ---------------------------------------------------------------------------
# Entropy
# ---------------------------------------------------------------------------

def compute_entropy(counts: list[int]) -> float:
    """Shannon entropy in bits (log base 2)."""
    total = sum(counts)
    if total == 0:
        return 0.0
    entropy = 0.0
    for c in counts:
        if c > 0:
            p = c / total
            entropy -= p * math.log2(p)
    return entropy


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    matrix_dir = Path(args.matrix_dir)
    features_path = matrix_dir / "remapped.features.tsv"
    matrix_path = matrix_dir / "remapped.matrix.mtx.gz"

    anchor_align_dir = (
        Path(args.anchor_align_dir)
        if args.anchor_align_dir
        else matrix_dir.parent / "anchor_align"
    )
    fastq_path = anchor_align_dir / "anchors.fastq.gz"
    sam_path = anchor_align_dir / "assigned_reads.tsv"

    output_path = Path(args.output) if args.output else matrix_dir / "anchor_entropy.csv"
    histogram_path = (
        Path(args.histogram) if args.histogram
        else matrix_dir / "anchor_entropy_histogram.png"
    )

    # --- Load anchor id <-> sequence <-> gene mappings ---
    print("[INFO] Loading anchor ID -> sequence mapping...", file=sys.stderr)
    anchor_id_to_seq = load_anchor_id_to_sequence(fastq_path)
    print(f"[INFO] Loaded {len(anchor_id_to_seq)} anchors from {fastq_path}", file=sys.stderr)

    print("[INFO] Loading anchor ID -> gene mapping...", file=sys.stderr)
    anchor_id_to_gene = load_anchor_gene_mapping(sam_path)
    print(f"[INFO] Loaded {len(anchor_id_to_gene)} entries from {sam_path}", file=sys.stderr)

    # Build anchor_sequence -> gene_id lookup
    seq_to_gene: dict[str, str] = {
        seq: anchor_id_to_gene.get(aid, "")
        for aid, seq in anchor_id_to_seq.items()
    }

    # --- Load features ---
    print("[INFO] Loading features...", file=sys.stderr)
    anchor_list, feature_anchor_ids = load_features(features_path)
    n_features = len(feature_anchor_ids)
    n_anchors = len(anchor_list)
    print(f"[INFO] Features: {n_features}, Anchors: {n_anchors}", file=sys.stderr)

    # anchor_feature_counts[anchor_id][feature_idx] -> total count across cells
    anchor_feature_counts: list[defaultdict[int, int]] = [
        defaultdict(int) for _ in range(n_anchors)
    ]

    # --- Parse MTX ---
    print("[INFO] Reading matrix...", file=sys.stderr)
    open_fn = gzip.open if str(matrix_path).endswith(".gz") else open

    feature_on_cols = False
    n_entries = 0

    with open_fn(matrix_path, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            n_rows, n_cols, n_entries = int(parts[0]), int(parts[1]), int(parts[2])
            if n_cols == n_features:
                feature_on_cols = True
            elif n_rows == n_features:
                feature_on_cols = False
            else:
                print(
                    f"[ERROR] Matrix dimensions ({n_rows}x{n_cols}) do not match",
                    file=sys.stderr,
                )
                sys.exit(1)
            print(
                f"[INFO] Matrix: {n_rows}x{n_cols}, nnz={n_entries}, "
                f"[INFO] orientation={'cells x features' if feature_on_cols else 'features x cells'}",
                file=sys.stderr,
            )
            break

        processed = 0
        report_every = 2_000_000
        for line in fh:
            parts = line.split()
            r, c, v = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2])
            feature_idx = c if feature_on_cols else r
            anchor_id = feature_anchor_ids[feature_idx]
            anchor_feature_counts[anchor_id][feature_idx] += v
            processed += 1
            if processed % report_every == 0:
                print(f"[INFO] Processed {processed}/{n_entries} entries...", file=sys.stderr)

    print(f"[INFO] Done. Processed {processed} entries.", file=sys.stderr)

    # --- Compute entropy per anchor ---
    print("[INFO] Computing entropy...", file=sys.stderr)
    rows: list[dict] = []
    for anchor_idx, anchor_seq in enumerate(anchor_list):
        counts = list(anchor_feature_counts[anchor_idx].values())
        total = sum(counts)
        if total < args.min_total_count:
            continue
        n_targets = len(counts)
        entropy = compute_entropy(counts)
        gene_id = seq_to_gene.get(anchor_seq, "")
        rows.append(
            {
                "anchor": anchor_seq,
                "gene_id": gene_id,
                "entropy": entropy,
                "n_targets": n_targets,
                "total_count": total,
            }
        )

    # Sort by entropy descending
    rows.sort(key=lambda x: x["entropy"], reverse=True)

    # --- Write CSV ---
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["anchor", "gene_id", "entropy", "n_targets", "total_count"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "anchor": row["anchor"],
                    "gene_id": row["gene_id"],
                    "entropy": f"{row['entropy']:.6f}",
                    "n_targets": row["n_targets"],
                    "total_count": row["total_count"],
                }
            )
    print(f"[INFO] CSV written to: {output_path}", file=sys.stderr)

    # --- Histogram ---
    print("[INFO] Generating histogram...", file=sys.stderr)
    entropies = [r["entropy"] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(entropies, bins=50, color="steelblue", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Shannon entropy (bits)", fontsize=12)
    ax.set_ylabel("Number of anchors", fontsize=12)
    ax.set_title(f"[INFO] Per-anchor entropy distribution (n={len(entropies)})", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(histogram_path, dpi=150)
    plt.close(fig)
    print(f"[INFO] Histogram written to: {histogram_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
