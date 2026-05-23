#!/usr/bin/env python3
"""Filter matrix features by anchor detected-cell percentage.

Input matrix is expected to be cell x anchor-target features:
- whitelist_filtered.matrix.mtx.gz
- whitelist_filtered.cells.tsv
- whitelist_filtered.features.tsv

This script groups features by anchor, computes in how many cells each anchor is
present (any non-zero value across its features), and keeps only anchors that
appear in at least min-anchor-cell-pct of cells.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

import numpy as np
from scipy.io import mmread, mmwrite
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter anchor-target matrix by anchor detected-cell percentage"
    )
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        required=True,
        help="Directory containing whitelist_filtered.{matrix,cells,features}",
    )
    parser.add_argument(
        "--input-prefix",
        default="whitelist_filtered",
        help="Input filename prefix inside --matrix-dir (default: whitelist_filtered)",
    )
    parser.add_argument(
        "--output-prefix",
        default="anchor_pct_filtered",
        help="Output filename prefix inside --matrix-dir (default: anchor_pct_filtered)",
    )
    parser.add_argument(
        "--min-anchor-cell-pct",
        type=float,
        default=10.0,
        help="Minimum percentage of cells where anchor must be detected (default: 10.0)",
    )
    return parser.parse_args()


def load_matrix(path: Path):
    with gzip.open(path, "rb") as fh:
        matrix = mmread(fh)
    return matrix.tocsr()


def load_cells(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def load_features(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def write_cells(path: Path, cells: list[dict[str, str]]) -> None:
    fieldnames = ["cell_id", "sample_id", "cell_barcode"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for cell in cells:
            writer.writerow({k: cell.get(k, "") for k in fieldnames})


def write_features(path: Path, features: list[dict[str, str]], keep_mask: np.ndarray) -> None:
    fieldnames = ["anchor_target", "anchor", "target"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for keep, feature in zip(keep_mask, features):
            if keep:
                writer.writerow({k: feature.get(k, "") for k in fieldnames})


def write_matrix(path: Path, matrix) -> None:
    with gzip.open(path, "wb") as fh:
        mmwrite(fh, matrix)


def build_distribution_plot(
    out_png: Path,
    feature_detected_cells: np.ndarray,
    anchor_detected_cells: np.ndarray,
    n_cells: int,
    min_cells: int,
    min_anchor_cell_fraction: float,
    title_suffix: str,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ax_feature_cells = axes[0, 0]
    ax_anchor_cells = axes[0, 1]
    ax_feature_frac = axes[1, 0]
    ax_anchor_frac = axes[1, 1]

    feature_cells = np.asarray(feature_detected_cells, dtype=float)
    feature_cells = feature_cells[np.isfinite(feature_cells)]
    ax_feature_cells.hist(feature_cells, bins=80, color="#1F77B4", alpha=0.85, edgecolor="white", linewidth=0.5)
    ax_feature_cells.set_xlabel("Detected cells per feature", fontsize=12)
    ax_feature_cells.set_ylabel("Feature count", fontsize=12)
    ax_feature_cells.set_title(f"Feature detected-cell distribution ({title_suffix})", fontsize=13, fontweight="bold")

    anchor_cells = np.asarray(anchor_detected_cells, dtype=float)
    anchor_cells = anchor_cells[np.isfinite(anchor_cells)]
    ax_anchor_cells.hist(anchor_cells, bins=80, color="#2CA02C", alpha=0.85, edgecolor="white", linewidth=0.5)
    ax_anchor_cells.axvline(min_cells, color="black", linestyle="--", linewidth=1.3, label=f"threshold={min_cells}")
    ax_anchor_cells.set_xlabel("Detected cells per anchor", fontsize=12)
    ax_anchor_cells.set_ylabel("Anchor count", fontsize=12)
    ax_anchor_cells.set_title(f"Anchor detected-cell distribution ({title_suffix})", fontsize=13, fontweight="bold")
    ax_anchor_cells.legend(frameon=False)

    feature_frac = feature_cells / n_cells if feature_cells.size > 0 and n_cells > 0 else np.array([])
    ax_feature_frac.hist(feature_frac, bins=80, color="#1F77B4", alpha=0.85, edgecolor="white", linewidth=0.5)
    ax_feature_frac.set_xlabel("Feature detected-cell fraction", fontsize=12)
    ax_feature_frac.set_ylabel("Feature count", fontsize=12)
    ax_feature_frac.set_title(f"Feature fraction distribution ({title_suffix})", fontsize=13, fontweight="bold")

    anchor_frac = anchor_cells / n_cells if anchor_cells.size > 0 and n_cells > 0 else np.array([])
    ax_anchor_frac.hist(anchor_frac, bins=80, color="#2CA02C", alpha=0.85, edgecolor="white", linewidth=0.5)
    ax_anchor_frac.axvline(
        min_anchor_cell_fraction,
        color="black",
        linestyle="--",
        linewidth=1.3,
        label=f"threshold={min_anchor_cell_fraction:.4f}",
    )
    ax_anchor_frac.set_xlabel("Anchor detected-cell fraction", fontsize=12)
    ax_anchor_frac.set_ylabel("Anchor count", fontsize=12)
    ax_anchor_frac.set_title(f"Anchor fraction distribution ({title_suffix})", fontsize=13, fontweight="bold")
    ax_anchor_frac.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    matrix_dir = args.matrix_dir.resolve()
    in_prefix = args.input_prefix
    out_prefix = args.output_prefix

    in_matrix = matrix_dir / f"{in_prefix}.matrix.mtx.gz"
    in_cells = matrix_dir / f"{in_prefix}.cells.tsv"
    in_features = matrix_dir / f"{in_prefix}.features.tsv"

    out_matrix = matrix_dir / f"{out_prefix}.matrix.mtx.gz"
    out_cells = matrix_dir / f"{out_prefix}.cells.tsv"
    out_features = matrix_dir / f"{out_prefix}.features.tsv"
    out_anchor_stats = matrix_dir / f"{out_prefix}.anchor_detected_cells.csv"
    out_dist_before = matrix_dir / f"{out_prefix}.detected_cells_fraction.before_filter.png"
    out_dist_after = matrix_dir / f"{out_prefix}.detected_cells_fraction.after_filter.png"

    for p in [in_matrix, in_cells, in_features]:
        if not p.is_file():
            raise FileNotFoundError(f"[ERROR] Missing required file: {p}")

    print(f"[INFO] Loading matrix from: {in_matrix}")
    matrix = load_matrix(in_matrix)
    print(f"[INFO] Loading cells from: {in_cells}")
    cells = load_cells(in_cells)
    print(f"[INFO] Loading features from: {in_features}")
    features = load_features(in_features)

    if matrix.shape[0] != len(cells):
        raise ValueError(
            f"[ERROR] Matrix rows {matrix.shape[0]} != cells {len(cells)}"
        )
    if matrix.shape[1] != len(features):
        raise ValueError(
            f"[ERROR] Matrix cols {matrix.shape[1]} != features {len(features)}"
        )

    n_cells = matrix.shape[0]
    min_cells = max(1, math.ceil(n_cells * args.min_anchor_cell_pct / 100.0))
    print(f"[INFO] Total cells: {n_cells}")
    print(f"[INFO] Anchor keep threshold: >= {min_cells} cells ({args.min_anchor_cell_pct:.2f}%)")

    anchor_to_indices: dict[str, list[int]] = {}
    for idx, feat in enumerate(features):
        anchor = feat.get("anchor", "").strip()
        if not anchor:
            continue
        if anchor not in anchor_to_indices:
            anchor_to_indices[anchor] = []
        anchor_to_indices[anchor].append(idx)

    anchor_records: list[dict[str, object]] = []
    kept_anchors: set[str] = set()

    for anchor in sorted(anchor_to_indices):
        col_idx = anchor_to_indices[anchor]
        sub = matrix[:, col_idx]
        detected_cells = int(np.count_nonzero(sub.getnnz(axis=1) > 0))
        keep = int(detected_cells >= min_cells)
        if keep:
            kept_anchors.add(anchor)
        anchor_records.append(
            {
                "anchor": anchor,
                "detected_cells": detected_cells,
                "detected_cell_fraction": detected_cells / n_cells,
                "n_anchor_targets": len(col_idx),
                "keep": keep,
            }
        )

    feature_detected_cells = matrix.getnnz(axis=0)
    anchor_detected_cells = np.array([int(x["detected_cells"]) for x in anchor_records], dtype=int)
    build_distribution_plot(
        out_png=out_dist_before,
        feature_detected_cells=feature_detected_cells,
        anchor_detected_cells=anchor_detected_cells,
        n_cells=n_cells,
        min_cells=min_cells,
        min_anchor_cell_fraction=args.min_anchor_cell_pct / 100.0,
        title_suffix="before filter",
    )

    feature_keep_mask = np.array(
        [feat.get("anchor", "").strip() in kept_anchors for feat in features],
        dtype=bool,
    )
    kept_matrix = matrix[:, feature_keep_mask]
    filtered_anchor_detected_cells = np.array(
        [int(x["detected_cells"]) for x in anchor_records if int(x["keep"]) == 1],
        dtype=int,
    )
    build_distribution_plot(
        out_png=out_dist_after,
        feature_detected_cells=kept_matrix.getnnz(axis=0),
        anchor_detected_cells=filtered_anchor_detected_cells,
        n_cells=n_cells,
        min_cells=min_cells,
        min_anchor_cell_fraction=args.min_anchor_cell_pct / 100.0,
        title_suffix="after filter",
    )
    print(f"[INFO] Kept {len(kept_anchors)} anchors out of {len(anchor_to_indices)} total anchors")
    print(f"[INFO] Kept {int(feature_keep_mask.sum())} features out of {len(features)} total features")

    print(f"[INFO] Writing filtered matrix to: {out_matrix}")
    write_matrix(out_matrix, kept_matrix)
    write_cells(out_cells, cells)
    write_features(out_features, features, feature_keep_mask)

    with out_anchor_stats.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "anchor",
                "detected_cells",
                "detected_cell_fraction",
                "n_anchor_targets",
                "keep",
            ],
        )
        writer.writeheader()
        writer.writerows(anchor_records)

    print(f"[INFO] Saved plot: {out_dist_before}")
    print(f"[INFO] Saved plot: {out_dist_after}")


if __name__ == "__main__":
    main()
