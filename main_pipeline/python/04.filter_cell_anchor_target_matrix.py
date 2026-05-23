#!/usr/bin/env python3
"""Filter cell x anchor-target matrix by whitelist and cell-level QC.

The script reads the matrix, cells metadata, features metadata, and a whitelist
TSV where column 2 is the original combo sequence. Cells are classified into:

1. in_whitelist_pass: in whitelist and passes QC thresholds
2. in_whitelist_fail: in whitelist but fails QC thresholds
3. not_in_whitelist_fail: not in whitelist

It writes a filtered matrix containing only in_whitelist_pass cells and draws a
QC scatter plot using total counts and detected features per cell.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread, mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter cell x anchor-target matrix by whitelist and QC thresholds"
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--whitelist", type=Path, required=True, help="Path to whitelist TSV")
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=20000,
        help="Minimum total counts per cell to pass QC. Default: 20000",
    )
    parser.add_argument(
        "--min-features",
        type=int,
        default=10000,
        help="Minimum detected anchor-target features per cell to pass QC. Default: 10000",
    )
    return parser.parse_args()


def load_whitelist(whitelist_path: Path) -> set[str]:
    whitelist: set[str] = set()
    with whitelist_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=1):
            if len(row) < 2:
                continue
            combo_sequence = row[1].strip().upper()
            if combo_sequence:
                whitelist.add(combo_sequence)
    if not whitelist:
        raise ValueError(f"[ERROR] No whitelist barcodes loaded from: {whitelist_path}")
    return whitelist


def load_cells(cells_path: Path) -> list[dict[str, str]]:
    with cells_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cells = list(reader)
    if not cells:
        raise ValueError(f"[ERROR] No cells loaded from: {cells_path}")
    return cells


def load_features(features_path: Path) -> list[dict[str, str]]:
    with features_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        features = list(reader)
    if not features:
        raise ValueError(f"[ERROR] No features loaded from: {features_path}")
    return features


def load_matrix(matrix_path: Path):
    if matrix_path.suffix == ".gz":
        with gzip.open(matrix_path, "rb") as handle:
            matrix = mmread(handle)
    else:
        matrix = mmread(matrix_path)
    return matrix.tocsr()


def classify_cells(
    matrix,
    cells: list[dict[str, str]],
    whitelist: set[str],
    min_total_count: int,
    min_features: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    detected_features = np.diff(matrix.indptr)
    in_whitelist = np.array(
        [cell["cell_barcode"].strip().upper() in whitelist for cell in cells],
        dtype=bool,
    )
    qc_pass = (total_counts >= min_total_count) & (detected_features >= min_features)
    in_whitelist_pass = in_whitelist & qc_pass
    in_whitelist_fail = in_whitelist & (~qc_pass)
    not_in_whitelist_fail = ~in_whitelist
    return total_counts, detected_features, in_whitelist_pass, in_whitelist_fail, not_in_whitelist_fail


def write_cells(path: Path, cells: list[dict[str, str]], keep_mask: np.ndarray) -> None:
    fieldnames = ["cell_id", "sample_id", "cell_barcode"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for keep, cell in zip(keep_mask, cells):
            if keep:
                writer.writerow({key: cell[key] for key in fieldnames})


def write_features(path: Path, features: list[dict[str, str]], keep_mask: np.ndarray) -> None:
    fieldnames = ["anchor_target", "anchor", "target"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for keep, feature in zip(keep_mask, features):
            if keep:
                writer.writerow({key: feature[key] for key in fieldnames})


def write_matrix(path: Path, matrix) -> None:
    with gzip.open(path, "wb") as handle:
        mmwrite(handle, matrix)


def plot_qc_scatter(
    total_counts: np.ndarray,
    detected_features: np.ndarray,
    in_whitelist_pass: np.ndarray,
    in_whitelist_fail: np.ndarray,
    not_in_whitelist_fail: np.ndarray,
    output_path: Path,
    args
) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))

    categories = [
        (not_in_whitelist_fail, "not_in_whitelist_fail", "#bdbdbd", 6, 0.45),
        (in_whitelist_fail, "in_whitelist_fail", "#d95f02", 8, 0.65),
        (in_whitelist_pass, "in_whitelist_pass", "#1b9e77", 8, 0.85),
    ]

    for mask, label, color, size, alpha in categories:
        if np.any(mask):
            ax.scatter(
                total_counts[mask],
                detected_features[mask],
                s=size,
                c=color,
                alpha=alpha,
                linewidths=0,
                label=f"{label} (n={int(mask.sum())})",
            )

    ax.axvline(x=args.min_total_count, color="black", linestyle="--", linewidth=1.2)
    ax.axhline(y=args.min_features, color="red", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(labelsize=8)
    ax.set_xlabel("Total anchor-target counts per cell", fontsize=10)
    ax.set_ylabel("Detected anchor-target features per cell", fontsize=10)
    ax.set_title("Cell Filtering by Whitelist and QC", fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=True, fontsize=8, borderpad=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_qc_histograms(
    total_counts: np.ndarray,
    detected_features: np.ndarray,
    output_path: Path,
    args,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    total_vals = np.asarray(total_counts, dtype=float)
    total_vals = total_vals[np.isfinite(total_vals)]
    total_vals = np.maximum(total_vals, 0.0)
    axes[0].hist(total_vals, bins=80, color="#1F77B4", alpha=0.8, edgecolor="white", linewidth=0.5)
    axes[0].axvline(args.min_total_count, color="black", linestyle="--", linewidth=1.2, label=f"cutoff={args.min_total_count}")
    axes[0].tick_params(labelsize=10)
    axes[0].set_xlabel("Total anchor-target counts per cell", fontsize=12)
    axes[0].set_ylabel("Cell count", fontsize=12)
    axes[0].set_title("Total counts distribution", fontsize=13, fontweight="bold")
    axes[0].legend(frameon=False)

    feature_vals = np.asarray(detected_features, dtype=float)
    feature_vals = feature_vals[np.isfinite(feature_vals)]
    feature_vals = np.maximum(feature_vals, 0.0)
    axes[1].hist(feature_vals, bins=80, color="#2CA02C", alpha=0.8, edgecolor="white", linewidth=0.5)
    axes[1].axvline(args.min_features, color="red", linestyle="--", linewidth=1.2, label=f"cutoff={args.min_features}")
    axes[1].tick_params(labelsize=10)
    axes[1].set_xlabel("Detected anchor-target features per cell", fontsize=12)
    axes[1].set_ylabel("Cell count", fontsize=12)
    axes[1].set_title("Detected features distribution", fontsize=13, fontweight="bold")
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    matrix_path = args.out_dir / "cell_anchor_target_counts.matrix.mtx.gz"
    cells_path = args.out_dir / "cell_anchor_target_counts.cells.tsv"
    features_path = args.out_dir / "cell_anchor_target_counts.features.tsv"

    os.makedirs(args.out_dir / "filtered", exist_ok=True)

    filtered_matrix_path = args.out_dir / "filtered" / "whitelist_filtered.matrix.mtx.gz"
    filtered_cells_path = args.out_dir / "filtered" / "whitelist_filtered.cells.tsv"
    filtered_features_path = args.out_dir / "filtered" / "whitelist_filtered.features.tsv"
    scatter_plot_path = args.out_dir / "filtered" / "whitelist_filter_scatter.png"
    hist_plot_path = args.out_dir / "filtered" / "whitelist_filter_hist.png"

    if not args.whitelist.is_file():
        raise FileNotFoundError(f"[ERROR] Whitelist file not found: {args.whitelist}")
    if not matrix_path.is_file():
        raise FileNotFoundError(f"[ERROR] Matrix file not found: {matrix_path}")
    if not cells_path.is_file():
        raise FileNotFoundError(f"[ERROR] Cells file not found: {cells_path}")
    if not features_path.is_file():
        raise FileNotFoundError(f"[ERROR] Features file not found: {features_path}")

    print(f"[INFO] Loading whitelist from {args.whitelist}")
    whitelist = load_whitelist(args.whitelist)

    print(f"[INFO] Loading cells from {cells_path}")
    cells = load_cells(cells_path)
    print(f"[INFO] Loaded {len(cells)} cells")

    print(f"[INFO] Loading features from {features_path}")
    features = load_features(features_path)
    print(f"[INFO] Loaded {len(features)} features")

    print(f"[INFO] Loading matrix from {matrix_path}")
    matrix = load_matrix(matrix_path)

    if matrix.shape[0] != len(cells):
        raise ValueError(f"[ERROR] Matrix row count {matrix.shape[0]} does not match cells count {len(cells)}")
    if matrix.shape[1] != len(features):
        raise ValueError(f"[ERROR] Matrix column count {matrix.shape[1]} does not match features count {len(features)}")

    total_counts, detected_features, in_whitelist_pass, in_whitelist_fail, not_in_whitelist_fail = classify_cells(
        matrix,
        cells,
        whitelist,
        args.min_total_count,
        args.min_features,
    )

    kept_matrix = matrix[in_whitelist_pass]
    kept_feature_mask = np.asarray(kept_matrix.sum(axis=0)).ravel() > 0
    kept_matrix = kept_matrix[:, kept_feature_mask]

    print(f"[INFO] Writing filtered cells to {filtered_cells_path}")
    write_cells(filtered_cells_path, cells, in_whitelist_pass)
    print(f"[INFO] Writing filtered features to {filtered_features_path}")
    write_features(filtered_features_path, features, kept_feature_mask)
    print(f"[INFO] Writing filtered matrix to {filtered_matrix_path}")
    write_matrix(filtered_matrix_path, kept_matrix)

    print(f"[INFO] Plotting QC scatter to {scatter_plot_path}")
    plot_qc_scatter(
        total_counts,
        detected_features,
        in_whitelist_pass,
        in_whitelist_fail,
        not_in_whitelist_fail,
        scatter_plot_path,
        args
    )

    print(f"[INFO] Plotting QC histograms to {hist_plot_path}")
    plot_qc_histograms(
        total_counts,
        detected_features,
        hist_plot_path,
        args,
    )


if __name__ == "__main__":
    main()
