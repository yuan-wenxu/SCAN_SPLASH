#!/usr/bin/env python3
"""Filter cell x anchor Cj matrix in a single pass.

Input files (under --matrix-dir):
- cell_anchor_cj.matrix.mtx.gz
- cell_anchor_cj.cells.tsv
- cell_anchor_cj.features.tsv

Filtering strategy (single-pass):
1. Mark cells as pass/fail by --min-anchors-per-cell.
2. Mark cells as in/out by whitelist barcode membership.
3. Keep only cells that are both in whitelist and pass filter.
4. Keep anchors observed in >= --min-cells-per-anchor among kept cells.

Output files (under --matrix-dir/filtered):
- cell_anchor_cj.filtered.matrix.mtx.gz
- cell_anchor_cj.filtered.cells.tsv
- cell_anchor_cj.filtered.features.tsv
- filter_cell_anchor_cj_matrix.summary.tsv
- filter_cell_anchor_cj_scatter.png
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread, mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter cell x anchor Cj matrix by nonzero support"
    )
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument(
        "--whitelist",
        type=Path,
        required=True,
        help="Whitelist TSV with combo_sequence in column 2",
    )
    parser.add_argument(
        "--min-anchors-per-cell",
        type=int,
        default=2000,
        help="Keep cells with at least this many nonzero anchors (default: 2000)",
    )
    parser.add_argument(
        "--min-cells-per-anchor",
        type=int,
        default=3,
        help="Keep anchors observed in at least this many cells (default: 3)",
    )
    parser.add_argument(
        "--scatter-png",
        type=Path,
        default=None,
        help="Output scatter PNG path (default: <matrix-dir>/filtered/filter_cell_anchor_cj_scatter.png)",
    )
    return parser.parse_args()


def load_whitelist(path: Path) -> set[str]:
    whitelist: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            combo = row[1].strip().upper()
            if combo:
                whitelist.add(combo)
    if not whitelist:
        raise ValueError(f"[ERROR] No whitelist barcodes loaded from: {path}")
    return whitelist


def load_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def load_matrix(path: Path):
    with gzip.open(path, "rb") as fh:
        m = mmread(fh)
    return m.tocsr()


def save_matrix(path: Path, matrix) -> None:
    with gzip.open(path, "wb") as fh:
        mmwrite(fh, matrix)


def plot_scatter(
    x_total_abs: np.ndarray,
    y_nnz: np.ndarray,
    in_whitelist: np.ndarray,
    pass_filter: np.ndarray,
    output_png: Path,
    min_anchors_per_cell: int,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))

    categories = [
        ((~in_whitelist) & (~pass_filter), "out_whitelist_fail", "#bdbdbd", 10, 0.45),
        ((~in_whitelist) & pass_filter, "out_whitelist_pass", "#7570b3", 12, 0.65),
        (in_whitelist & (~pass_filter), "in_whitelist_fail", "#d95f02", 14, 0.65),
        (in_whitelist & pass_filter, "in_whitelist_pass", "#1b9e77", 16, 0.85),
    ]

    for mask, label, color, size, alpha in categories:
        if np.any(mask):
            ax.scatter(
                x_total_abs[mask],
                y_nnz[mask],
                s=size,
                c=color,
                alpha=alpha,
                linewidths=0,
                label=f"{label} (n={int(mask.sum())})",
            )

    ax.axhline(y=min_anchors_per_cell, color="black", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total absolute Cj per cell")
    ax.set_ylabel("Detected anchors per cell (nonzero)")
    ax.set_title("Cell filtering by whitelist and nonzero-anchor threshold")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    args = parse_args()

    matrix_path = args.matrix_dir / "cell_anchor_cj.matrix.mtx.gz"
    cells_path = args.matrix_dir / "cell_anchor_cj.cells.tsv"
    features_path = args.matrix_dir / "cell_anchor_cj.features.tsv"

    if not matrix_path.is_file():
        raise FileNotFoundError(f"[ERROR] Missing matrix file: {matrix_path}")
    if not cells_path.is_file():
        raise FileNotFoundError(f"[ERROR] Missing cells file: {cells_path}")
    if not features_path.is_file():
        raise FileNotFoundError(f"[ERROR] Missing features file: {features_path}")
    if not args.whitelist.is_file():
        raise FileNotFoundError(f"[ERROR] Missing whitelist file: {args.whitelist}")

    print(f"[INFO] Loading matrix: {matrix_path}")
    matrix = load_matrix(matrix_path)
    cells = load_tsv_rows(cells_path)
    features = load_tsv_rows(features_path)
    whitelist = load_whitelist(args.whitelist)

    if matrix.shape[0] != len(cells):
        raise ValueError(f"[ERROR] Matrix rows ({matrix.shape[0]}) != cells ({len(cells)})")
    if matrix.shape[1] != len(features):
        raise ValueError(f"[ERROR] Matrix cols ({matrix.shape[1]}) != features ({len(features)})")

    cell_nnz = np.asarray((matrix != 0).sum(axis=1)).ravel()
    cell_total_abs = np.asarray(np.abs(matrix).sum(axis=1)).ravel()
    in_whitelist = np.array(
        [str(c.get("cell_barcode", "")).strip().upper() in whitelist for c in cells],
        dtype=bool,
    )
    pass_filter = cell_nnz >= args.min_anchors_per_cell

    keep_cells = in_whitelist & pass_filter
    sub = matrix[keep_cells]
    feat_nnz = np.asarray((sub != 0).sum(axis=0)).ravel()
    keep_feats = feat_nnz >= args.min_cells_per_anchor

    filtered = sub[:, keep_feats]
    kept_cells = [c for i, c in enumerate(cells) if keep_cells[i]]
    kept_features = [f for i, f in enumerate(features) if keep_feats[i]]

    out_dir = args.matrix_dir / "filtered"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_matrix = out_dir / "cell_anchor_cj.filtered.matrix.mtx.gz"
    out_cells = out_dir / "cell_anchor_cj.filtered.cells.tsv"
    out_features = out_dir / "cell_anchor_cj.filtered.features.tsv"
    out_summary = out_dir / "filter_cell_anchor_cj_matrix.summary.tsv"
    out_scatter = args.scatter_png if args.scatter_png else (out_dir / "filter_cell_anchor_cj_scatter.png")

    save_matrix(out_matrix, filtered)

    cell_fields = ["cell_id", "sample_id", "cell_barcode"]
    feat_fields = ["feature_id", "anchor"]
    write_tsv_rows(out_cells, kept_cells, cell_fields)
    write_tsv_rows(out_features, kept_features, feat_fields)

    plot_scatter(
        x_total_abs=cell_total_abs,
        y_nnz=cell_nnz,
        in_whitelist=in_whitelist,
        pass_filter=pass_filter,
        output_png=out_scatter,
        min_anchors_per_cell=args.min_anchors_per_cell,
    )

    with out_summary.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["n_cells_before", len(cells)])
        writer.writerow(["n_anchors_before", len(features)])
        writer.writerow(["nnz_before", matrix.nnz])
        writer.writerow(["n_cells_after", len(kept_cells)])
        writer.writerow(["n_anchors_after", len(kept_features)])
        writer.writerow(["nnz_after", filtered.nnz])
        writer.writerow(["min_anchors_per_cell", args.min_anchors_per_cell])
        writer.writerow(["min_cells_per_anchor", args.min_cells_per_anchor])
        writer.writerow(["n_in_whitelist", int(in_whitelist.sum())])
        writer.writerow(["n_pass_filter", int(pass_filter.sum())])
        writer.writerow(["n_in_whitelist_pass", int((in_whitelist & pass_filter).sum())])
        writer.writerow(["n_in_whitelist_fail", int((in_whitelist & (~pass_filter)).sum())])
        writer.writerow(["n_out_whitelist_pass", int(((~in_whitelist) & pass_filter).sum())])
        writer.writerow(["n_out_whitelist_fail", int(((~in_whitelist) & (~pass_filter)).sum())])

    print(f"[OK] n_cells_before={len(cells)} n_cells_after={len(kept_cells)}")
    print(f"[OK] n_anchors_before={len(features)} n_anchors_after={len(kept_features)}")
    print(f"[OK] nnz_before={matrix.nnz} nnz_after={filtered.nnz}")
    print(f"[OK] filtered_matrix={out_matrix}")
    print(f"[OK] scatter_png={out_scatter}")


if __name__ == "__main__":
    main()
