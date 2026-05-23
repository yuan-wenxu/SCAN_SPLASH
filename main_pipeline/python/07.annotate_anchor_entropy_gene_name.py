#!/usr/bin/env python3
"""
Annotate an anchor entropy CSV with gene names from a GTF file.

Input CSV is expected to contain a gene_id column, for example:
anchor,gene_id,entropy,n_targets,total_count

Output CSV appends a gene_name column while preserving row order.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s+"([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add gene_name to an anchor entropy CSV using a GTF file"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input anchor entropy CSV containing a gene_id column",
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        required=True,
        help="Reference GTF file containing gene_id and gene_name attributes",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (default: <input stem>.with_gene_name.csv)",
    )
    return parser.parse_args()


def parse_attributes(attr_text: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE_RE.findall(attr_text)}


def load_gene_names(gtf_path: Path) -> dict[str, str]:
    gene_names: dict[str, str] = {}
    with open(gtf_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            if gene_id and gene_id not in gene_names:
                gene_names[gene_id] = gene_name
    return gene_names


def main() -> None:
    args = parse_args()
    output_path = args.output
    if output_path is None:
        output_path = args.input.with_name(f"{args.input.stem}.with_gene_name.csv")

    gene_names = load_gene_names(args.gtf)

    with open(args.input, "r", encoding="utf-8", newline="") as in_fh:
        reader = csv.DictReader(in_fh)
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames:
            raise ValueError("[ERROR] Input CSV must contain a gene_id column")

        output_fields = list(reader.fieldnames)
        if "gene_name" not in output_fields:
            output_fields.append("gene_name")

        with open(output_path, "w", encoding="utf-8", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=output_fields)
            writer.writeheader()

            for row in reader:
                gene_id = row.get("gene_id", "")
                row["gene_name"] = gene_names.get(gene_id, "")
                writer.writerow(row)

    print(f"[INFO] Annotated CSV written to: {output_path}")
    print(f"[INFO] Loaded {len(gene_names)} gene_id -> gene_name mappings from {args.gtf}")


if __name__ == "__main__":
    main()