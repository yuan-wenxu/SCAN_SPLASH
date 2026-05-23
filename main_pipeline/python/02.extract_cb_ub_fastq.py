#!/usr/bin/env python3
"""Build pseudo-R1 FASTQ from CB/UB tags stored in R2 read IDs.

Input header format supported:
    @read1+CB:ACGT+UB:TTAA

Output:
    FASTQ whose sequence is CB+UB (or mapped CB when mapping is provided)
    and whose quality is a pseudo quality string.
    The original R2 sequence and quality are not modified because this script only
    writes the pseudo-R1 file.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from typing import Optional, TextIO, Tuple


def _open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def _parse_header(h: str) -> Tuple[str, str, str]:
    if not h.startswith("@"):
        raise ValueError(f"[ERROR] Invalid FASTQ header (does not start with @): {h.rstrip()}")

    content = h[1:].rstrip("\n")

    parts = content.split("+")
    if len(parts) != 3:
        raise ValueError(
            f"[ERROR] Invalid header format, expected '@<readid>+CB:<cellbarcode>+UB:<umi>': {h.rstrip()}"
        )

    read_name, cb_part, ub_part = parts

    if not cb_part.startswith("CB:"):
        raise ValueError(f"[ERROR] Invalid CB field in header: {h.rstrip()}")
    if not ub_part.startswith("UB:"):
        raise ValueError(f"[ERROR] Invalid UB field in header: {h.rstrip()}")

    cb = cb_part[len("CB:"):]
    ub = ub_part[len("UB:"):]

    if cb == "":
        raise ValueError(f"[ERROR] Empty CB value in header: {h.rstrip()}")
    if ub == "":
        raise ValueError(f"[ERROR] Empty UB value in header: {h.rstrip()}")

    return read_name, cb, ub


def load_barcode_mapping(mapping_csv: str) -> dict[str, str]:
    """Load mapping table as combo_sequence -> mapped_sequence.

    Supported formats:
    - Header CSV with columns including combo_sequence,mapped_sequence
      (also supports extra columns like well).
    - Headerless CSV where col2=combo_sequence, last_col=mapped_sequence.
    """
    mapping: dict[str, str] = {}
    with _open_text(mapping_csv, "r") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        raise ValueError(f"[ERROR] Empty mapping CSV: {mapping_csv}")

    first = [x.strip().lower() for x in rows[0]]
    has_header = ("combo_sequence" in first and "mapped_sequence" in first)

    if has_header:
        header = first
        seq_idx = header.index("combo_sequence")
        mapped_idx = header.index("mapped_sequence")
        data_rows = rows[1:]
    else:
        if len(rows[0]) < 3:
            raise ValueError(f"[ERROR] Invalid mapping CSV format: {mapping_csv}")
        seq_idx = 1
        mapped_idx = len(rows[0]) - 1
        data_rows = rows

    for line_no, parts in enumerate(data_rows, start=2 if has_header else 1):
        if len(parts) <= max(seq_idx, mapped_idx):
            continue
        combo_seq = parts[seq_idx].strip().upper()
        mapped_seq = parts[mapped_idx].strip().upper()
        if not combo_seq or not mapped_seq:
            continue
        if combo_seq in mapping and mapping[combo_seq] != mapped_seq:
            raise ValueError(f"[ERROR] Conflicting mapping for combo_sequence at line {line_no}: {combo_seq}")
        mapping[combo_seq] = mapped_seq

    if not mapping:
        raise ValueError(f"[ERROR] No valid mapping entries loaded from: {mapping_csv}")
    return mapping


def run(
    input_fastq: str,
    mapping_csv: Optional[str],
    output_fastq: str,
    qual_char: str,
) -> None:
    reads = 0
    mapped_reads = 0

    mapping: Optional[dict[str, str]] = None
    if mapping_csv:
        mapping = load_barcode_mapping(mapping_csv)

    os.makedirs(os.path.dirname(os.path.abspath(output_fastq)), exist_ok=True)

    with _open_text(input_fastq, "r") as fin, _open_text(output_fastq, "w") as fout:
        while True:
            h = fin.readline()
            if h == "":
                break

            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if seq == "" or plus == "" or qual == "":
                raise ValueError("[ERROR] Truncated FASTQ record detected at file end")

            read_name, cb, ub = _parse_header(h)

            if mapping is not None:
                mapped = mapping.get(cb)
                if mapped is None:
                    raise ValueError(f"[ERROR] CB combination not found in mapping: {cb}")
                combo_seq = mapped + ub
                mapped_reads += 1
            else:
                combo_seq = cb + ub

            r1_qual = qual_char * len(combo_seq)

            fout.write(f"@{read_name}\n")
            fout.write(f"{combo_seq}\n")
            fout.write("+\n")
            fout.write(f"{r1_qual}\n")
            reads += 1

    if mapping is None:
        print(f"[OK] reads={reads} input_r2={input_fastq} output_r1={output_fastq}")
    else:
        print(
            f"[OK] reads={reads} mapped_reads={mapped_reads} "
            f"[OK] input_r2={input_fastq} mapping_csv={mapping_csv} output_r1={output_fastq}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build pseudo-R1 FASTQ from CB/UB tags stored in R2 read IDs")
    parser.add_argument("--input", required=True, help="Input R2 FASTQ(.gz) containing CB/UB tags in read ID")
    parser.add_argument("--mapping-csv", help="Optional CSV mapping file to remap CB/UB combinations to new sequences")
    parser.add_argument("--output-fastq", required=True, help="Output pseudo-R1 FASTQ(.gz)")
    parser.add_argument("--qual-char", default="I", help="Single character used to build pseudo quality values for R1")
    args = parser.parse_args()

    if len(args.qual_char) != 1:
        print("[ERROR] --qual-char must be a single character", file=sys.stderr)
        sys.exit(1)

    try:
        if args.mapping_csv:
            print(f"[INFO] Loading mapping from {args.mapping_csv}")
            print(f"[INFO] Running with mapping...")
            run(args.input, args.mapping_csv, args.output_fastq, args.qual_char)
        else:
            print(f"[INFO] Running without mapping...")
            run(args.input, None, args.output_fastq, args.qual_char)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
