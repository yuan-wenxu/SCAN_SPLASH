#!/usr/bin/env python3
"""Convert SPLASH compactors TSV to FASTQ.GZ.

Input TSV is expected to include a column named 'compactor'.
Each non-empty compactor sequence becomes one FASTQ read.
- Read ID: random UUID-like string with optional prefix.
- Quality: fixed single character repeated to sequence length.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import uuid
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert compactors TSV to FASTQ.GZ with random read IDs and fixed quality"
    )
    parser.add_argument("--input-tsv", type=Path, required=True, help="Input compactors TSV")
    parser.add_argument(
        "--scores-tsv",
        type=Path,
        required=True,
        help="Input scores TSV used to annotate read IDs with anchor-level pval_opt and effect_size_bin",
    )
    parser.add_argument("--output-fastq-gz", type=Path, required=True, help="Output FASTQ.GZ path")
    parser.add_argument(
        "--compactor-column",
        type=str,
        default="compactor",
        help="Column name that stores sequence (default: compactor)",
    )
    parser.add_argument(
        "--quality-char",
        type=str,
        default="I",
        help="Single fixed quality character (default: I)",
    )
    parser.add_argument(
        "--readid-prefix",
        type=str,
        default="cmp",
        help="Prefix for random read IDs (default: cmp)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible read IDs",
    )
    return parser.parse_args()


def normalize_seq(seq: str) -> str:
    seq = seq.strip().upper()
    # Keep only canonical bases and N; drop others to avoid invalid FASTQ bases.
    return "".join(ch for ch in seq if ch in {"A", "C", "G", "T", "N"})


def make_read_id(prefix: str) -> str:
    return f"{prefix}_{uuid.uuid4().hex}_{random.getrandbits(32):08x}"


def sanitize_readid_part(value: object) -> str:
    text = str(value).strip()
    if not text:
        return "NA"
    # Keep scientific notation intact (e.g. 5.4e-08) and only replace risky separators.
    return text.replace("\t", "_").replace(" ", "_")


def load_anchor_scores(scores_tsv: Path) -> dict[str, tuple[str, str]]:
    anchor_scores: dict[str, tuple[str, str]] = {}
    with scores_tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("[ERROR] Scores TSV has no header")

        field_map = {field.strip(): field for field in reader.fieldnames}
        if "anchor" not in field_map:
            raise ValueError(f"[ERROR] Column 'anchor' not found. Available: {reader.fieldnames}")
        if "pval_opt" not in field_map:
            raise ValueError(f"[ERROR] Column 'pval_opt' not found. Available: {reader.fieldnames}")
        if "effect_size_bin" not in field_map:
            raise ValueError(f"[ERROR] Column 'effect_size_bin' not found. Available: {reader.fieldnames}")

        anchor_field = field_map["anchor"]
        pval_field = field_map["pval_opt"]
        effect_field = field_map["effect_size_bin"]

        for row in reader:
            anchor = (row.get(anchor_field) or "").strip()
            if not anchor:
                continue
            pval_opt = sanitize_readid_part(row.get(pval_field, "NA"))
            effect_size_bin = sanitize_readid_part(row.get(effect_field, "NA"))
            if anchor not in anchor_scores:
                anchor_scores[anchor] = (pval_opt, effect_size_bin)

    return anchor_scores


def main() -> None:
    args = parse_args()

    if not args.input_tsv.is_file():
        raise FileNotFoundError(f"[ERROR] Input TSV not found: {args.input_tsv}")
    if not args.scores_tsv.is_file():
        raise FileNotFoundError(f"[ERROR] Scores TSV not found: {args.scores_tsv}")

    if len(args.quality_char) != 1:
        raise ValueError("[ERROR] --quality-char must be exactly one character")

    if args.seed is not None:
        random.seed(args.seed)

    args.output_fastq_gz.parent.mkdir(parents=True, exist_ok=True)
    anchor_scores = load_anchor_scores(args.scores_tsv)

    total_rows = 0
    written_reads = 0
    empty_or_invalid = 0
    missing_scores = 0

    with args.input_tsv.open("r", encoding="utf-8", newline="") as in_fh, gzip.open(
        args.output_fastq_gz, "wt", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("[ERROR] Input TSV has no header")
        if args.compactor_column not in reader.fieldnames:
            raise ValueError(
                f"[ERROR] Column '{args.compactor_column}' not found. Available: {reader.fieldnames}"
            )

        for row in reader:
            total_rows += 1
            raw_seq = row.get(args.compactor_column, "")
            anchor = row.get("anchor", "NA")
            seq = normalize_seq(raw_seq)
            if not seq:
                empty_or_invalid += 1
                continue

            read_id = make_read_id(args.readid_prefix)
            pval_opt, effect_size_bin = anchor_scores.get(anchor, ("NA", "NA"))
            if anchor not in anchor_scores:
                missing_scores += 1
            read_id = f"{read_id}_{anchor}_pvalOpt_{pval_opt}_effectSizeBin_{effect_size_bin}"
            qual = args.quality_char * len(seq)

            out_fh.write(f"@{read_id}\n")
            out_fh.write(f"{seq}\n")
            out_fh.write("+\n")
            out_fh.write(f"{qual}\n")
            written_reads += 1

    print(f"[OK] input={args.input_tsv}")
    print(f"[OK] output={args.output_fastq_gz}")
    print(f"[OK] total_rows={total_rows}")
    print(f"[OK] written_reads={written_reads}")
    print(f"[OK] empty_or_invalid={empty_or_invalid}")
    print(f"[OK] missing_scores={missing_scores}")


if __name__ == "__main__":
    main()
