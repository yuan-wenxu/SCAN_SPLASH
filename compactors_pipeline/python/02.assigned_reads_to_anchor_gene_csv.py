#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert assigned_reads txt to CSV with anchor(31bp), flag, chr, pos, gene_id, gene_name; merge close pos within gene_id"
    )
    parser.add_argument("--input", type=Path, required=True, help="Path to assigned_reads txt")
    parser.add_argument("--output", type=Path, required=True, help="Path to output CSV")
    parser.add_argument(
        "--merge-distance",
        type=int,
        default=15,
        help="Merge positions within the same gene_id when pos difference to group start <= this value (default: 15, non-transitive)",
    )
    parser.add_argument(
        "--gtf-gene-id-attr",
        type=str,
        default="gene_id",
        help="GTF attribute key for gene ID (default: gene_id)",
    )
    parser.add_argument(
        "--gtf-gene-name-attr",
        type=str,
        default="gene_name",
        help="GTF attribute key for gene name (default: gene_name)",
    )
    return parser.parse_args()


def load_gene_name_map(gtf_path: Path, gene_id_attr: str = "gene_id", gene_name_attr: str = "gene_name"):
    gene_id_re = re.compile(rf'{gene_id_attr} "([^"]+)"')
    gene_name_re = re.compile(rf'{gene_name_attr} "([^"]+)"')
    gene_name_map = {}
    with gtf_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue

            attrs = fields[8]
            gid_match = gene_id_re.search(attrs)
            if not gid_match:
                continue

            gene_id = gid_match.group(1)
            gname_match = gene_name_re.search(attrs)
            gene_name = gname_match.group(1) if gname_match else ""
            gene_name_map[gene_id] = gene_name

    return gene_name_map


def parse_gene_id(last_field: str):
    if last_field.startswith("GX:Z:"):
        return last_field[5:]
    elif last_field.startswith("GN:Z:"):
        return last_field[5:]
    return last_field


def parse_flag(flag_text: str):
    try:
        return int(flag_text)
    except (TypeError, ValueError):
        return 0


def extract_anchor_from_read_id(read_id: str) -> str:
    # Preferred read ID format:
    #   <prefix>_<uuid>_<rand>_<anchor>_pvalOpt_<pval>_effectSizeBin_<effect_size>
    # Parse using marker tokens to avoid breakage when pval contains '-' or '.'.
    marker = "_pvalOpt_"
    marker2 = "_effectSizeBin_"

    if marker in read_id and marker2 in read_id:
        left, right = read_id.split(marker, 1)
        if marker2 in right:
            pval, effect_size = right.split(marker2, 1)
            anchor = left.rsplit("_", 1)[-1] if "_" in left else left
            return anchor, pval, effect_size

    # Backward-compatible fallback
    parts = read_id.split("_")
    if len(parts) > 5:
        return parts[-5], parts[-3], parts[-1]
    return "", "", ""


def parse_pos(pos_text: str):
    try:
        return int(pos_text)
    except (TypeError, ValueError):
        return None


def merge_records_non_transitive(records, merge_distance: int):
    grouped = {}
    for rec in records:
        grouped.setdefault(rec["gene_id"], []).append(rec)

    merged = []
    for items in grouped.values():
        valid = sorted(
            [r for r in items if r["pos_int"] is not None],
            key=lambda x: x["pos_int"],
        )
        merged.extend(r for r in items if r["pos_int"] is None)

        i = 0
        while i < len(valid):
            base_pos = valid[i]["pos_int"]
            merged.append(valid[i])
            j = i + 1
            while j < len(valid) and valid[j]["pos_int"] - base_pos <= merge_distance:
                j += 1
            i = j

    merged.sort(key=lambda x: (x["gene_id"], x["chr"], x["pos_int"] if x["pos_int"] is not None else -1, x["anchor"]))
    return merged


def main():
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"[ERROR] Input file not found: {args.input}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    records = []
    with args.input.open("r", encoding="utf-8") as fin:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            read_id = fields[0] if fields else ""
            flag_text = fields[1] if len(fields) >= 2 else "0"
            flag = parse_flag(flag_text)
            chr_name = fields[2] if len(fields) >= 3 else ""
            pos = fields[3] if len(fields) >= 4 else ""
            anchor, pval, effect_size = extract_anchor_from_read_id(read_id)
            seq = fields[9] if len(fields) >= 10 else ""

            gene_id_field = fields[-2] if fields else ""
            gene_name_field = fields[-1] if fields else ""
            gene_id = parse_gene_id(gene_id_field)
            gene_name = parse_gene_id(gene_name_field)

            records.append(
                {
                    "anchor": anchor,
                    "pval": pval,
                    "effect_size": effect_size,
                    "flag": flag,
                    "chr": chr_name,
                    "pos": pos,
                    "pos_int": parse_pos(pos),
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "seq": seq,
                }
            )

    merged_records = merge_records_non_transitive(records, args.merge_distance)

    with args.output.open("w", encoding="utf-8", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["anchor", "pval", "effect_size", "flag", "chr", "pos", "gene_id", "gene_name", "seq"])
        for rec in merged_records:
            writer.writerow([rec["anchor"], rec["pval"], rec["effect_size"], rec["flag"], rec["chr"], rec["pos"], rec["gene_id"], rec["gene_name"], rec["seq"]])


if __name__ == "__main__":
    main()
