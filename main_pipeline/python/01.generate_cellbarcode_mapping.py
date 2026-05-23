#!/usr/bin/env python3
"""Generate barcode mappings with large pairwise distances.

Input: headerless TSV, col1=combo_name, col2=combo_sequence, col3=well.
Output CSV columns: combo_name, combo_sequence, well, mapped_sequence.
"""

from __future__ import annotations

import argparse
import csv
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import matplotlib.pyplot as plt


DNA_ALPHABET = "ACGT"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate 16bp barcode mapping CSV with large Hamming distance")
    parser.add_argument("--input", type=Path, help="Input CSV/TSV with combo name and combo sequence")
    parser.add_argument("--output", type=Path, help="Output CSV path")
    parser.add_argument("--code-len", type=int, default=16, help="Length of mapped sequence")
    parser.add_argument("--pool-size", type=int, default=800, help="Candidate pool size per barcode; larger gives better distances")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible mapping")
    return parser.parse_args()


@dataclass(frozen=True)
class Combo:
    name: str
    sequence: str
    well: str


def detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".tsv":
        return "\t"
    if path.suffix.lower() == ".csv":
        return ","
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    return "\t" if sample.count("\t") > sample.count(",") else ","


def looks_like_dna(seq: str) -> bool:
    s = seq.strip().upper()
    return bool(s) and all(base in DNA_ALPHABET for base in s)


def load_combos(path: Path) -> list[Combo]:
    delimiter = detect_delimiter(path)
    combos: list[Combo] = []
    seen_names: set[str] = set()

    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter=delimiter))

    if not rows:
        raise ValueError(f"[ERROR] Input file is empty: {path}")

    for i, row in enumerate(rows, start=1):
        if len(row) != 3:
            continue
        name = row[0].strip()
        seq = row[1].strip().upper()
        well = row[2].strip()

        if not name or not seq:
            continue
        if not looks_like_dna(seq):
            raise ValueError(f"[ERROR] Invalid DNA bases in line {i}: {seq}")
        if name in seen_names:
            raise ValueError(f"[ERROR] Duplicate combo name '{name}' in input")
        seen_names.add(name)
        combos.append(Combo(name=name, sequence=seq, well=well))

    if not combos:
        raise ValueError(f"[ERROR] No valid combos loaded from {path}")
    return combos


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def random_barcode(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA_ALPHABET) for _ in range(length))


def choose_next_code(
    selected: list[str],
    code_len: int,
    pool_size: int,
    rng: random.Random,
) -> str:
    max_possible = code_len
    best_code = ""
    best_min_dist = -1

    tries = 0
    while tries < pool_size:
        tries += 1
        cand = random_barcode(code_len, rng)
        if cand in selected:
            continue

        min_dist = max_possible
        for s in selected:
            d = hamming(cand, s)
            if d < min_dist:
                min_dist = d

        if min_dist > best_min_dist or (min_dist == best_min_dist):
            best_code = cand
            best_min_dist = min_dist
            if best_min_dist >= 4:  # early stop if we find a code with distance >= 4
                break

    if not best_code:
        while True:
            cand = random_barcode(code_len, rng)
            if cand not in selected:
                return cand
    return best_code


def build_mapping(combos: Iterable[Combo], code_len: int, pool_size: int, seed: int) -> list[tuple[Combo, str]]:
    combo_list = list(combos)
    if len(combo_list) > (4**code_len):
        raise ValueError(f"[ERROR] Need {len(combo_list)} unique codes but only {4**code_len} possible for length={code_len}")

    rng = random.Random(seed)
    selected: list[str] = []
    out: list[tuple[Combo, str]] = []

    for combo in combo_list:
        if not selected:
            code = random_barcode(code_len, rng)
            selected.append(code)
            out.append((combo, code))
            continue

        code = choose_next_code(
            selected=selected,
            code_len=code_len,
            pool_size=pool_size,
            rng=rng,
        )
        selected.append(code)
        out.append((combo, code))

    return out


def pairwise_distance_stats(codes: list[str]) -> tuple[int, float, int]:
    if len(codes) < 2:
        return (0, 0.0, 0)
    dists: list[int] = []
    for i in range(len(codes)):
        for j in range(i + 1, len(codes)):
            dists.append(hamming(codes[i], codes[j]))
    return dists


def plot_distance_histogram(dists: list[int], output_path: Path) -> None:
    plt.figure(figsize=(5, 4))
    plt.hist(dists, bins=range(min(dists), max(dists) + 2), align="left", rwidth=0.8)
    plt.xlabel("Hamming Distance")
    plt.ylabel("Count")
    plt.title("Pairwise Hamming Distance Distribution")
    plt.xticks(range(min(dists), max(dists) + 1))
    plt.grid(axis="y", alpha=0.75)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")


def write_mapping_csv(path: Path, rows: list[tuple[Combo, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["combo_name", "combo_sequence", "well", "mapped_sequence"])
        for combo, mapped in rows:
            writer.writerow([combo.name, combo.sequence, combo.well, mapped])


def main() -> None:
    args = parse_args()

    if args.code_len <= 0:
        raise ValueError("[ERROR] --code-len must be > 0")
    if args.pool_size <= 0:
        raise ValueError("[ERROR] --pool-size must be > 0")
    if not args.input.is_file():
        raise ValueError(f"[ERROR] Input file not found: {args.input}")

    print(f"[INFO] Loading combos from {args.input}")
    combos = load_combos(args.input)
    print(f"[INFO] Loaded {len(combos)} combos. Generating mappings with code_len={args.code_len}, pool_size={args.pool_size}, seed={args.seed}")
    print(f"[INFO] Total possible codes of length {args.code_len}: {4**args.code_len}")
    print(f"[INFO] Starting mapping generation...")
    mapped_rows = build_mapping(
        combos=combos,
        code_len=args.code_len,
        pool_size=args.pool_size,
        seed=args.seed,
    )
    write_mapping_csv(args.output, mapped_rows)

    codes = [x[1] for x in mapped_rows]
    dists = pairwise_distance_stats(codes)
    output_hist = args.output.with_suffix(".distance_histogram.png")
    plot_distance_histogram(dists, output_hist)
    print(f"[DONE] combos={len(mapped_rows)} output={args.output} distance_histogram={output_hist}")


if __name__ == "__main__":
    main()