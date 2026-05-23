#!/usr/bin/env python3
"""
Compare targets within each anchor file and mark indels.

For each input anchor file:
- Use the first target as the reference sequence.
- Align every target to the reference.
- Represent deletions as '-' in aligned target.
- If a target is considered unalignable, treat all leading bases as missing:
  aligned_target = '-' * len(reference) + target
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from html import escape
from pathlib import Path
import sys


CPP_ALIGN_SRC = Path(__file__).resolve().parents[2] / "fast_align_cpp" / "src"

if CPP_ALIGN_SRC.is_dir() and str(CPP_ALIGN_SRC) not in sys.path:
    sys.path.insert(0, str(CPP_ALIGN_SRC))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare targets per anchor and output indel-formatted alignments"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing per-anchor target files (*.txt)"
    )
    parser.add_argument(
        "--pattern",
        type=str,
        default="*.txt",
        help="Glob pattern for input files (default: *.txt)"
    )
    parser.add_argument(
        "--html-out-dir",
        type=Path,
        default=None,
        help="Optional output directory for HTML visualization files"
    )
    parser.add_argument(
        "--min-common-substr",
        type=int,
        default=5,
        help="Minimum longest common substring length to consider alignable (default: 5)"
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.5,
        help="Minimum alignment identity to consider alignable (default: 0.5)"
    )
    parser.add_argument(
        "--max-gap-runs",
        type=int,
        default=2,
        help="Maximum number of gap runs allowed in pairwise alignment before fallback (default: 2)"
    )
    parser.add_argument(
        "--max-fragments",
        type=int,
        default=3,
        help="Maximum number of nucleotide fragments allowed before fallback (default: 3)"
    )
    parser.add_argument(
        "--max-errors",
        type=int,
        default=6,
        help="Maximum allowed errors per target (mismatch/indel) before fallback (default: 6)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of worker processes for parallel per-file processing (default: 8)"
    )
    return parser.parse_args()


def _python_longest_common_substring_len(a, b):
    """Return the length of the longest common substring."""
    if not a or not b:
        return 0

    prev = [0] * (len(b) + 1)
    best = 0
    for i, ca in enumerate(a, start=1):
        curr = [0]
        for j, cb in enumerate(b, start=1):
            if ca == cb:
                val = prev[j - 1] + 1
                curr.append(val)
                if val > best:
                    best = val
            else:
                curr.append(0)
        prev = curr
    return best


def _python_needleman_wunsch(ref, seq, match=2, mismatch=-1, gap=-2):
    """Global alignment for short sequences."""
    n = len(ref)
    m = len(seq)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    # trace: 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        trace[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        trace[0][j] = 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if ref[i - 1] == seq[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = diag
            t = 0
            if up > best:
                best = up
                t = 1
            if left > best:
                best = left
                t = 2

            score[i][j] = best
            trace[i][j] = t

    i, j = n, m
    aligned_ref = []
    aligned_seq = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == 0:
            aligned_ref.append(ref[i - 1])
            aligned_seq.append(seq[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace[i][j] == 1):
            aligned_ref.append(ref[i - 1])
            aligned_seq.append("-")
            i -= 1
        else:
            aligned_ref.append("-")
            aligned_seq.append(seq[j - 1])
            j -= 1

    aligned_ref.reverse()
    aligned_seq.reverse()

    return "".join(aligned_ref), "".join(aligned_seq), score[n][m]


longest_common_substring_len = _python_longest_common_substring_len
needleman_wunsch = _python_needleman_wunsch
ALIGNMENT_BACKEND = "python"

try:
    from scan_splash_align import (  # type: ignore
        longest_common_substring_len as _cpp_longest_common_substring_len,
        needleman_wunsch as _cpp_needleman_wunsch,
    )
except ImportError:
    pass
else:
    longest_common_substring_len = _cpp_longest_common_substring_len
    needleman_wunsch = _cpp_needleman_wunsch
    ALIGNMENT_BACKEND = "cpp"


def alignment_identity(aligned_ref, aligned_seq):
    if not aligned_ref:
        return 0.0
    matches = 0
    valid = 0
    for a, b in zip(aligned_ref, aligned_seq):
        if a == "-" and b == "-":
            continue
        valid += 1
        if a == b:
            matches += 1
    if valid == 0:
        return 0.0
    return matches / valid


def alignment_error_count(aligned_ref, aligned_seq):
    """Count total mismatches and indels for one pairwise alignment."""
    errors = 0
    for r_ch, t_ch in zip(aligned_ref, aligned_seq):
        if r_ch == t_ch:
            continue
        errors += 1
    return errors


def count_runs(seq, is_gap):
    runs = 0
    in_run = False
    for ch in seq:
        flag = is_gap(ch)
        if flag and not in_run:
            runs += 1
            in_run = True
        elif not flag:
            in_run = False
    return runs


def is_fragmented_alignment(aligned_target, max_gap_runs, max_fragments):
    gap_runs = count_runs(aligned_target, lambda c: c == "-")
    base_fragments = count_runs(aligned_target, lambda c: c != "-")
    return gap_runs > max_gap_runs or base_fragments > max_fragments


def parse_alignment_to_segments(aligned_ref, aligned_target, ref_len):
    insertions = {}
    bases = []
    ref_idx = 0

    for r_ch, t_ch in zip(aligned_ref, aligned_target):
        if r_ch == "-":
            insertions[ref_idx] = insertions.get(ref_idx, "") + t_ch
        else:
            bases.append(t_ch)
            ref_idx += 1

    if len(bases) < ref_len:
        bases.extend(["-"] * (ref_len - len(bases)))
    elif len(bases) > ref_len:
        bases = bases[:ref_len]

    return insertions, bases


def choose_reference_insertions(members, ref_len):
    ref_insertions = {}
    for bp in range(ref_len + 1):
        best = ""
        for m in members:
            insertions, _ = parse_alignment_to_segments(m["aligned_ref"], m["aligned_target"], ref_len)
            ins = insertions.get(bp, "")
            if len(ins) > len(best):
                best = ins
        if best:
            ref_insertions[bp] = best
    return ref_insertions


def build_projected_target(aligned_ref, aligned_target, ref_insertions, ref_len):
    insertions, bases = parse_alignment_to_segments(aligned_ref, aligned_target, ref_len)
    out = []
    for bp in range(ref_len + 1):
        slot_len = len(ref_insertions.get(bp, ""))
        if slot_len > 0:
            ins = insertions.get(bp, "")
            if ins:
                out.append(ins)
                if slot_len > len(ins):
                    out.append(" " * (slot_len - len(ins)))
            else:
                out.append(" " * slot_len)
        if bp < ref_len:
            out.append(bases[bp])
    return "".join(out)


def build_reference_with_gap_slots(reference, ref_insertions):
    out = []
    ref_len = len(reference)
    for bp in range(ref_len + 1):
        slot_len = len(ref_insertions.get(bp, ""))
        if slot_len > 0:
            out.append("-" * slot_len)
        if bp < ref_len:
            out.append(reference[bp])
    return "".join(out)


def render_colored_line(ref_line, target_line):
    spans = []
    for r_ch, t_ch in zip(ref_line, target_line):
        ch = escape(t_ch)
        if t_ch == " ":
            cls = "blank"
        elif t_ch == "-":
            cls = "gap"
        elif r_ch == "-":
            cls = "ins"
        elif t_ch == r_ch:
            cls = "match"
        else:
            cls = "sub"
        spans.append(f'<span class="{cls}">{ch}</span>')
    return "".join(spans)


def write_clustered_html_output(html_file, clusters):
    with open(html_file, "w") as out:
        out.write("<!doctype html>\n")
        out.write("<html><head><meta charset=\"utf-8\"><title>Target Alignment</title>\n")
        out.write("<style>\n")
        out.write("body { font-family: Menlo, Consolas, monospace; margin: 16px; }\n")
        out.write(".line { white-space: pre; margin: 1px 0; }\n")
        out.write(".ref { margin-bottom: 10px; font-weight: bold; }\n")
        out.write(".match { color: #1f2937; }\n")
        out.write(".sub { color: #b91c1c; font-weight: 700; }\n")
        out.write(".ins { color: #1d4ed8; font-weight: 700; }\n")
        out.write(".gap { color: #9ca3af; }\n")
        out.write(".blank { color: transparent; }\n")
        out.write(".cluster { border: 1px solid #e5e7eb; border-radius: 8px; padding: 10px; margin-bottom: 14px; }\n")
        out.write(".meta { color: #6b7280; margin-bottom: 6px; }\n")
        out.write("</style></head><body>\n")
        out.write(f"<h2>Reference Clusters ({len(clusters)})</h2>\n")

        for idx, cluster in enumerate(clusters, start=1):
            out.write('<div class="cluster">\n')
            out.write(f"<h3>Cluster {idx}</h3>\n")
            out.write(f'<div class="meta">targets: {len(cluster["members"])}</div>\n')
            out.write("<div>Cluster Reference (first target)</div>\n")
            out.write(f'<div class="line ref">{escape(cluster["reference"])}</div>\n')
            out.write("<div>Aligned Targets</div>\n")

            for t_line in cluster["projected_targets"]:
                out.write(f'<div class="line">{render_colored_line(cluster["render_reference"], t_line)}</div>\n')

            out.write("</div>\n")

        out.write("</body></html>\n")


def write_cluster_mapping_output(mapping_file, clusters):
    with open(mapping_file, "w") as out:
        out.write("target_sequence\treference_sequence\tcluster_id\tidentity\terrors\n")
        for idx, cluster in enumerate(clusters, start=1):
            for member in cluster["members"]:
                out.write(
                    f'{member["target"]}\t{cluster["reference"]}\t{idx}\t'
                    f'{member["identity"]:.6f}\t{member["errors"]}\n'
                )


def evaluate_alignment(reference, target, min_common_substr, min_identity, max_gap_runs, max_fragments, max_errors):
    lcs = longest_common_substring_len(reference, target)
    if lcs < min_common_substr:
        return {
            "ok": False,
            "identity": 0.0,
            "errors": 0,
            "aligned_ref": reference,
            "aligned_target": target,
        }

    aligned_ref, aligned_target, _ = needleman_wunsch(reference, target)
    ident = alignment_identity(aligned_ref, aligned_target)
    if ident < min_identity:
        return {
            "ok": False,
            "identity": ident,
            "errors": 0,
            "aligned_ref": aligned_ref,
            "aligned_target": aligned_target,
        }

    fragmented = is_fragmented_alignment(aligned_target, max_gap_runs, max_fragments)
    if fragmented:
        return {
            "ok": False,
            "identity": ident,
            "errors": 0,
            "aligned_ref": aligned_ref,
            "aligned_target": aligned_target,
        }

    err_cnt = alignment_error_count(aligned_ref, aligned_target)
    if err_cnt > max_errors:
        return {
            "ok": False,
            "identity": ident,
            "errors": err_cnt,
            "aligned_ref": aligned_ref,
            "aligned_target": aligned_target,
        }

    return {
        "ok": True,
        "identity": ident,
        "errors": err_cnt,
        "aligned_ref": aligned_ref,
        "aligned_target": aligned_target,
    }


def load_targets(file_path):
    targets = []
    with open(file_path, "r") as f:
        for line in f:
            seq = line.strip()
            if seq:
                targets.append(seq)
    return targets


def compare_one_file(
    in_file,
    min_common_substr,
    min_identity,
    max_gap_runs,
    max_fragments,
    max_errors,
    html_out_file,
    mapping_out_file,
):
    targets = load_targets(in_file)
    if not targets:
        return 0, 0, 0

    # Dynamic multi-reference clustering.
    clusters = []

    def new_cluster(seq):
        ref_len = len(seq)
        cluster = {
            "reference": seq,
            "ref_len": ref_len,
            "members": [{
                "target": seq,
                "errors": 0,
                "identity": 1.0,
                "aligned_ref": seq,
                "aligned_target": seq,
            }],
        }
        clusters.append(cluster)

    new_cluster(targets[0])

    aligned_count = 1
    fallback_count = 0

    for target in targets[1:]:
        candidates = []
        for c_idx, cluster in enumerate(clusters):
            ev = evaluate_alignment(
                cluster["reference"],
                target,
                min_common_substr,
                min_identity,
                max_gap_runs,
                max_fragments,
                max_errors,
            )
            if ev["ok"]:
                candidates.append((c_idx, ev))

        if not candidates:
            new_cluster(target)
            fallback_count += 1
            continue

        # Choose the most similar cluster: fewer errors, then higher identity.
        best_idx, best_ev = min(candidates, key=lambda x: (x[1]["errors"], -x[1]["identity"]))
        cluster = clusters[best_idx]

        cluster["members"].append({
            "target": target,
            "errors": best_ev["errors"],
            "identity": best_ev["identity"],
            "aligned_ref": best_ev["aligned_ref"],
            "aligned_target": best_ev["aligned_target"],
        })
        aligned_count += 1

    # Build projected lines per cluster and sort by similarity.
    for cluster in clusters:
        members = sorted(cluster["members"], key=lambda m: (m["errors"], -m["identity"], m["target"]))
        ref_insertions = choose_reference_insertions(members, cluster["ref_len"])
        cluster["render_reference"] = build_reference_with_gap_slots(cluster["reference"], ref_insertions)
        cluster["projected_targets"] = [
            build_projected_target(m["aligned_ref"], m["aligned_target"], ref_insertions, cluster["ref_len"])
            for m in members
        ]

    # Larger/similar groups first.
    clusters.sort(key=lambda c: len(c["members"]), reverse=True)
    write_clustered_html_output(html_out_file, clusters)
    write_cluster_mapping_output(mapping_out_file, clusters)

    return len(targets), aligned_count, fallback_count


def process_one_file(task):
    """Worker entry for multiprocessing."""
    (
        in_file,
        html_out_file,
        mapping_out_file,
        min_common_substr,
        min_identity,
        max_gap_runs,
        max_fragments,
        max_errors,
    ) = task

    n_targets, n_aligned, n_fallback = compare_one_file(
        in_file,
        min_common_substr,
        min_identity,
        max_gap_runs,
        max_fragments,
        max_errors,
        html_out_file,
        mapping_out_file,
    )

    return in_file.name, n_targets, n_aligned, n_fallback


def main():
    args = parse_args()

    input_dir = args.input_dir
    html_out_dir = args.html_out_dir
    if html_out_dir is not None:
        html_out_dir.mkdir(parents=True, exist_ok=True)
    else:
        html_out_dir = input_dir / "html"
        html_out_dir.mkdir(parents=True, exist_ok=True)

    in_files = sorted(input_dir.glob(args.pattern))
    if not in_files:
        raise FileNotFoundError(f"[ERROR] No files found in {input_dir} with pattern {args.pattern}")

    print(f"[INFO] Alignment backend: {ALIGNMENT_BACKEND}")

    tasks = [
        (
            in_file,
            html_out_dir / f"{in_file.stem}.html",
            html_out_dir / f"{in_file.stem}.cluster_mapping.tsv",
            args.min_common_substr,
            args.min_identity,
            args.max_gap_runs,
            args.max_fragments,
            args.max_errors,
        )
        for in_file in in_files
    ]

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_one_file, task) for task in tasks]
        total = len(futures)
        finished = 0
        for fut in as_completed(futures):
            fut.result()
            finished += 1
            print(f"\r[INFO] Processed {finished}/{total} files", end="", flush=True)

    print()

if __name__ == "__main__":
    main()
