#!/usr/bin/env python3
"""
Extract wanted anchors from assigned_reads.txt and get their corresponding targets from features.tsv.
"""
import argparse
from collections import defaultdict
from pathlib import Path


COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq):
    return seq.translate(COMP)[::-1]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract wanted anchors from assigned_reads.txt and filter targets from features.tsv"
    )
    parser.add_argument(
        "--assigned-reads",
        type=Path,
        required=True,
        help="Path to assigned_reads.txt file (SAM format)"
    )
    parser.add_argument(
        "--features-file",
        type=Path,
        required=True,
        help="Path to features.tsv file with anchor-target mapping"
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for organized anchor-target files"
    )
    parser.add_argument(
        "--merge-distance",
        type=int,
        default=31,
        help="Merge anchors when POS difference is within this value (default: 31); no transitive/connected merge"
    )
    return parser.parse_args()


def extract_wanted_anchors(assigned_reads_file):
    """
    Extract unique anchor sequences from assigned_reads.txt (column 10).
    Returns: set of anchor sequences
    """
    wanted_anchors = {}
    
    with open(assigned_reads_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # SAM format: column 10 contains the sequence
            fields = line.split('\t')
            
            if len(fields) >= 10:
                anchor_id = fields[0].split("|")[0].strip()  # Column 1 (0-indexed: 0) is the read ID
                rname = fields[2].strip() if len(fields) > 2 else "*"
                try:
                    pos = int(fields[3])
                except ValueError:
                    pos = 0
                try:
                    flag = int(fields[1])
                except ValueError:
                    flag = 0
                strand = "-" if (flag & 16) else "+"
                anchor_seq = fields[9].strip()  # Column 10 (0-indexed: 9)
                if anchor_seq and len(anchor_seq) > 0:
                    wanted_anchors[anchor_id] = {
                        "seq": anchor_seq,
                        "rname": rname,
                        "pos": pos,
                        "strand": strand,
                    }
    
    return wanted_anchors


def filter_targets_by_wanted_anchors(features_file, wanted_anchors):
    """
    Parse features.tsv and extract targets for wanted anchors.
    Supports exact and reverse-complement anchor matching.
    Returns:
        anchor_targets: dict of wanted_anchor -> set of targets
        exact_anchor_hits: number of wanted anchors matched exactly
        rc_anchor_hits: number of wanted anchors matched by reverse complement
    """
    anchor_targets = defaultdict(set)
    exact_seen = set()
    rc_seen = set()
    wanted_anchor_seqs = {meta["seq"] for meta in wanted_anchors.values()}
    
    with open(features_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # TSV format: columns are anchor_target, anchor, target
            fields = line.split('\t')
            
            if len(fields) < 3:
                continue
            
            # Extract anchor and target
            anchor = fields[1].strip()
            target = fields[2].strip()
            
            if not target:
                continue

            # Exact orientation match
            if anchor in wanted_anchor_seqs:
                anchor_targets[anchor].add(target)
                exact_seen.add(anchor)
                continue

            # Reverse-complement orientation match
            anchor_rc = reverse_complement(anchor)
            if anchor_rc in wanted_anchor_seqs:
                anchor_targets[anchor_rc].add(target)
                rc_seen.add(anchor_rc)
    
    return anchor_targets, len(exact_seen), len(rc_seen)


def merge_anchor_ids_by_position(wanted_anchors, anchor_targets, merge_distance):
    """
    Merge anchors by position on SAM POS (column 4) within same (rname, strand).
    Use non-transitive grouping: anchors are merged only if they are within merge_distance
    from the first anchor in the group (no connected-component chaining).
    Keep only one anchor_id per component: the one with most targets.
    """
    groups = defaultdict(list)
    for anchor_id, meta in wanted_anchors.items():
        key = (meta["rname"], meta["strand"])
        groups[key].append((anchor_id, meta))

    kept_ids = set()
    merged_components = 0
    removed_ids = 0

    for _, items in groups.items():
        # Unmapped/invalid positions are treated as standalone anchors.
        valid_items = []
        for anchor_id, meta in items:
            if meta["rname"] == "*" or meta["pos"] <= 0:
                kept_ids.add(anchor_id)
            else:
                valid_items.append((anchor_id, meta))

        valid_items.sort(key=lambda x: x[1]["pos"])
        i = 0
        while i < len(valid_items):
            comp = [valid_items[i]]
            base_pos = valid_items[i][1]["pos"]
            j = i + 1
            while j < len(valid_items):
                if valid_items[j][1]["pos"] - base_pos <= merge_distance:
                    comp.append(valid_items[j])
                    j += 1
                else:
                    break

            if len(comp) == 1:
                kept_ids.add(comp[0][0])
            else:
                merged_components += 1

                def score(item):
                    anchor_id, meta = item
                    return len(anchor_targets.get(meta["seq"], set()))

                best = max(comp, key=lambda it: (score(it), it[0]))
                kept_ids.add(best[0])
                removed_ids += len(comp) - 1

            i = j

    return kept_ids, merged_components, removed_ids


def write_targets_by_anchor(anchor_targets, output_dir, wanted_anchors, kept_anchor_ids):
    """
    Write each anchor's targets to a separate file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    kept_anchor_ids = set(kept_anchor_ids)

    # Remove stale output files from previous runs.
    valid_names = {f"{anchor_id}.txt" for anchor_id in kept_anchor_ids}
    for old_file in output_dir.glob("*.txt"):
        if old_file.name not in valid_names:
            old_file.unlink()
    
    file_count = 0
    total_records = 0
    
    for anchor_id in sorted(kept_anchor_ids):
        anchor_seq = wanted_anchors[anchor_id]["seq"]
        targets = sorted(anchor_targets.get(anchor_seq, set()))
        output_file = output_dir / f"{anchor_id}.txt"
        
        with open(output_file, 'w') as f:
            for target in targets:
                f.write(f"{target}\n")
        
        file_count += 1
        total_records += len(targets)
        #print(f"  {anchor_id}: {len(targets)} target(s)")
    
    return file_count, total_records


def main():
    args = parse_args()

    assigned_reads_file = args.assigned_reads
    if not assigned_reads_file.exists():
        raise FileNotFoundError(f"[ERROR] assigned_reads.txt not found: {assigned_reads_file}")

    features_file = args.features_file
    if not features_file.exists():
        raise FileNotFoundError(f"[ERROR] Features file not found: {features_file}")

    output_dir = args.out_dir
    merge_distance = args.merge_distance

    print(f"[INFO] Extracting wanted anchors from: {assigned_reads_file}")
    wanted_anchors = extract_wanted_anchors(assigned_reads_file)

    print(f"[INFO] Filtering targets from features file: {features_file}")
    anchor_targets, exact_hits, rc_hits = filter_targets_by_wanted_anchors(features_file, wanted_anchors)
    print(f"[INFO] Total wanted anchors: {len(wanted_anchors)}")
    print(f"[INFO] Exact orientation matches: {exact_hits}")
    print(f"[INFO] Reverse-complement orientation matches: {rc_hits}")

    print(f"[INFO] Merging anchors by position with merge_distance={merge_distance}")
    kept_anchor_ids, merged_components, removed_ids = merge_anchor_ids_by_position(
        wanted_anchors,
        anchor_targets,
        merge_distance,
    )
    print(f"[INFO] Merged {merged_components} components, removed {removed_ids} redundant anchors")
    
    print(f"[INFO] Writing targets for {len(kept_anchor_ids)} kept anchors to directory: {output_dir}")
    file_count, total_records = write_targets_by_anchor(
        anchor_targets,
        output_dir,
        wanted_anchors,
        kept_anchor_ids,
    )
    print(f"[INFO] Total files written: {file_count}")
    print(f"[INFO] Total anchor-target records: {total_records}")

if __name__ == "__main__":
    main()
