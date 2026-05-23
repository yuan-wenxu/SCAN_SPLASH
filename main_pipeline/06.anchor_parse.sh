#!/usr/bin/env bash
#SBATCH -J anchor_pipeline
#SBATCH -c 8
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARSE_PY="$SCRIPT_DIR/python/06.parse_assigned_reads.py"
COMPARE_PY="$SCRIPT_DIR/python/06.compare_targets_indel.py"

# Required
ASSIGNED_READS=""
FEATURES_FILE=""
OUT_DIR=""

# Optional: parse_assigned_reads.py
MERGE_DISTANCE="15"

# Optional: compare_targets_indel.py
HTML_OUT_DIR=""
PATTERN="*.txt"
MIN_COMMON_SUBSTR="5"
MIN_IDENTITY="0.5"
MAX_GAP_RUNS="2"
MAX_FRAGMENTS="3"
MAX_ERRORS="6"
WORKERS="8"

PIXI_ENV="main"

usage() {
    cat <<'EOF'
Usage:
    anchor.sh --assigned-reads FILE --features-file FILE --out-dir DIR [options]

Required arguments:
    --assigned-reads FILE   Path to assigned_reads.txt (SAM format)
    --features-file FILE    Path to features.tsv with anchor-target mapping
    --out-dir DIR           Output directory for anchor-target files and HTML

Optional arguments (parse_assigned_reads):
    --merge-distance INT    Merge distance for position-based anchor dedup (default: 15)

Optional arguments (compare_targets_indel):
    --html-out-dir DIR      Output directory for HTML files (default: OUT_DIR/html)
    --pattern STR           Glob pattern for input files (default: *.txt)
    --min-common-substr INT Minimum LCS length to consider alignable (default: 5)
    --min-identity FLOAT    Minimum alignment identity (default: 0.5)
    --max-gap-runs INT      Maximum gap runs in alignment (default: 2)
    --max-fragments INT     Maximum nucleotide fragments in alignment (default: 3)
    --max-errors INT        Maximum allowed errors per target (default: 6)
    --workers INT           Number of parallel workers (default: 8)

General:
    --env STR               Pixi environment name (default: main)
    -h, --help              Show this help message

Examples:
    bash anchor.sh \
        --assigned-reads /path/to/assigned_reads.txt \
        --features-file /path/to/features.tsv \
        --out-dir /path/to/output
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --assigned-reads)
            ASSIGNED_READS="$2"
            shift 2
            ;;
        --features-file)
            FEATURES_FILE="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --merge-distance)
            MERGE_DISTANCE="$2"
            shift 2
            ;;
        --html-out-dir)
            HTML_OUT_DIR="$2"
            shift 2
            ;;
        --pattern)
            PATTERN="$2"
            shift 2
            ;;
        --min-common-substr)
            MIN_COMMON_SUBSTR="$2"
            shift 2
            ;;
        --min-identity)
            MIN_IDENTITY="$2"
            shift 2
            ;;
        --max-gap-runs)
            MAX_GAP_RUNS="$2"
            shift 2
            ;;
        --max-fragments)
            MAX_FRAGMENTS="$2"
            shift 2
            ;;
        --max-errors)
            MAX_ERRORS="$2"
            shift 2
            ;;
        --workers)
            WORKERS="$2"
            shift 2
            ;;
        --env)
            PIXI_ENV="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$ASSIGNED_READS" || -z "$FEATURES_FILE" || -z "$OUT_DIR" ]]; then
    echo "[ERROR] --assigned-reads, --features-file, and --out-dir are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$ASSIGNED_READS" ]]; then
    echo "[ERROR] assigned-reads file not found: $ASSIGNED_READS" >&2
    exit 1
fi

if [[ ! -f "$FEATURES_FILE" ]]; then
    echo "[ERROR] features-file not found: $FEATURES_FILE" >&2
    exit 1
fi

if [[ ! -f "$PARSE_PY" ]]; then
    echo "[ERROR] Python script not found: $PARSE_PY" >&2
    exit 1
fi

if [[ ! -f "$COMPARE_PY" ]]; then
    echo "[ERROR] Python script not found: $COMPARE_PY" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

if [[ -z "$HTML_OUT_DIR" ]]; then
    HTML_OUT_DIR="$OUT_DIR/html"
fi

echo "[INFO] Step 1: Parsing assigned reads and extracting anchor targets"
pixi run -e "$PIXI_ENV" python "$PARSE_PY" \
    --assigned-reads "$ASSIGNED_READS" \
    --features-file "$FEATURES_FILE" \
    --out-dir "$OUT_DIR" \
    --merge-distance "$MERGE_DISTANCE"

echo "[INFO] Step 2: Comparing targets and generating HTML visualizations"
pixi run -e "$PIXI_ENV" python "$COMPARE_PY" \
    --input-dir "$OUT_DIR" \
    --pattern "$PATTERN" \
    --html-out-dir "$HTML_OUT_DIR" \
    --min-common-substr "$MIN_COMMON_SUBSTR" \
    --min-identity "$MIN_IDENTITY" \
    --max-gap-runs "$MAX_GAP_RUNS" \
    --max-fragments "$MAX_FRAGMENTS" \
    --max-errors "$MAX_ERRORS" \
    --workers "$WORKERS"

echo "[OK] Anchor pipeline completed. HTML output: $HTML_OUT_DIR"
