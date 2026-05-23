#!/usr/bin/env bash
#SBATCH -J build_filter_matrix
#SBATCH -c 1
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_PY="$SCRIPT_DIR/python/04.build_cell_anchor_target_matrix.py"
FILTER_PY="$SCRIPT_DIR/python/04.filter_cell_anchor_target_matrix.py"
ANCHOR_FILTER_PY="$SCRIPT_DIR/python/04.filter_anchor_by_cell_pct.py"

OUT_DIR=""
SATC_DUMP_BIN=""
WHITELIST=""

CBC_MAPPING_CSV=""
ANCHOR_LEN="31"
TARGET_LEN="31"
KEEP_DUMP_TEXT="false"
MIN_TOTAL_COUNT="20000"
MIN_FEATURES="10000"
MIN_ANCHOR_CELL_PCT="10"
PIXI_ENV="main"

usage() {
    cat <<'EOF'
Usage:
    matrix.sh --out-dir DIR --satc-dump-bin FILE --whitelist FILE [options]

Required arguments:
    --out-dir DIR          Output directory used by the matrix scripts
    --satc-dump-bin FILE   Path to satc_dump binary
    --whitelist FILE       Whitelist TSV for filtering

Optional arguments:
    --cbc-mapping-csv FILE   Optional CBC reverse-mapping CSV
    --anchor-len INT         Anchor length for SATC decoding (default: 31)
    --target-len INT         Target length for SATC decoding (default: 31)
    --keep-dump-text         Keep intermediate SATC text dump
    --min-total-count INT    Cell QC threshold (default: 20000)
    --min-features INT       Feature QC threshold (default: 10000)
    --min-anchor-cell-pct F  Anchor keep threshold in percent (default: 10)
    --env STR                Pixi environment name (default: main)
    -h, --help               Show this help message

Examples:
    bash matrix.sh \
        --out-dir /path/to/matrix \
        --satc-dump-bin /mnt/splicing/splash/satc_dump \
        --whitelist /path/to/whitelist.tsv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --satc-dump-bin)
            SATC_DUMP_BIN="$2"
            shift 2
            ;;
        --whitelist)
            WHITELIST="$2"
            shift 2
            ;;
        --cbc-mapping-csv)
            CBC_MAPPING_CSV="$2"
            shift 2
            ;;
        --anchor-len)
            ANCHOR_LEN="$2"
            shift 2
            ;;
        --target-len)
            TARGET_LEN="$2"
            shift 2
            ;;
        --keep-dump-text)
            KEEP_DUMP_TEXT="true"
            shift
            ;;
        --min-total-count)
            MIN_TOTAL_COUNT="$2"
            shift 2
            ;;
        --min-features)
            MIN_FEATURES="$2"
            shift 2
            ;;
        --min-anchor-cell-pct)
            MIN_ANCHOR_CELL_PCT="$2"
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

if [[ -z "$OUT_DIR" || -z "$SATC_DUMP_BIN" || -z "$WHITELIST" ]]; then
    echo "[ERROR] --out-dir, --satc-dump-bin, and --whitelist are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$BUILD_PY" ]]; then
    echo "[ERROR] Python script not found: $BUILD_PY" >&2
    exit 1
fi

if [[ ! -f "$FILTER_PY" ]]; then
    echo "[ERROR] Python script not found: $FILTER_PY" >&2
    exit 1
fi

if [[ ! -f "$ANCHOR_FILTER_PY" ]]; then
    echo "[ERROR] Python script not found: $ANCHOR_FILTER_PY" >&2
    exit 1
fi

if [[ ! -f "$SATC_DUMP_BIN" ]]; then
    echo "[ERROR] satc_dump binary not found: $SATC_DUMP_BIN" >&2
    exit 1
fi

if [[ ! -f "$WHITELIST" ]]; then
    echo "[ERROR] whitelist file not found: $WHITELIST" >&2
    exit 1
fi

if [[ -n "$CBC_MAPPING_CSV" && ! -f "$CBC_MAPPING_CSV" ]]; then
    echo "[ERROR] CBC mapping CSV not found: $CBC_MAPPING_CSV" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

build_cmd=(
    pixi run -e "$PIXI_ENV" python "$BUILD_PY"
    --out-dir "$OUT_DIR"
    --satc-dump-bin "$SATC_DUMP_BIN"
    --anchor-len "$ANCHOR_LEN"
    --target-len "$TARGET_LEN"
)

if [[ -n "$CBC_MAPPING_CSV" ]]; then
    build_cmd+=(--cbc-mapping-csv "$CBC_MAPPING_CSV")
fi

if [[ "$KEEP_DUMP_TEXT" == "true" ]]; then
    build_cmd+=(--keep-dump-text)
fi

echo "[INFO] Building cell x anchor-target matrix"
"${build_cmd[@]}"

echo "[INFO] Filtering matrix by whitelist and QC"
pixi run -e "$PIXI_ENV" python "$FILTER_PY" \
    --out-dir "$OUT_DIR/matrix" \
    --whitelist "$WHITELIST" \
    --min-total-count "$MIN_TOTAL_COUNT" \
    --min-features "$MIN_FEATURES"

echo "[INFO] Filtering anchors by detected-cell percentage"
pixi run -e "$PIXI_ENV" python "$ANCHOR_FILTER_PY" \
    --matrix-dir "$OUT_DIR/matrix/filtered" \
    --input-prefix "whitelist_filtered" \
    --output-prefix "anchor_pct_filtered" \
    --min-anchor-cell-pct "$MIN_ANCHOR_CELL_PCT"

echo "[OK] Matrix build and filter completed"
