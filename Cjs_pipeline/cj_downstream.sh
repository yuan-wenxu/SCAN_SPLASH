#!/usr/bin/env bash
#SBATCH -J cj_downstream
#SBATCH -c 1
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_PY="$SCRIPT_DIR/python/build_cell_anchor_cj_matrix.py"
FILTER_PY="$SCRIPT_DIR/python/filter_cell_anchor_cj_matrix.py"

CJS_TSV_GZ=""
BARCODE_MAPPING_CSV=""
MATRIX_DIR=""
WHITELIST=""

MIN_ANCHORS_PER_CELL="2000"
MIN_CELLS_PER_ANCHOR="3"
DROP_UNMAPPED="false"
SCATTER_PNG=""
PIXI_ENV="main"

usage() {
    cat <<'EOF'
Usage:
    cj_downstream.sh --cjs-tsv-gz FILE --barcode-mapping-csv FILE --matrix-dir DIR --whitelist FILE [options]

Required arguments:
    --cjs-tsv-gz FILE         Path to SPLASH result_Cjs/cjs.tsv.gz
    --barcode-mapping-csv FILE
                              Path to barcode_mapping.csv
    --matrix-dir DIR          Output directory for cell x anchor Cj matrix
    --whitelist FILE          Whitelist TSV for filtering

Optional arguments:
    --min-anchors-per-cell INT  Keep cells with >= this many nonzero anchors (default: 2000)
    --min-cells-per-anchor INT  Keep anchors in >= this many kept cells (default: 3)
    --drop-unmapped             Drop barcodes not found in mapping table
    --scatter-png FILE          Custom scatter PNG output path
    --env STR                   Pixi environment name (default: main)
    -h, --help                  Show this help message

Examples:
    bash cj_downstream.sh \
        --cjs-tsv-gz /path/to/result_Cjs/cjs.tsv.gz \
        --barcode-mapping-csv /path/to/barcode_mapping.csv \
        --matrix-dir /path/to/result_Cjs/matrix_from_cj \
        --whitelist /path/to/whitelist.tsv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cjs-tsv-gz)
            CJS_TSV_GZ="$2"
            shift 2
            ;;
        --barcode-mapping-csv)
            BARCODE_MAPPING_CSV="$2"
            shift 2
            ;;
        --matrix-dir)
            MATRIX_DIR="$2"
            shift 2
            ;;
        --whitelist)
            WHITELIST="$2"
            shift 2
            ;;
        --min-anchors-per-cell)
            MIN_ANCHORS_PER_CELL="$2"
            shift 2
            ;;
        --min-cells-per-anchor)
            MIN_CELLS_PER_ANCHOR="$2"
            shift 2
            ;;
        --drop-unmapped)
            DROP_UNMAPPED="true"
            shift
            ;;
        --scatter-png)
            SCATTER_PNG="$2"
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

if [[ -z "$CJS_TSV_GZ" || -z "$BARCODE_MAPPING_CSV" || -z "$MATRIX_DIR" || -z "$WHITELIST" ]]; then
    echo "[ERROR] --cjs-tsv-gz, --barcode-mapping-csv, --matrix-dir, and --whitelist are required." >&2
    usage >&2
    exit 1
fi

for _path in "$BUILD_PY" "$FILTER_PY" "$CJS_TSV_GZ" "$BARCODE_MAPPING_CSV" "$WHITELIST"; do
    if [[ ! -f "$_path" ]]; then
        echo "[ERROR] File not found: $_path" >&2
        exit 1
    fi
done

mkdir -p "$MATRIX_DIR"

build_cmd=(
    pixi run -e "$PIXI_ENV" python "$BUILD_PY"
    --cjs-tsv-gz "$CJS_TSV_GZ"
    --barcode-mapping-csv "$BARCODE_MAPPING_CSV"
    --out-dir "$MATRIX_DIR"
)
if [[ "$DROP_UNMAPPED" == "true" ]]; then
    build_cmd+=(--drop-unmapped)
fi

echo "[INFO] Step 1: Build cell x anchor Cj matrix"
"${build_cmd[@]}"

echo "[INFO] Step 2: Filter Cj matrix"
filter_cmd=(
    pixi run -e "$PIXI_ENV" python "$FILTER_PY"
    --matrix-dir "$MATRIX_DIR"
    --whitelist "$WHITELIST"
    --min-anchors-per-cell "$MIN_ANCHORS_PER_CELL"
    --min-cells-per-anchor "$MIN_CELLS_PER_ANCHOR"
)
if [[ -n "$SCATTER_PNG" ]]; then
    filter_cmd+=(--scatter-png "$SCATTER_PNG")
fi
"${filter_cmd[@]}"

echo "[OK] Cj downstream completed"
echo "[OK] Matrix dir: $MATRIX_DIR"
echo "[OK] Filtered dir: $MATRIX_DIR/filtered"
