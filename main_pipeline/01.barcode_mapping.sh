#!/usr/bin/env bash
#SBATCH -J barcode_mapping
#SBATCH -c 1
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/python/01.generate_cellbarcode_mapping.py"

CODE_LEN="16"
POOL_SIZE="800"
SEED="42"
PIXI_ENV="main"

usage() {
    cat <<'EOF'
Usage:
    barcode_mapping.sh --input FILE [--output FILE] [options]

Required arguments:
    --input  FILE        Input CSV/TSV with combo name, combo_sequence, well
    --output FILE        Output CSV path

Optional arguments:
    --code-len INT      Length of mapped barcode (default: 16)
    --pool-size INT     Candidate pool size (default: 800)
    --seed INT          Random seed (default: 42)
    --env STR           Pixi environment name (default: main)

Other:
    -h, --help          Show this help message

Examples:
    bash barcode_mapping.sh --input 4CL.tsv --output barcode_mapping.csv
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --code-len)
            CODE_LEN="$2"
            shift 2
            ;;
        --pool-size)
            POOL_SIZE="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
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

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "[ERROR] --input and --output are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "[ERROR] Input file not found: $INPUT" >&2
    exit 1
fi

if [[ ! -f "$PY_SCRIPT" ]]; then
    echo "[ERROR] Python script not found: $PY_SCRIPT" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

echo "[INFO] input=$INPUT"
echo "[INFO] output=$OUTPUT"
pixi run -e "$PIXI_ENV" python "$PY_SCRIPT" \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --code-len "$CODE_LEN" \
    --pool-size "$POOL_SIZE" \
    --seed "$SEED"
