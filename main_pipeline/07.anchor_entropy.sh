#!/usr/bin/env bash
#SBATCH -J remap_anchor_entropy
#SBATCH -c 4
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMAP_PY="$SCRIPT_DIR/python/07.remap_filtered_matrix_targets.py"
ENTROPY_PY="$SCRIPT_DIR/python/07.compute_anchor_entropy.py"
ANNOTATE_PY="$SCRIPT_DIR/python/07.annotate_anchor_entropy_gene_name.py"

MATRIX_DIR=""
MAPPING_DIR=""
GTF=""

REMAPPED_OUT_DIR=""
ANCHOR_ALIGN_DIR=""
MIN_TOTAL_COUNT="1"
PIXI_ENV="main"

FEATURES_FILE="whitelist_filtered.features.tsv"
CELLS_FILE="whitelist_filtered.cells.tsv"
MATRIX_FILE="whitelist_filtered.matrix.mtx.gz"
ANCHORS_FASTQ="anchor_align/anchors.fastq.gz"

ENTROPY_CSV=""
ANNOTATED_CSV=""
HISTOGRAM_PNG=""

usage() {
    cat <<'EOF'
Usage:
    anchor_entropy_similarity.sh --matrix-dir DIR --mapping-dir DIR --gtf FILE [options]

Required arguments:
    --matrix-dir DIR        Directory containing whitelist_filtered.{features,cells,matrix}
    --mapping-dir DIR       Directory containing anchor_*.cluster_mapping.tsv
    --gtf FILE              GTF file containing gene_id and gene_name

Optional arguments:
    --remapped-out-dir DIR  Output directory for remapped matrix files
                            (default: <matrix-dir>/remapped)
    --features-file FILE    Features filename in --matrix-dir (default: whitelist_filtered.features.tsv)
    --cells-file FILE       Cells filename in --matrix-dir (default: whitelist_filtered.cells.tsv)
    --matrix-file FILE      Matrix filename in --matrix-dir (default: whitelist_filtered.matrix.mtx.gz)
    --anchors-fastq FILE    Anchor FASTQ path relative to --matrix-dir
                            (default: anchor_align/anchors.fastq.gz)
    --anchor-align-dir DIR  Directory containing anchors.fastq.gz and assigned_reads.tsv
                            (default: <matrix-dir>/anchor_align)
    --min-total-count INT   Min count filter for entropy output (default: 1)
    --entropy-csv FILE      Output path for anchor entropy CSV
                            (default: <remapped-out-dir>/anchor_entropy.csv)
    --annotated-csv FILE    Output path for gene_name-annotated entropy CSV
                            (default: <remapped-out-dir>/anchor_entropy.with_gene_name.csv)
    --histogram FILE        Output path for entropy histogram PNG
                            (default: <remapped-out-dir>/anchor_entropy_histogram.png)
    --env STR               Pixi environment name (default: main)
    -h, --help              Show this help message

Examples:
    bash anchor_entropy_similarity.sh \
        --matrix-dir /path/to/filtered \
        --mapping-dir /path/to/anchors_targets_wanted_html \
        --gtf /mnt/ref/genes.gtf
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --matrix-dir)
            MATRIX_DIR="$2"
            shift 2
            ;;
        --mapping-dir)
            MAPPING_DIR="$2"
            shift 2
            ;;
        --gtf)
            GTF="$2"
            shift 2
            ;;
        --remapped-out-dir)
            REMAPPED_OUT_DIR="$2"
            shift 2
            ;;
        --features-file)
            FEATURES_FILE="$2"
            shift 2
            ;;
        --cells-file)
            CELLS_FILE="$2"
            shift 2
            ;;
        --matrix-file)
            MATRIX_FILE="$2"
            shift 2
            ;;
        --anchors-fastq)
            ANCHORS_FASTQ="$2"
            shift 2
            ;;
        --anchor-align-dir)
            ANCHOR_ALIGN_DIR="$2"
            shift 2
            ;;
        --min-total-count)
            MIN_TOTAL_COUNT="$2"
            shift 2
            ;;
        --entropy-csv)
            ENTROPY_CSV="$2"
            shift 2
            ;;
        --annotated-csv)
            ANNOTATED_CSV="$2"
            shift 2
            ;;
        --histogram)
            HISTOGRAM_PNG="$2"
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

if [[ -z "$MATRIX_DIR" || -z "$MAPPING_DIR" || -z "$GTF" ]]; then
    echo "[ERROR] --matrix-dir, --mapping-dir, and --gtf are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "[ERROR] matrix-dir not found: $MATRIX_DIR" >&2
    exit 1
fi

if [[ ! -d "$MAPPING_DIR" ]]; then
    echo "[ERROR] mapping-dir not found: $MAPPING_DIR" >&2
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "[ERROR] gtf not found: $GTF" >&2
    exit 1
fi

for _py in "$REMAP_PY" "$ENTROPY_PY" "$ANNOTATE_PY"; do
    if [[ ! -f "$_py" ]]; then
        echo "[ERROR] Python script not found: $_py" >&2
        exit 1
    fi
done

if [[ -z "$REMAPPED_OUT_DIR" ]]; then
    REMAPPED_OUT_DIR="$MATRIX_DIR/remapped"
fi

if [[ -z "$ANCHOR_ALIGN_DIR" ]]; then
    ANCHOR_ALIGN_DIR="$MATRIX_DIR/anchor_align"
fi

if [[ -z "$ENTROPY_CSV" ]]; then
    ENTROPY_CSV="$REMAPPED_OUT_DIR/anchor_entropy.csv"
fi

if [[ -z "$ANNOTATED_CSV" ]]; then
    ANNOTATED_CSV="$REMAPPED_OUT_DIR/anchor_entropy.with_gene_name.csv"
fi

if [[ -z "$HISTOGRAM_PNG" ]]; then
    HISTOGRAM_PNG="$REMAPPED_OUT_DIR/anchor_entropy_histogram.png"
fi

echo "[INFO] Step 1: Remap filtered matrix targets"
pixi run -e "$PIXI_ENV" python3 "$REMAP_PY" \
    --matrix-dir "$MATRIX_DIR" \
    --mapping-dir "$MAPPING_DIR" \
    --output-dir "$REMAPPED_OUT_DIR" \
    --features-file "$FEATURES_FILE" \
    --cells-file "$CELLS_FILE" \
    --matrix-file "$MATRIX_FILE" \
    --anchors-fastq "$ANCHORS_FASTQ"

echo "[INFO] Step 2: Compute anchor entropy on remapped matrix"
entropy_cmd=(
    pixi run -e "$PIXI_ENV" python3 "$ENTROPY_PY"
    --matrix-dir "$REMAPPED_OUT_DIR"
    --output "$ENTROPY_CSV"
    --histogram "$HISTOGRAM_PNG"
    --min-total-count "$MIN_TOTAL_COUNT"
)
if [[ -n "$ANCHOR_ALIGN_DIR" ]]; then
    entropy_cmd+=(--anchor-align-dir "$ANCHOR_ALIGN_DIR")
fi
"${entropy_cmd[@]}"

echo "[INFO] Step 3: Annotate entropy CSV with gene_name"
pixi run -e "$PIXI_ENV" python3 "$ANNOTATE_PY" \
    --input "$ENTROPY_CSV" \
    --gtf "$GTF" \
    --output "$ANNOTATED_CSV"

echo "[OK] Remap + anchor entropy pipeline completed"
echo "[OK] Remapped matrix:    $REMAPPED_OUT_DIR"
echo "[OK] Entropy CSV:       $ENTROPY_CSV"
echo "[OK] Annotated CSV:     $ANNOTATED_CSV"
echo "[OK] Entropy histogram: $HISTOGRAM_PNG"