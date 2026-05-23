#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="$SCRIPT_DIR/config.sh"
MODE="local"
STEPS_TO_RUN="1,2,3,4,5,6,7"
TOTAL_STEPS="7"

usage() {
    cat <<'EOF'
Usage:
    main.sh [--config FILE] [--mode MODE] [--steps LIST]

Options:
    --config FILE     Config file path (default: config.sh next to this script)
    --mode MODE       Execution mode: local or hpc (default: local)
                        local  - run steps sequentially in current shell
                        hpc    - submit each step as a SLURM job with dependency chaining
    --steps LIST      Comma-separated step numbers to run (default: 1,2,3,4,5,6,7)
    -h, --help        Show this help message

Pipeline steps:
    1) barcode_mapping      Generate cell barcode combo -> sequence mapping
    2) batch_extract_cb_ub  Extract CB/UB from R2 to produce pseudo-R1 FASTQ
    3) run_splash           Run SPLASH on pseudo-R1 + R2 pairs
    4) matrix               Build and filter cell x anchor-target count matrix
    5) anchor_align_genome  Align anchors to genome and assign reads via featureCounts
    6) anchor               Parse assigned reads and generate HTML target visualizations
    7) anchor_entropy       Remap matrix, compute anchor entropy, and annotate gene names

Examples:
    # Run full pipeline locally
    bash main.sh --config my_config.sh --mode local

    # Submit all steps to SLURM with dependency chaining
    bash main.sh --config my_config.sh --mode hpc
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --mode)
            MODE="$2"
            shift 2
            ;;
        --steps)
            STEPS_TO_RUN="$2"
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

if [[ "$MODE" != "local" && "$MODE" != "hpc" ]]; then
    echo "[ERROR] --mode must be 'local' or 'hpc', got: $MODE" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "$CONFIG_PATH"

for _var in BARCODE_INPUT R2_DIR SPLASH_DIR GENOME_FA GENOME_GTF; do
    if [[ -z "${!_var:-}" ]]; then
        echo "[ERROR] Missing required config variable: $_var" >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------
mkdir -p "$PIPELINE_OUT_DIR"
MAPPING_CSV="$PIPELINE_OUT_DIR/barcode_mapping.csv"
R1_DIR="$PIPELINE_OUT_DIR/r1"
SPLASH_OUT_DIR="$PIPELINE_OUT_DIR/splash_results"
MATRIX_OUT_DIR="$SPLASH_OUT_DIR/matrix"
MATRIX_FILTERED_OUT_DIR="$MATRIX_OUT_DIR/filtered"
ANCHOR_ALIGN_OUT_DIR="$MATRIX_FILTERED_OUT_DIR/anchor_align"
PARSED_OUT_DIR="$MATRIX_FILTERED_OUT_DIR/anchors_targets_wanted"
HTML_OUT_DIR="$MATRIX_FILTERED_OUT_DIR/anchors_targets_wanted_html"
REMAPPED_OUT_DIR="$MATRIX_FILTERED_OUT_DIR/remapped"
SATC_DUMP_BIN="$SPLASH_DIR/satc_dump"
MATRIX_FILTERED_FEATURES="$MATRIX_FILTERED_OUT_DIR/whitelist_filtered.features.tsv"
MATRIX_FILTERED_COUNTS_FEATURES="$MATRIX_FILTERED_OUT_DIR/anchor_pct_filtered.features.tsv"
ASSIGNED_READS="$ANCHOR_ALIGN_OUT_DIR/assigned_reads.tsv"
LOG_DIR="$PIPELINE_OUT_DIR/logs"

mkdir -p "$R1_DIR" "$SPLASH_OUT_DIR" "$MATRIX_OUT_DIR" \
         "$ANCHOR_ALIGN_OUT_DIR" "$PARSED_OUT_DIR" "$HTML_OUT_DIR" "$REMAPPED_OUT_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
step_enabled() {
    echo ",$STEPS_TO_RUN," | grep -q ",$1,"
}

PREV_JID=""

# run_or_submit <step_num> <step_label> <script_path> [script_args...]
run_or_submit() {
    local step_num="$1"
    local step_label="$2"
    local script="$3"
    shift 3

    if ! step_enabled "$step_num"; then
        echo "[SKIP] Step $step_num: $step_label"
        return
    fi

    if [[ "$MODE" == "local" ]]; then
        echo "[STEP $step_num/$TOTAL_STEPS] $step_label"
        bash "$script" "$@"
        echo "[DONE $step_num/$TOTAL_STEPS] $step_label"
    else
        local sbatch_opts=(
            --parsable
            --job-name="scan_splash_s${step_num}_${SAMPLE_NAME}"
            --output="$LOG_DIR/step${step_num}_${step_label}_%j.out"
            --error="$LOG_DIR/step${step_num}_${step_label}_%j.err"
        )
        if [[ -n "$PREV_JID" ]]; then
            sbatch_opts+=(--dependency="afterok:$PREV_JID")
        fi
        local jid
        jid=$(sbatch "${sbatch_opts[@]}" "$script" "$@")
        echo "[SUBMITTED] Step $step_num: $step_label → job ID $jid"
        PREV_JID="$jid"
    fi
}

# ---------------------------------------------------------------------------
# Step 1: barcode_mapping
# ---------------------------------------------------------------------------
run_or_submit 1 "barcode_mapping" "$SCRIPT_DIR/01.barcode_mapping.sh" \
    --input  "$BARCODE_INPUT" \
    --output "$MAPPING_CSV" \
    --code-len "$MAPPING_CODE_LEN" \
    --pool-size "$MAPPING_POOL_SIZE" \
    --seed "$MAPPING_SEED" \
    --env "$PIXI_ENV"

# ---------------------------------------------------------------------------
# Step 2: batch_extract_cb_ub
# ---------------------------------------------------------------------------
step2_args=(
    --input-dir  "$R2_DIR"
    --output-dir "$R1_DIR"
    --qual-char  "$EXTRACT_QUAL_CHAR"
    --mapping-csv "$MAPPING_CSV"
    --env "$PIXI_ENV"
)
run_or_submit 2 "batch_extract_cb_ub" "$SCRIPT_DIR/02.batch_extract_cb_ub.sh" \
    "${step2_args[@]}"

# ---------------------------------------------------------------------------
# Step 3: run_splash
# ---------------------------------------------------------------------------
run_or_submit 3 "run_splash" "$SCRIPT_DIR/03.run_splash.sh" \
    --splash-dir  "$SPLASH_DIR" \
    --r1-dir      "$R1_DIR" \
    --r2-dir      "$R2_DIR" \
    --out-dir     "$SPLASH_OUT_DIR" \
    --threads     "$SPLASH_THREADS" \
    --n-bins      "$SPLASH_N_BINS" \
    --cb-len      "$SPLASH_CB_LEN" \
    --ub-len      "$SPLASH_UB_LEN" \
    --compactors-top-n "$SPLASH_COMPACTORS_TOP_N" \
    --compactors-max-length "$SPLASH_COMPACTORS_MAX_LENGTH" \
    --compactors-reads-buffer-gb "$SPLASH_COMPACTORS_READS_BUFFER_GB" \
    --env         "$PIXI_ENV"

# ---------------------------------------------------------------------------
# Step 4: matrix (build + filter)
# ---------------------------------------------------------------------------
step4_args=(
    --out-dir       "$SPLASH_OUT_DIR"
    --satc-dump-bin "$SATC_DUMP_BIN"
    --whitelist     "$WHITELIST"
    --cbc-mapping-csv "$MAPPING_CSV"
    --anchor-len    "$ANCHOR_LEN"
    --target-len    "$TARGET_LEN"
    --min-total-count "$MIN_TOTAL_COUNT"
    --min-features  "$MIN_FEATURES"
)
if [[ "$KEEP_DUMP_TEXT" == "true" ]]; then
    step4_args+=(--keep-dump-text)
fi
run_or_submit 4 "matrix" "$SCRIPT_DIR/04.build_filter_matrix.sh" "${step4_args[@]}" \
    --env "$PIXI_ENV"

# ---------------------------------------------------------------------------
# Step 5: anchor_align_genome
# ---------------------------------------------------------------------------
step5_args=(
    --features     "$MATRIX_FILTERED_COUNTS_FEATURES"
    --out-dir      "$ANCHOR_ALIGN_OUT_DIR"
    --bowtie-index "$BOWTIE_INDEX"
    --genome       "$GENOME_FA"
    --genome-gtf   "$GENOME_GTF"
    --threads      "$ALIGN_THREADS"
    --kmer         "$ALIGN_KMER"
    --window       "$ALIGN_WINDOW"
    --qual-char    "$ALIGN_QUAL_CHAR"
    --max-alignments-per-read "$MAX_ALIGNMENTS_PER_READ"
)
if [[ "$KEEP_DUPLICATE_TARGETS" == "true" ]]; then
    step5_args+=(--keep-duplicate-targets)
fi
run_or_submit 5 "anchor_align_genome" "$SCRIPT_DIR/05.anchor_align_genome.sh" \
    "${step5_args[@]}" \
    --env "$PIXI_ENV"

# ---------------------------------------------------------------------------
# Step 6: anchor (parse + HTML) (it is not suitable for step05)
# ---------------------------------------------------------------------------
# run_or_submit 6 "anchor" "$SCRIPT_DIR/06.anchor_parse.sh" \
#     --assigned-reads  "$ASSIGNED_READS" \
#     --features-file   "$MATRIX_FILTERED_FEATURES" \
#     --out-dir         "$PARSED_OUT_DIR" \
#     --merge-distance  "$MERGE_DISTANCE" \
#     --html-out-dir    "$HTML_OUT_DIR" \
#     --pattern         "$COMPARE_PATTERN" \
#     --min-common-substr "$COMPARE_MIN_COMMON_SUBSTR" \
#     --min-identity    "$COMPARE_MIN_IDENTITY" \
#     --max-gap-runs    "$COMPARE_MAX_GAP_RUNS" \
#     --max-fragments   "$COMPARE_MAX_FRAGMENTS" \
#     --max-errors      "$COMPARE_MAX_ERRORS" \
#     --workers         "$COMPARE_WORKERS" \
#     --env             "$PIXI_ENV"

# ---------------------------------------------------------------------------
# Step 7: anchor_entropy (remap + entropy + annotate)
# ---------------------------------------------------------------------------
# run_or_submit 7 "anchor_entropy" "$SCRIPT_DIR/07.anchor_entropy.sh" \
#     --matrix-dir       "$MATRIX_FILTERED_OUT_DIR" \
#     --mapping-dir      "$HTML_OUT_DIR" \
#     --gtf              "$GENOME_GTF" \
#     --remapped-out-dir "$REMAPPED_OUT_DIR" \
#     --anchor-align-dir "$ANCHOR_ALIGN_OUT_DIR" \
#     --min-total-count  "$ANCHOR_ENTROPY_MIN_TOTAL_COUNT" \
#     --env              "$PIXI_ENV"

# ---------------------------------------------------------------------------
if [[ "$MODE" == "local" ]]; then
    echo "[OK] Pipeline finished."
else
    echo "[OK] All steps submitted to SLURM."
    echo "[OK] Monitor with: squeue -u \$USER"
    echo "[OK] Logs: $LOG_DIR"
fi
