#!/usr/bin/env bash
#SBATCH -J splash_pseudo_r1
#SBATCH -c 8
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

# Build SPLASH input from pseudo-R1 files plus original barcode-tagged R2 files,
# then run SPLASH in 10x mode if CBC/UMI lengths are supported.

SPLASH_DIR=""
R1_DIR=""
R2_DIR=""
OUT_DIR=""

THREADS="${SLURM_CPUS_PER_TASK:-8}"
N_BINS="64"
CB_LEN="16"
UB_LEN="8"
PIXI_ENV="main"
COMPACTORS_TOP_N="1000"
COMPACTORS_NUM_KMERS="2"
COMPACTORS_KMER_LEN="27"
COMPACTORS_MAX_LENGTH="186"
COMPACTORS_READS_BUFFER_GB="8"
COMPACTORS_MIN_EXTENDER_SPECIFICITY="0.9"

usage() {
  cat <<'EOF'
Usage:
  run_splash.sh --splash-dir DIR --r1-dir DIR --r2-dir DIR --out-dir DIR

Required path arguments:
  --splash-dir DIR   Directory containing SPLASH binaries and splash script
  --r1-dir DIR       Directory containing pseudo-R1 FASTQ files
  --r2-dir DIR       Directory containing original barcode-tagged R2 FASTQ files
  --out-dir DIR      Output directory

Optional arguments (with defaults):
  --threads INT      Number of threads (default: SLURM_CPUS_PER_TASK or 8)
  --n-bins INT       Number of bins (default: 64)
  --cb-len INT       CBC length (default: 16)
  --ub-len INT       UMI length (default: 8)
  --compactors-top-n INT
                     Number of anchors used for each compactors run
                     (top target entropy and top effect size; default: 1000)
  --compactors-max-length INT
                     Maximum compactor length in bases (default: 186)
  --compactors-reads-buffer-gb INT
                     Read buffer size for compactors (default: 8)
  -h, --help         Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --splash-dir)
      SPLASH_DIR="$2"
      shift 2
      ;;
    --r1-dir)
      R1_DIR="$2"
      shift 2
      ;;
    --r2-dir)
      R2_DIR="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --n-bins)
      N_BINS="$2"
      shift 2
      ;;
    --cb-len)
      CB_LEN="$2"
      shift 2
      ;;
    --ub-len)
      UB_LEN="$2"
      shift 2
      ;;
    --env)
      PIXI_ENV="$2"
      shift 2
      ;;
    --compactors-top-n)
      COMPACTORS_TOP_N="$2"
      shift 2
      ;;
    --compactors-max-length)
      COMPACTORS_MAX_LENGTH="$2"
      shift 2
      ;;
    --compactors-reads-buffer-gb)
      COMPACTORS_READS_BUFFER_GB="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$SPLASH_DIR" || -z "$R1_DIR" || -z "$R2_DIR" || -z "$OUT_DIR" ]]; then
  echo "[ERROR] --splash-dir, --r1-dir, --r2-dir, and --out-dir are required." >&2
  usage >&2
  exit 1
fi

if [[ ! -d "$R1_DIR" ]]; then
  echo "[ERROR] R1 directory not found: $R1_DIR" >&2
  exit 1
fi

if [[ ! -d "$R2_DIR" ]]; then
  echo "[ERROR] R2 directory not found: $R2_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/logs"

COMPACTORS_CONFIG="$OUT_DIR/compactors_config.json"
cat > "$COMPACTORS_CONFIG" <<EOF
{
  "num_threads": $THREADS,
  "epsilon": 0.001,
  "beta": 0.5,
  "num_kmers": $COMPACTORS_NUM_KMERS,
  "kmer_len": $COMPACTORS_KMER_LEN,
  "max_length": $COMPACTORS_MAX_LENGTH,
  "min_extender_specificity": $COMPACTORS_MIN_EXTENDER_SPECIFICITY,
  "reads_buffer_gb": $COMPACTORS_READS_BUFFER_GB
}
EOF

mapfile -t R2_FILES < <(find "$R2_DIR" -maxdepth 1 -type f -name '*.fastq.gz' | LC_ALL=C sort)
if [[ ${#R2_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No R2 FASTQ files found in $R2_DIR" >&2
  exit 1
fi

PAIR_LIST="$OUT_DIR/pairs.txt"
INPUT_LIST="$OUT_DIR/input.txt"

: > "$PAIR_LIST"
for r2 in "${R2_FILES[@]}"; do
  bn="$(basename "$r2")"
  r1="$R1_DIR/${bn%.fastq.gz}.R1.fastq.gz"
  if [[ ! -f "$r1" ]]; then
    echo "[ERROR] Matching R1 file not found for $r2" >&2
    echo "       Expected: $r1" >&2
    echo "[TIP] First generate pseudo-R1 files." >&2
    exit 1
  fi
  printf '%s,%s\n' "$r1" "$r2" >> "$PAIR_LIST"
done

printf '%s %s\n' sample1 "$PAIR_LIST" > "$INPUT_LIST"

FIRST_R2="$(sed -n '1s/^[^,]*,//p' "$PAIR_LIST")"
if [[ -z "$FIRST_R2" ]] || [[ ! -f "$FIRST_R2" ]]; then
  echo "[ERROR] Cannot determine first R2 from $PAIR_LIST" >&2
  exit 1
fi

pixi run -e "$PIXI_ENV" python "$SPLASH_DIR/splash" \
  "$INPUT_LIST" \
  --outname_prefix "$OUT_DIR/result" \
  --gap_len 0 \
  --technology 10x \
  --poly_ACGT_len 6 \
  --anchor_count_threshold 50 \
  --anchor_samples_threshold 1 \
  --anchor_sample_counts_threshold 5 \
  --fdr_threshold 0.05 \
  --keep_significant_anchors_satc \
  --export_filtered_input \
  --keep_top_n_target_entropy "$COMPACTORS_TOP_N" \
  --keep_top_n_effect_size_bin "$COMPACTORS_TOP_N" \
  --compactors_config "$COMPACTORS_CONFIG" \
  --dump_Cjs \
  --bin_path "$SPLASH_DIR" \
  --n_threads_stage_1 1 \
  --n_threads_stage_1_internal "$THREADS" \
  --n_threads_stage_2 "$THREADS" \
  --n_bins "$N_BINS" \
  --logs_dir "$OUT_DIR/logs" \
  --cbc_len "$CB_LEN" \
  --umi_len "$UB_LEN" \
  --cbc_filtering_thr 100 \
  --export_cbc_logs > "$OUT_DIR/splash.log" 2>&1
