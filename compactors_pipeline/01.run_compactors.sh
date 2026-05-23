#!/usr/bin/env bash
#SBATCH -J compactors
#SBATCH -c 8
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

# Run SPLASH compactors on the current filtered/deduplicated input reads.

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-8}"
NUM_KMERS="2"
KMER_LEN="27"
MAX_LENGTH="186"
EPSILON="0.001"
BETA="0.5"
LOWER_BOUND="10"
MAX_MISMATCH="4"
MIN_EXTENDER_SPECIFICITY="0.9"
NUM_EXTENDERS="1"
EXTENDERS_SHIFT="1"
MAX_ANCHOR_COMPACTORS="1000"
MAX_CHILD_COMPACTORS="20"
READS_BUFFER_GB="64"
RUN_ALL_SCORES=0

usage() {
  cat >&2 <<EOF
Usage: $0 [options]

Options:
  --splash-dir DIR                 SPLASH binary directory.
  --result-dir DIR                 SPLASH result directory containing result.after_correction.*.tsv.
  --filtered-input-dir DIR         Directory containing result_filtered_input FASTQ files.
  --out-dir DIR                    Output directory for manual compactors results.
  --threads INT                    Threads for compactors. Default: SLURM_CPUS_PER_TASK or 8.
  --num-kmers INT                  Number of downstream k-mers per extension segment. Default: 2.
  --kmer-len INT                   Length of each downstream k-mer. Default: 27.
  --max-length INT                 Maximum compactor length in bases. Default: 186.
  --epsilon FLOAT                  Sequencing error rate used by compactors. Default: 0.001.
  --beta FLOAT                     Error-model/stringency multiplier. Default: 0.5.
  --lower-bound INT                Minimum k-mer abundance for active-set candidates. Default: 10.
  --max-mismatch INT               Maximum accumulated mismatches for candidate support. Default: 4.
  --min-extender-specificity FLOAT Minimum specificity required for recursive extension. Default: 0.9.
  --reads-buffer-gb INT            Read buffer size in GB. Default: 8.
  -h, --help                       Show this help message.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --splash-dir)
      SPLASH_DIR="$2"
      shift 2
      ;;
    --result-dir)
      RESULT_DIR="$2"
      FILTERED_INPUT_DIR="$RESULT_DIR/result_filtered_input"
      OUT_DIR="$RESULT_DIR/result_compactors_manual"
      shift 2
      ;;
    --filtered-input-dir)
      FILTERED_INPUT_DIR="$2"
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
    --num-kmers)
      NUM_KMERS="$2"
      shift 2
      ;;
    --kmer-len)
      KMER_LEN="$2"
      shift 2
      ;;
    --max-length)
      MAX_LENGTH="$2"
      shift 2
      ;;
    --epsilon)
      EPSILON="$2"
      shift 2
      ;;
    --beta)
      BETA="$2"
      shift 2
      ;;
    --lower-bound)
      LOWER_BOUND="$2"
      shift 2
      ;;
    --max-mismatch)
      MAX_MISMATCH="$2"
      shift 2
      ;;
    --min-extender-specificity)
      MIN_EXTENDER_SPECIFICITY="$2"
      shift 2
      ;;
    --reads-buffer-gb)
      READS_BUFFER_GB="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

COMPACTORS="$SPLASH_DIR/compactors"
if [[ ! -x "$COMPACTORS" ]]; then
  echo "[ERROR] compactors binary not found or not executable: $COMPACTORS" >&2
  exit 1
fi

if [[ ! -d "$FILTERED_INPUT_DIR" ]]; then
  echo "[ERROR] filtered input directory not found: $FILTERED_INPUT_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR/logs"

FASTQ_LIST="$OUT_DIR/filtered_fastq.list"
find "$FILTERED_INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) | LC_ALL=C sort > "$FASTQ_LIST"

if [[ ! -s "$FASTQ_LIST" ]]; then
  echo "[ERROR] no FASTQ files found in: $FILTERED_INPUT_DIR" >&2
  exit 1
fi

run_compactors_one() {
  local label="$1"
  local anchors_tsv="$2"
  local out_tsv="$OUT_DIR/${label}.compactors.tsv"
  local out_fasta="$OUT_DIR/${label}.compactors.fasta"
  local log="$OUT_DIR/logs/${label}.compactors.log"

  if [[ ! -f "$anchors_tsv" ]]; then
    echo "[WARN] skip $label, anchors TSV not found: $anchors_tsv" >&2
    return
  fi

  echo "[INFO] Running compactors: $label" >&2
  echo "[INFO] anchors: $anchors_tsv" >&2
  echo "[INFO] output:  $out_tsv" >&2

  "$COMPACTORS" \
    --input_format fastq \
    --num_threads "$THREADS" \
    --num_kmers "$NUM_KMERS" \
    --kmer_len "$KMER_LEN" \
    --epsilon "$EPSILON" \
    --beta "$BETA" \
    --lower_bound "$LOWER_BOUND" \
    --max_mismatch "$MAX_MISMATCH" \
    --max_length "$MAX_LENGTH" \
    --min_extender_specificity "$MIN_EXTENDER_SPECIFICITY" \
    --num_extenders "$NUM_EXTENDERS" \
    --extenders_shift "$EXTENDERS_SHIFT" \
    --max_anchor_compactors "$MAX_ANCHOR_COMPACTORS" \
    --max_child_compactors "$MAX_CHILD_COMPACTORS" \
    --reads_buffer_gb "$READS_BUFFER_GB" \
    --out_fasta "$out_fasta" \
    --log "$log" \
    "$FASTQ_LIST" \
    "$anchors_tsv" \
    "$out_tsv" \
    > "$OUT_DIR/logs/${label}.stdout.log" 2>&1
}

run_compactors_one \
  "all_scores" \
  "$RESULT_DIR/result.after_correction.scores.tsv"

echo "[DONE] compactors outputs written to: $OUT_DIR" >&2
