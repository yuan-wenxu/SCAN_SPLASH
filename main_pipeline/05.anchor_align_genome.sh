#!/usr/bin/env bash
#SBATCH -J anchor_align
#SBATCH -c 8
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/python/05.extract_anchor_fastq.py"

FEATURES=""
OUT_DIR=""
GENOME=""
GENOME_GTF=""
BOWTIE_INDEX=""
THREADS="${SLURM_CPUS_PER_TASK:-6}"
KMER="11"
WINDOW="5"
QUAL_CHAR="I"
KEEP_DUPLICATE_ANCHORS="false"
PIXI_ENV="main"
MAX_ALIGNMENTS_PER_READ="20"
JOBS="2"

usage() {
  cat <<'EOF'
Usage:
  anchor_align_genome.sh --features FILE --genome FILE --genome-gtf FILE [options]

Required arguments:
  --features FILE      features.tsv produced by filtered matrix
  --bowtie-index FILE  Bowtie index prefix for genome alignment
  --genome FILE        Reference genome FASTA
  --genome-gtf FILE    Reference annotation GTF

Optional arguments (with defaults):
  --out-dir DIR        Output directory (default: dirname of --features)
  --threads INT        Threads (default: SLURM_CPUS_PER_TASK or 8)
  --kmer INT           minimap2 -k value (default: 11)
  --window INT         minimap2 -w value (default: 5)
  --qual-char CHAR     FASTQ quality character (default: I)
  --keep-duplicate-targets
                       Keep duplicate targets in per-anchor target FASTQ
  --max-alignments-per-read INT
                       Keep up to N alignments per read in bowtie (default: 20)
  --jobs INT           Number of FASTQ files processed in parallel (default: 12)
  --env STR           Pixi environment name (default: main)
  -h, --help           Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --features)
      FEATURES="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --bowtie-index)
      BOWTIE_INDEX="$2"
      shift 2
      ;;
    --genome)
      GENOME="$2"
      shift 2
      ;;
    --genome-gtf)
      GENOME_GTF="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --kmer)
      KMER="$2"
      shift 2
      ;;
    --window)
      WINDOW="$2"
      shift 2
      ;;
    --qual-char)
      QUAL_CHAR="$2"
      shift 2
      ;;
    --keep-duplicate-targets)
      KEEP_DUPLICATE_ANCHORS="true"
      shift
      ;;
    --max-alignments-per-read)
      MAX_ALIGNMENTS_PER_READ="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
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
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$FEATURES" || -z "$BOWTIE_INDEX" || -z "$GENOME_GTF" ]]; then
  echo "[ERROR] --features, --bowtie-index, and --genome-gtf are required." >&2
  usage >&2
  exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$(dirname "$FEATURES")"
fi

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [[ "$JOBS" -lt 1 ]]; then
  echo "[ERROR] --jobs must be a positive integer: $JOBS" >&2
  exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "[ERROR] --threads must be a positive integer: $THREADS" >&2
  exit 1
fi

if [[ ! -f "$FEATURES" ]]; then
  echo "[ERROR] features file not found: $FEATURES" >&2
  exit 1
fi
if [[ ! -f "$BOWTIE_INDEX".1.ebwt ]]; then
  echo "[ERROR] Bowtie index not found: $BOWTIE_INDEX" >&2
  exit 1
fi
if [[ ! -f "$GENOME_GTF" ]]; then
  echo "[ERROR] genome GTF file not found: $GENOME_GTF" >&2
  exit 1
fi
if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "[ERROR] Python extractor not found: $PY_SCRIPT" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
FASTQ_GZ="$OUT_DIR/anchors.fastq.gz"
TARGET_FASTQ_DIR="$OUT_DIR/anchor_targets_fastq"

extract_cmd=(
  pixi run -e "$PIXI_ENV" python "$PY_SCRIPT"
  --features "$FEATURES"
  --out-dir "$OUT_DIR"
  --qual-char "$QUAL_CHAR"
)
if [[ "$KEEP_DUPLICATE_ANCHORS" == "true" ]]; then
  extract_cmd+=(--keep-duplicate-targets)
fi

"${extract_cmd[@]}"

if [[ ! -s "$FASTQ_GZ" ]]; then
  echo "[ERROR] extracted FASTQ not found or empty: $FASTQ_GZ" >&2
  exit 1
fi
if [[ ! -d "$TARGET_FASTQ_DIR" ]]; then
  echo "[ERROR] target FASTQ directory not found: $TARGET_FASTQ_DIR" >&2
  exit 1
fi

FASTQ_FILES=("$FASTQ_GZ")
while IFS= read -r -d '' fq; do
  FASTQ_FILES+=("$fq")
done < <(find "$TARGET_FASTQ_DIR" -maxdepth 1 -type f -name '*.fastq.gz' -print0 | sort -z)

if [[ "${#FASTQ_FILES[@]}" -eq 0 ]]; then
  echo "[ERROR] no FASTQ files found for alignment" >&2
  exit 1
fi

MERGED_FASTQ_GZ="$OUT_DIR/merged_anchor_reads.fastq.gz"
MERGED_BAM="$OUT_DIR/merged_anchor_reads.bam"
MERGED_FEATURECOUNTS="$OUT_DIR/merged_anchor_reads.featurecounts.txt"

echo "[INFO] Merging ${#FASTQ_FILES[@]} FASTQ files -> $MERGED_FASTQ_GZ"
cat "${FASTQ_FILES[@]}" > "$MERGED_FASTQ_GZ"

if [[ ! -s "$MERGED_FASTQ_GZ" ]]; then
  echo "[ERROR] merged FASTQ not found or empty: $MERGED_FASTQ_GZ" >&2
  exit 1
fi

echo "[INFO] Aligning merged FASTQ with bowtie (k=$MAX_ALIGNMENTS_PER_READ)"
{
pixi run -e "$PIXI_ENV" bowtie \
  -q \
  -v 1 \
  -p "$THREADS" \
  -k "$MAX_ALIGNMENTS_PER_READ" \
  --best \
  --strata \
  -S \
  "$BOWTIE_INDEX" \
  "$MERGED_FASTQ_GZ" |
pixi run -e "$PIXI_ENV" samtools view -bS -F 4 - | \
pixi run -e "$PIXI_ENV" samtools sort -o "$MERGED_BAM" -
} >> "$OUT_DIR/bowtie.log" 2>&1

pixi run -e "$PIXI_ENV" samtools index "$MERGED_BAM" >> "$OUT_DIR/bowtie.log" 2>&1

pixi run -e "$PIXI_ENV" featureCounts -a "$GENOME_GTF" -o "$MERGED_FEATURECOUNTS" -T "$THREADS" -R BAM "$MERGED_BAM" >> "$OUT_DIR/featureCounts.log" 2>&1

echo "[OK] fastq=$FASTQ_GZ"
echo "[OK] merged_fastq=$MERGED_FASTQ_GZ"
echo "[OK] merged_bam=$MERGED_BAM"
echo "[OK] featureCounts=$MERGED_FEATURECOUNTS"
