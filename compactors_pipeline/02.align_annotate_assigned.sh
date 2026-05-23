#!/usr/bin/env bash
#SBATCH -J compactors_align_annotate
#SBATCH -c 8
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAPQ_MIN="30"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_DIR="${SCRIPT_DIR}/python"
TO_FASTQ_SCRIPT="${PY_DIR}/02.compactors_tsv_to_fastq.py"
ASSIGNED_SCRIPT="${PY_DIR}/02.assigned_reads_to_anchor_gene_csv.py"
PIXI_ENV="main"
MERGE_DISTANCE="15"

INPUT_TSV=""
SCORES_TSV=""
STAR_GENOME_DIR=""
OUT_DIR=""

usage() {
    cat >&2 <<EOF
Usage: $0 --input-tsv FILE --star-genome-dir DIR --gtf FILE --out-dir DIR [options]

Pipeline:
  1) Run 02.compactors_tsv_to_fastq.py (TSV -> FASTQ.GZ)
  2) STAR alignment
  3) featureCounts gene annotation
  4) Extract MAPQ > threshold and Assigned records to assigned_reads_<MAPQ>.txt
  5) Run 02.assigned_reads_to_anchor_gene_csv.py

Required:
  --input-tsv FILE         Compactors TSV from 01.run_compactors.sh
  --scores-tsv FILE        Scores TSV for input tsv
  --star-genome-dir DIR    STAR genome directory for alignment
  --out-dir DIR            Output directory

Optional:
  --threads INT            Threads (default: SLURM_CPUS_PER_TASK or 8)
  --mapq-min INT           Minimum MAPQ for assigned read extraction (default: 30)
  --to-fastq-script FILE   Path to 02.compactors_tsv_to_fastq.py
  --assigned-script FILE   Path to 02.assigned_reads_to_anchor_gene_csv.py
  -h, --help               Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-tsv)
            INPUT_TSV="$2"
            shift 2
            ;;
        --scores-tsv)
            SCORES_TSV="$2"
            shift 2
            ;;
        --star-genome-dir)
            STAR_GENOME_DIR="$2"
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
        --mapq-min)
            MAPQ_MIN="$2"
            shift 2
            ;;
        --to-fastq-script)
            TO_FASTQ_SCRIPT="$2"
            shift 2
            ;;
        --assigned-script)
            ASSIGNED_SCRIPT="$2"
            shift 2
            ;;
        --merge-distance)
            MERGE_DISTANCE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_TSV" || -z "$SCORES_TSV" || -z "$STAR_GENOME_DIR" || -z "$OUT_DIR" ]]; then
    echo "[ERROR] Missing required arguments" >&2
    usage
    exit 1
fi

if [[ ! -f "$INPUT_TSV" ]]; then
    echo "[ERROR] input TSV not found: $INPUT_TSV" >&2
    exit 1
fi
if [[ ! -f "$SCORES_TSV" ]]; then
    echo "[ERROR] scores TSV not found: $SCORES_TSV" >&2
    exit 1
fi
if [[ ! -d "$STAR_GENOME_DIR" ]]; then
    echo "[ERROR] STAR genome directory not found: $STAR_GENOME_DIR" >&2
    exit 1
fi
if [[ ! -f "$TO_FASTQ_SCRIPT" ]]; then
    echo "[ERROR] 02 script not found: $TO_FASTQ_SCRIPT" >&2
    exit 1
fi
if [[ ! -f "$ASSIGNED_SCRIPT" ]]; then
    echo "[ERROR] assigned script not found: $ASSIGNED_SCRIPT" >&2
    exit 1
fi

mkdir -p "$OUT_DIR" "$OUT_DIR/logs"

BASE_NAME="$(basename "$INPUT_TSV" .tsv)"
FASTQ_GZ="$OUT_DIR/${BASE_NAME}.fastq.gz"
BAM_PATH="$OUT_DIR/star/Aligned.sortedByCoord.out.bam"
ASSIGNED_TXT="$OUT_DIR/assigned_reads.txt"

echo "[INFO] Step 1/3: compactors TSV -> FASTQ.GZ" >&2
pixi run -e "$PIXI_ENV" "$TO_FASTQ_SCRIPT" \
    --input-tsv "$INPUT_TSV" \
    --scores-tsv "$SCORES_TSV" \
    --output-fastq-gz "$FASTQ_GZ" \

echo "[INFO] Step 2/3: star alignment" >&2

pixi run -e "$PIXI_ENV" star \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_GENOME_DIR" \
    --readFilesCommand zcat \
    --readFilesIn "$FASTQ_GZ" \
    --outFileNamePrefix "$OUT_DIR/star/" \
    --outTmpDir "$OUT_DIR/solotmp" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM GX GN > "$OUT_DIR/logs/02.star.log" 2>&1

pixi run -e "$PIXI_ENV" samtools view -@ "$THREADS" -q "$MAPQ_MIN" "$BAM_PATH" | \
awk 'BEGIN { FS="\t" }
{
    if ($3 == "chrM") {
        next
    }

    gx = ""
    gn = ""
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^GX:Z:/) {
            gx = substr($i, 6)
        } else if ($i ~ /^GN:Z:/) {
            gn = substr($i, 6)
        }
    }

    if (gn ~ /^DARLIN-(CA|RA|TA)$/) {
        next
    }

    if (gx != "" && gx != "-" && gn != "" && gn != "-") {
        print
    }
}' > "$ASSIGNED_TXT"

echo "[INFO] assigned records: $(wc -l < "$ASSIGNED_TXT")" >&2

echo "[INFO] Step 3/3: assigned CSV" >&2
pixi run -e "$PIXI_ENV" "$ASSIGNED_SCRIPT" \
    --input "$ASSIGNED_TXT" \
    --output "$ASSIGNED_TXT.csv" \
    --merge-distance "$MERGE_DISTANCE" \
