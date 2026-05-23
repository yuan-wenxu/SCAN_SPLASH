#!/usr/bin/env bash
#SBATCH -J build_pseudo_r1
#SBATCH -c 1
#SBATCH -p amd-ep2,intel-sc3
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --requeue

# writes pseudo-R1 FASTQ outputs. The workload is I/O-bound and single-threaded,
# so 1 CPU and 8G RAM are sufficient; walltime is padded to 24h for cluster I/O.

set -euo pipefail

# Batch wrapper for extract_cb_ub_fastq.py
#
# Usage:
#   bash batch_extract_cb_ub.sh -i <input_dir> -o <output_dir> [-m <mapping_csv>] [-q <qual_char>]
#
# It scans only top-level files in <input_dir> with extensions:
#   *.fastq, *.fq, *.fastq.gz, *.fq.gz
#
# For each R2 file, outputs:
#   <output_dir>/<basename>.R1.fastq.gz

usage() {
  cat >&2 <<EOF
Usage: $0 -i <input_dir> -o <output_dir> [-m <mapping_csv>] [-q <qual_char>]

Options:
  -i, --input-dir     Directory containing input FASTQ files.
  -m, --mapping-csv   Optional CSV file to remap CB combinations.
  -o, --output-dir    Directory for generated pseudo-R1 FASTQ files.
  -q, --qual-char     Quality character used for synthetic CB/UB bases. Default: I
  -e, --pixi-env      Pixi environment name. Default: main
      --env           Alias of --pixi-env
  -h, --help          Show this help message.
EOF
}

mapping_csv=""
qual_char="I"
pixi_env="main"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input-dir)
      input_dir="${2:-}"
      shift 2
      ;;
    -m|--mapping-csv)
      mapping_csv="${2:-}"
      shift 2
      ;;
    -o|--output-dir)
      output_dir="${2:-}"
      shift 2
      ;;
    -q|--qual-char)
      qual_char="${2:-}"
      shift 2
      ;;
    -e|--pixi-env|--env)
      pixi_env="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$input_dir" || -z "$output_dir" ]]; then
  echo "[ERROR] -i and -o are required." >&2
  usage
  exit 1
fi

if [[ ! -d "$input_dir" ]]; then
  echo "[ERROR] input_dir not found: $input_dir" >&2
  exit 1
fi

if [[ -n "$mapping_csv" && ! -f "$mapping_csv" ]]; then
  echo "[ERROR] mapping_csv not found: $mapping_csv" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "$0")" && pwd)"
py_script="$script_dir/python/02.extract_cb_ub_fastq.py"

if [[ ! -f "$py_script" ]]; then
  echo "[ERROR] Python script not found: $py_script" >&2
  exit 1
fi

mkdir -p "$output_dir"

mapfile -t fastq_files < <(
  find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) | LC_ALL=C sort
)

if [[ ${#fastq_files[@]} -eq 0 ]]; then
  echo "[ERROR] No FASTQ files found in: $input_dir" >&2
  exit 1
fi

echo "[INFO] Found ${#fastq_files[@]} FASTQ files in $input_dir" >&2

for f in "${fastq_files[@]}"; do
  bn="$(basename "$f")"
  stem="$bn"
  stem="${stem%.fastq.gz}"
  stem="${stem%.fq.gz}"
  stem="${stem%.fastq}"
  stem="${stem%.fq}"

  out_fastq="$output_dir/${stem}.R1.fastq.gz"

  echo "[RUN] $bn" >&2
  cmd=(
    pixi run -e "$pixi_env" python "$py_script"
    --input "$f"
    --output-fastq "$out_fastq"
    --qual-char "$qual_char"
  )
  if [[ -n "$mapping_csv" ]]; then
    cmd+=(--mapping-csv "$mapping_csv")
  fi

  "${cmd[@]}"
done

echo "[INFO] Outputs written to: $output_dir" >&2
