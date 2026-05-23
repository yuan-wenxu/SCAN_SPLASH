#!/usr/bin/env bash
# config.sh — SCAN-SPLASH pipeline configuration
# Edit this file and pass it with: bash main.sh --config config.sh
#
# Pipeline steps (main.sh --steps):
#   1) barcode_mapping      barcode_mapping.sh
#   2) batch_extract_cb_ub  batch_extract_cb_ub.sh
#   3) run_splash           run_splash.sh
#   4) matrix               matrix.sh  (build + filter)
#   5) anchor_align_genome  anchor_align_genome.sh
#   6) anchor               anchor.sh  (parse + HTML)
#   7) anchor_entropy       anchor_entropy.sh (remap + entropy + annotate)

# ---------------------------------------------------------------------------
# Required — must be set before running
# ---------------------------------------------------------------------------
BARCODE_INPUT=""                   # Input barcode combo TSV (e.g. 4CL.tsv)
WHITELIST=""                       # Whitelist TSV for matrix filtering
R2_DIR=""                          # Directory containing R2 FASTQ files
SPLASH_DIR=""                      # Directory containing the SPLASH binaries (satc_dump, splash, …)
GENOME_FA=""                       # Reference genome FASTA (not used)
BOWTIE_INDEX=""                    # Bowtie index prefix
BOWTIE_TRANSCRIPTOME_INDEX=""      # Bowtie index prefix for transcriptome
GENOME_GTF=""                      # Reference genome GTF / DARLIN GTF

# ---------------------------------------------------------------------------
# Output and naming
# ---------------------------------------------------------------------------
PIPELINE_OUT_DIR=""
PIXI_ENV="main"

# ---------------------------------------------------------------------------
# Step 1: barcode_mapping
# ---------------------------------------------------------------------------
MAPPING_CODE_LEN="16"
MAPPING_POOL_SIZE="800"
MAPPING_SEED="42"

# ---------------------------------------------------------------------------
# Step 2: batch_extract_cb_ub
# ---------------------------------------------------------------------------
EXTRACT_QUAL_CHAR="I"

# ---------------------------------------------------------------------------
# Step 3: run_splash
# ---------------------------------------------------------------------------
SPLASH_THREADS="${SLURM_CPUS_PER_TASK:-8}"
SPLASH_N_BINS="64"
SPLASH_CB_LEN="16"
SPLASH_UB_LEN="8"
SPLASH_COMPACTORS_TOP_N="1000"
SPLASH_COMPACTORS_MAX_LENGTH="186"
SPLASH_COMPACTORS_READS_BUFFER_GB="8"

# ---------------------------------------------------------------------------
# Step 4: matrix (build_cell_anchor_target_matrix + filter)
# ---------------------------------------------------------------------------
ANCHOR_LEN="31"
TARGET_LEN="31"
KEEP_DUMP_TEXT="false"
MIN_TOTAL_COUNT="20000"
MIN_FEATURES="10000"

# ---------------------------------------------------------------------------
# Step 5: anchor_align_genome
# ---------------------------------------------------------------------------
ALIGN_THREADS="${SLURM_CPUS_PER_TASK:-8}"
ALIGN_KMER="11"
ALIGN_WINDOW="5"
ALIGN_QUAL_CHAR="I"
MAX_ALIGNMENTS_PER_READ="20"
KEEP_DUPLICATE_TARGETS="false"

# ---------------------------------------------------------------------------
# Step 6: anchor (parse_assigned_reads + compare_targets_indel)
# ---------------------------------------------------------------------------
MERGE_DISTANCE="15"
COMPARE_PATTERN="*.txt"
COMPARE_MIN_COMMON_SUBSTR="5"
COMPARE_MIN_IDENTITY="0.5"
COMPARE_MAX_GAP_RUNS="2"
COMPARE_MAX_FRAGMENTS="3"
COMPARE_MAX_ERRORS="6"
COMPARE_WORKERS="8"

# ---------------------------------------------------------------------------
# Step 7: anchor_entropy (remap_filtered_matrix_targets + entropy + annotate)
# ---------------------------------------------------------------------------
ANCHOR_ENTROPY_MIN_TOTAL_COUNT="1"

export BARCODE_INPUT WHITELIST R2_DIR SPLASH_DIR GENOME_FA GENOME_GTF
export BOWTIE_INDEX BOWTIE_TRANSCRIPTOME_INDEX
export PIPELINE_OUT_DIR PIXI_ENV
export MAPPING_CODE_LEN MAPPING_POOL_SIZE MAPPING_SEED
export EXTRACT_QUAL_CHAR
export SPLASH_THREADS SPLASH_N_BINS SPLASH_CB_LEN SPLASH_UB_LEN
export SPLASH_COMPACTORS_TOP_N SPLASH_COMPACTORS_MAX_LENGTH SPLASH_COMPACTORS_READS_BUFFER_GB
export ANCHOR_LEN TARGET_LEN KEEP_DUMP_TEXT
export MIN_TOTAL_COUNT MIN_FEATURES
export ALIGN_THREADS ALIGN_KMER ALIGN_WINDOW ALIGN_QUAL_CHAR MAX_ALIGNMENTS_PER_READ KEEP_DUPLICATE_TARGETS
export MERGE_DISTANCE
export COMPARE_PATTERN COMPARE_MIN_COMMON_SUBSTR COMPARE_MIN_IDENTITY COMPARE_MAX_GAP_RUNS
export COMPARE_MAX_FRAGMENTS COMPARE_MAX_ERRORS COMPARE_WORKERS
export ANCHOR_ENTROPY_MIN_TOTAL_COUNT