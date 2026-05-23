# SCAN_SPLASH

> Status: this is an incomplete pipeline. The main workflow has been updated
> through Step 5, but Steps 6 and 7 are from an earlier pipeline version and
> have not been updated to match the current Step 5 outputs. As a result, Steps
> 6 and 7 cannot be connected directly after Step 5. The older version of the 
> pipeline can be found in `Alternative_Splicing` ​​repository.Due to project reasons,
> there are no current plans for further updates.

`SCAN_SPLASH` is a SPLASH-based downstream workflow for barcode-aware long-read
data in the `Alternative_Splicing` project. `Alternative_Splicing` is currently
a private, unpublished repository. This workflow starts from barcode-tagged R2
FASTQ files produced by the upstream pipeline and builds cell-level
anchor/target matrices and anchor alignments.

Alternatively, you can use a custom-constructed fastq file with readID: `readID+CB:ACGT+UB:ACGT`. CB means cell barcode and UB means UMI.

The main workflow lives in `main_pipeline/`. Two auxiliary workflows are also
provided:

- `Cjs_pipeline/`: build and filter a cell x anchor matrix from SPLASH
  `result_Cjs/cjs.tsv.gz`.
- `compactors_pipeline/`: rerun SPLASH compactors and align/annotate compactor
  sequences.

## When to Run

Run this folder after the upstream `Alternative_Splicing` barcode extraction
step has produced the Step 3 barcode-tagged R2 FASTQ directory, for example:

```text
.../pipeline/step03_barcode/
or custom-constructed fastq file 
```

Minimum inputs for the main workflow:

- Full barcode-combo whitelist table, used to generate `barcode_mapping.csv` (All possible barcode combinations).
- Selected/experiment whitelist table, used for matrix filtering (Barcodes in this experiment).
- Barcode-tagged R2 FASTQ directory.
- SPLASH installation directory containing `splash`, `satc_dump`, and
  `compactors`.
- Bowtie genome index and GTF for anchor/target annotation.

## Directory Layout

```text
SCAN_SPLASH/
|-- main_pipeline/
|   |-- 00.config.sh
|   |-- 00.main.sh
|   |-- 01.barcode_mapping.sh
|   |-- 02.batch_extract_cb_ub.sh
|   |-- 03.run_splash.sh
|   |-- 04.build_filter_matrix.sh
|   |-- 05.anchor_align_genome.sh
|   |-- 06.anchor_parse.sh
|   |-- 07.anchor_entropy.sh
|   `-- python/
|-- Cjs_pipeline/
|   |-- cj_downstream.sh
|   `-- python/
|-- compactors_pipeline/
|   |-- 01.run_compactors.sh
|   |-- 02.align_annotate_assigned.sh
|   `-- python/
`-- fast_align_cpp/
```

## Main Pipeline

Edit a config file first. The template is:

```bash
main_pipeline/00.config.sh
```

For a real run, copy or edit a config such as `main_pipeline/my_config.sh`, then
run from the `main_pipeline/` directory:

```bash
cd /mnt/work/project/SCAN_SPLASH/main_pipeline
bash 00.main.sh --config my_config.sh --mode local
```

Submit the same workflow to SLURM with dependency chaining:

```bash
bash 00.main.sh --config my_config.sh --mode hpc
```

Run selected steps only:

```bash
bash 00.main.sh --config my_config.sh --mode local --steps 1,2,3,4
```

### Required Config Variables

Set these variables before running:

```bash
BARCODE_INPUT=""      # Full barcode combo TSV/CSV
WHITELIST=""          # Selected whitelist TSV for filtering
R2_DIR=""             # Barcode-tagged R2 FASTQ directory
SPLASH_DIR=""         # SPLASH binary/script directory
GENOME_FA=""          # Reference FASTA; kept for compatibility
BOWTIE_INDEX=""       # Bowtie genome index prefix
GENOME_GTF=""         # Reference/DARLIN GTF
PIPELINE_OUT_DIR=""   # Output directory
PIXI_ENV="main"       # Pixi environment
```

Common tunables are defined in the same config: SPLASH thread/bin settings,
barcode length, UMI length, matrix QC thresholds, anchor length, target length,
and anchor annotation thresholds.

### Step Summary

`00.main.sh` defines these seven logical steps:

1. `barcode_mapping`: generate `barcode_mapping.csv` from the barcode-combo
   input table.
2. `batch_extract_cb_ub`: build pseudo-R1 FASTQ files from barcode-tagged R2
   reads.
3. `run_splash`: run SPLASH in 10x mode on pseudo-R1/R2 pairs.
4. `matrix`: dump SATC, build a cell x anchor-target matrix, filter cells by
   whitelist/QC, and filter anchors by detected-cell percentage.
5. `anchor_align_genome`: align anchor and target FASTQ records to the genome
   with Bowtie and annotate with featureCounts.
6. `anchor`: legacy step from an earlier pipeline version; parses assigned
   reads and generates per-anchor target comparison HTML.
7. `anchor_entropy`: legacy step from an earlier pipeline version; remaps
   clustered targets, computes entropy, and annotates entropy rows with gene
   names.

Current script status:

- `00.main.sh` currently executes Steps 1-5. The calls for Steps 6 and 7 are
  present but commented out.
- Steps 6 and 7 are retained as historical scripts from the previous pipeline
  design. They have not been updated after the Step 5 redesign, so their
  expected inputs do not line up with the files currently produced by Step 5.
  They should not be treated as a working continuation of the main pipeline.
- In `05.anchor_align_genome.sh`, the FASTQ extraction command is currently
  commented out. The script therefore expects these files/directories to already
  exist under the Step 5 output directory:
  - `anchors.fastq.gz`
  - `anchor_targets_fastq/`
- `05.anchor_align_genome.sh` currently writes `merged_anchor_reads.bam` and
  `merged_anchor_reads.featurecounts.txt`, but it does not create
  `assigned_reads.tsv`. Steps 6 and 7 require a SAM-format assigned reads file,
  so that file must be generated separately before those steps are used.


### Main Outputs

Given `PIPELINE_OUT_DIR=/path/to/splash`, the main workflow writes:

```text
/path/to/splash/
|-- barcode_mapping.csv
|-- r1/
|   `-- *.R1.fastq.gz
|-- logs/
`-- splash_results/
    |-- input.txt
    |-- pairs.txt
    |-- splash.log
    |-- result*.tsv / result*.satc / result_* outputs from SPLASH
    `-- matrix/
        |-- filtered/
        |   |-- whitelist_filtered.features.tsv
        |   |-- whitelist_filtered.cells.tsv
        |   |-- whitelist_filtered.matrix.mtx.gz
        |   |-- anchor_pct_filtered.features.tsv
        |   |-- anchor_pct_filtered.cells.tsv
        |   |-- anchor_pct_filtered.matrix.mtx.gz
        |   |-- anchor_align/
        |   `-- remapped/
        `-- ...
```

## Cj Downstream Pipeline

Use `Cjs_pipeline/cj_downstream.sh` to build a cell x anchor matrix directly
from SPLASH Cj output:

```bash
cd ./Cjs_pipeline

bash cj_downstream.sh \
  --cjs-tsv-gz /path/to/splash_results/result_Cjs/cjs.tsv.gz \
  --barcode-mapping-csv /path/to/splash/barcode_mapping.csv \
  --matrix-dir /path/to/splash/result_Cjs/matrix_from_cj \
  --whitelist /path/to/whitelist.tsv \
  --env main
```

Outputs:

```text
matrix_from_cj/
`-- filtered/
```

Useful options:

- `--min-anchors-per-cell INT` default `2000`
- `--min-cells-per-anchor INT` default `3`
- `--drop-unmapped`
- `--scatter-png FILE`

## Compactors Pipeline

Use `compactors_pipeline/01.run_compactors.sh` to run SPLASH compactors on the
filtered input reads exported by the main SPLASH run:

```bash
cd ./compactors_pipeline

bash 01.run_compactors.sh \
  --splash-dir /path/to/splash \
  --result-dir /path/to/splash_results
```

This reads:

```text
/path/to/splash_results/result_filtered_input/
/path/to/splash_results/result.after_correction.scores.tsv
```

and writes compactors output to:

```text
/path/to/splash_results/result_compactors_manual/
```

Then align and annotate compactors:

```bash
bash 02.align_annotate_assigned.sh \
  --input-tsv /path/to/result_compactors_manual/all_scores.compactors.tsv \
  --scores-tsv /path/to/splash_results/result.after_correction.scores.tsv \
  --star-genome-dir /path/to/star/genomeDir \
  --out-dir /path/to/result_compactors_manual/all_scores_annotated
```

The second script converts compactors TSV to FASTQ, aligns with STAR, extracts
assigned non-mitochondrial/non-DARLIN records, and writes:

```text
assigned_reads.txt
assigned_reads.txt.csv
```

## Environment

The scripts invoke tools through Pixi:

```bash
pixi run -e "$PIXI_ENV" ...
```

The default environment is `main`. Make sure it includes the Python packages and
external tools used by the requested steps, including SPLASH, Bowtie, STAR,
samtools, featureCounts/subread, and pybind11/C++ build tools where needed.

## C++ Alignment Acceleration

### This part is for the old version's sixth step and can be deleted.

`fast_align_cpp/` contains the local Python package `scan_splash_align`, a
pybind11/C++ backend used by alignment-heavy Python code such as
`main_pipeline/python/06.compare_targets_indel.py`.

Install/update the Pixi environment from the parent project if the local package
is declared there:

```bash
cd ./SCAN_SPLASH
pixi install -e main
```

During development, rebuild the extension from source if needed:

```bash
cd ./fast_align_cpp
pixi run -e main python -m pip install -e .
```

At runtime, the Python scripts try the C++ backend first and fall back to the
pure-Python implementation if the extension is unavailable.

## External References

- SPLASH: https://github.com/refresh-bio/SPLASH
- https://www.nature.com/articles/s41587-026-03084-6
- https://www.nature.com/articles/s41587-024-02381-2
- https://www.cell.com/cell/fulltext/S0092-8674(23)01179-0
