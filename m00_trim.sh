#!/usr/bin/env bash
# =============================================================================
# 00_trim.sh — Step 00: Adapter trimming with fastp (SLURM array job)
# =============================================================================
# Called by run_pipeline.sh — do not run directly unless testing
# SLURM_ARRAY_TASK_ID maps to line number in sample list (1-indexed)
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REFERENCE=""; INPUT_DIR=""; OUTPUT_DIR=""; SAMPLE_LIST=""
PATTERN="R1R2_001"; THREADS=8; ADAPTER="auto"; CONDA_ENV=""
KEEP_INTERMEDIATES=false

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":r:i:o:s:p:t:a:ke:" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        p) PATTERN="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        a) ADAPTER="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        e) CONDA_ENV="$OPTARG" ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$SAMPLE_LIST" ]] && err "Sample list required (-s)"
[[ -z "$INPUT_DIR" ]]   && err "Input directory required (-i)"
[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"

# ── Resolve suffix from pattern ───────────────────────────────────────────────
case "$PATTERN" in
    R1R2_001) R1_SUFFIX="_R1_001.fastq.gz"; R2_SUFFIX="_R2_001.fastq.gz" ;;
    R1R2)     R1_SUFFIX="_R1.fastq.gz";     R2_SUFFIX="_R2.fastq.gz"     ;;
    12)       R1_SUFFIX="_1.fastq.gz";      R2_SUFFIX="_2.fastq.gz"      ;;
    *) err "Unknown pattern: $PATTERN" ;;
esac

# ── Adapter args ──────────────────────────────────────────────────────────────
case "$ADAPTER" in
    auto)    ADAPTER_ARGS="" ;;
    nextera) ADAPTER_ARGS="--adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT" ;;
    truseq)  ADAPTER_ARGS="--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" ;;
    none)    ADAPTER_ARGS="--disable_adapter_trimming" ;;
    *) err "Unknown adapter: $ADAPTER" ;;
esac

# ── Activate conda environment ────────────────────────────────────────────────
[[ -z "$CONDA_ENV" ]] && err "Conda environment path required (-e)"
set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u
src=$PWD
log "Conda env: ${BOLD}${CONDA_ENV}${RESET}"

# ── Get sample for this array task ───────────────────────────────────────────
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(sed -n "${TASK_ID}p" "$SAMPLE_LIST")
[[ -z "$SAMPLE" ]] && err "No sample at line ${TASK_ID} in ${SAMPLE_LIST}"

log "Array task ${TASK_ID} → Sample: ${BOLD}${SAMPLE}${RESET}"

R1="${INPUT_DIR}/${SAMPLE}${R1_SUFFIX}"
R2="${INPUT_DIR}/${SAMPLE}${R2_SUFFIX}"
[[ ! -f "$R1" ]] && err "R1 not found: $R1"
[[ ! -f "$R2" ]] && err "R2 not found: $R2"

TRIM_R1="${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
TRIM_R2="${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
FASTP_JSON="${OUTPUT_DIR}/qc/${SAMPLE}_fastp.json"
FASTP_HTML="${OUTPUT_DIR}/qc/${SAMPLE}_fastp.html"

mkdir -p "${OUTPUT_DIR}/trimmed" "${OUTPUT_DIR}/qc"

log "Trimming with fastp..."
fastp \
    --in1 "$R1" --in2 "$R2" \
    --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
    $ADAPTER_ARGS \
    --length_required 36 \
    --cut_front --cut_tail \
    --cut_mean_quality 20 \
    --thread "$THREADS" \
    --json "$FASTP_JSON" \
    --html "$FASTP_HTML"

ok "Trimming complete → ${TRIM_R1}"