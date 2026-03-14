#!/usr/bin/env bash
# =============================================================================
# m02_dedup_filter.sh — Step 02: fixmate, markdup, quality filter (SLURM array)
# =============================================================================

set -euo pipefail

REFERENCE=""; INPUT_DIR=""; OUTPUT_DIR=""; SAMPLE_LIST=""
THREADS=8; MIN_MAPQ=20; CONDA_ENV=""; KEEP_INTERMEDIATES=false

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":r:i:o:s:t:q:e:k" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"
[[ -z "$SAMPLE_LIST" ]] && err "Sample list required (-s)"

# ── Activate conda environment ────────────────────────────────────────────────
[[ -z "$CONDA_ENV" ]] && err "Conda environment path required (-e)"
set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u
src=$PWD
log "Conda env: ${BOLD}${CONDA_ENV}${RESET}"

# ── Get sample for this array task ────────────────────────────────────────────
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(sed -n "${TASK_ID}p" "$SAMPLE_LIST")
[[ -z "$SAMPLE" ]] && err "No sample at line ${TASK_ID} in ${SAMPLE_LIST}"
log "Array task ${TASK_ID} → Sample: ${BOLD}${SAMPLE}${RESET}"

BAM_SORTED="${OUTPUT_DIR}/bam/${SAMPLE}_sorted.bam"
BAM_NAMESORT="${OUTPUT_DIR}/bam/${SAMPLE}_namesort.bam"
BAM_FIXMATE="${OUTPUT_DIR}/bam/${SAMPLE}_fixmate.bam"
BAM_COORDSORT="${OUTPUT_DIR}/bam/${SAMPLE}_coordsort.bam"
BAM_DEDUP="${OUTPUT_DIR}/bam/${SAMPLE}_dedup.bam"
BAM_FINAL="${OUTPUT_DIR}/bam/${SAMPLE}_final.bam"
STATS="${OUTPUT_DIR}/qc/${SAMPLE}_markdup_stats.txt"

[[ ! -f "$BAM_SORTED" ]] && err "Sorted BAM not found: $BAM_SORTED"
mkdir -p "${OUTPUT_DIR}/qc"

# ── Name sort (required for fixmate) ─────────────────────────────────────────
log "Sorting by queryname..."
samtools sort \
    -n \
    -@ "$THREADS" \
    -o "$BAM_NAMESORT" \
    "$BAM_SORTED"

# ── Fixmate (adds ms score tag required by markdup) ───────────────────────────
log "Running fixmate..."
samtools fixmate \
    -m \
    -@ "$THREADS" \
    "$BAM_NAMESORT" \
    "$BAM_FIXMATE"
rm -f "$BAM_NAMESORT"

# ── Re-sort by coordinate ─────────────────────────────────────────────────────
log "Re-sorting by coordinate..."
samtools sort \
    -@ "$THREADS" \
    -o "$BAM_COORDSORT" \
    "$BAM_FIXMATE"
rm -f "$BAM_FIXMATE"

# ── Mark and remove duplicates ────────────────────────────────────────────────
log "Marking duplicates..."
samtools markdup \
    -r \
    -@ "$THREADS" \
    -f "$STATS" \
    "$BAM_COORDSORT" \
    "$BAM_DEDUP"
samtools index "$BAM_DEDUP"
rm -f "$BAM_COORDSORT"

# ── Filter: mapped, primary, MAPQ threshold ───────────────────────────────────
log "Filtering (MAPQ >= ${MIN_MAPQ})..."
samtools view \
    -b \
    -F 4 \
    -F 256 \
    -F 2048 \
    -q "$MIN_MAPQ" \
    -@ "$THREADS" \
    -o "$BAM_FINAL" \
    "$BAM_DEDUP"
samtools index "$BAM_FINAL"

[[ "$KEEP_INTERMEDIATES" == false ]] && rm -f "$BAM_DEDUP" "${BAM_DEDUP}.bai"

# ── Coverage summary ──────────────────────────────────────────────────────────
log "Computing coverage..."
samtools coverage "$BAM_FINAL" \
    > "${OUTPUT_DIR}/qc/${SAMPLE}_coverage.txt" 2>/dev/null || true

# ── Stats ─────────────────────────────────────────────────────────────────────
MAPPED=$(samtools view -c -F 4 "$BAM_FINAL" 2>/dev/null || true)
TOTAL=$(samtools view -c "$BAM_SORTED" 2>/dev/null || true)
log "Total reads: ${TOTAL} | Mapped after filtering: ${MAPPED}"
[[ "$MAPPED" -eq 0 ]] && warn "0 reads remain after filtering for ${SAMPLE}"

ok "Dedup + filter complete → ${BAM_FINAL}"