#!/usr/bin/env bash
# =============================================================================
# m01_map.sh — Step 01: bwa mem mapping + samtools sort (SLURM array job)
# =============================================================================

set -euo pipefail

REFERENCE=""; INPUT_DIR=""; OUTPUT_DIR=""; SAMPLE_LIST=""
THREADS=8; CONDA_ENV=""; KEEP_INTERMEDIATES=false
NO_MERGE=false
FWD_PATTERN="*_f_paired.fastq.gz"
REV_PATTERN="*_r_paired.fastq.gz"

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":r:i:o:s:t:e:F:V:kM" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        F) FWD_PATTERN="$OPTARG" ;;
        V) REV_PATTERN="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        M) NO_MERGE=true ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$REFERENCE" ]]   && err "Reference required (-r)"
[[ -z "$INPUT_DIR" ]]   && err "Input directory required (-i)"
[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"
[[ -z "$SAMPLE_LIST" ]] && err "Sample list required (-s)"
[[ ! -f "$REFERENCE" ]] && err "Reference not found: $REFERENCE"

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

# ── Find all lane files for this sample ───────────────────────────────────────
mapfile -t F_FILES < <(find "$INPUT_DIR" -maxdepth 1 \
    -name "${SAMPLE}${FWD_PATTERN#\*}" | sort)
mapfile -t R_FILES < <(find "$INPUT_DIR" -maxdepth 1 \
    -name "${SAMPLE}${REV_PATTERN#\*}" | sort)

[[ ${#F_FILES[@]} -eq 0 ]] && err "No forward FASTQ files found for sample: ${SAMPLE}"
[[ ${#R_FILES[@]} -eq 0 ]] && err "No reverse FASTQ files found for sample: ${SAMPLE}"
[[ ${#F_FILES[@]} -ne ${#R_FILES[@]} ]] && \
    err "Mismatched lane counts for ${SAMPLE}: ${#F_FILES[@]} forward, ${#R_FILES[@]} reverse"

N_LANES=${#F_FILES[@]}
LANE_LABELS=$(for f in "${F_FILES[@]}"; do basename "$f"; done | tr '\n' ' ')
log "Lanes found (${N_LANES}): ${CYAN}${LANE_LABELS}${RESET}"

BAM_SORTED="${OUTPUT_DIR}/bam/${SAMPLE}_sorted.bam"
mkdir -p "${OUTPUT_DIR}/bam"

# ── Merge lanes into temp files if needed ─────────────────────────────────────
TMP_F="${OUTPUT_DIR}/bam/${SAMPLE}_tmp_R1.fastq.gz"
TMP_R="${OUTPUT_DIR}/bam/${SAMPLE}_tmp_R2.fastq.gz"

if [[ "$NO_MERGE" == true ]]; then
    # Single-lane mode: use first file pair directly, no concatenation
    [[ ${N_LANES} -gt 1 ]] && warn "${N_LANES} files found but --no-merge set; using only: ${F_FILES[0]}"
    F_INPUT="${F_FILES[0]}"
    R_INPUT="${R_FILES[0]}"
elif [[ ${#F_FILES[@]} -gt 1 ]]; then
    log "Merging ${N_LANES} lanes..."
    cat "${F_FILES[@]}" > "$TMP_F"
    cat "${R_FILES[@]}" > "$TMP_R"
    F_INPUT="$TMP_F"
    R_INPUT="$TMP_R"
else
    F_INPUT="${F_FILES[0]}"
    R_INPUT="${R_FILES[0]}"
fi

# ── Map ───────────────────────────────────────────────────────────────────────
log "Mapping with bwa mem (${THREADS} threads, ${N_LANES} lane(s))..."
bwa mem \
    -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_lib1" \
    "$REFERENCE" \
    "$F_INPUT" \
    "$R_INPUT" \
    | samtools sort \
        -@ "$THREADS" \
        -m 2G \
        -o "$BAM_SORTED" \
        -T "${OUTPUT_DIR}/bam/${SAMPLE}_sort_tmp"

samtools index "$BAM_SORTED"

# ── Clean up temp merged files ────────────────────────────────────────────────
[[ -f "$TMP_F" ]] && rm -f "$TMP_F" "$TMP_R"
log "Temp files cleaned up"

# ── Sanity check ──────────────────────────────────────────────────────────────
MAPPED=$(samtools view -c -F 4 "$BAM_SORTED" 2>/dev/null || true)
TOTAL=$(samtools view -c "$BAM_SORTED" 2>/dev/null || true)
log "Total reads: ${TOTAL} | Mapped: ${MAPPED}"
[[ "$MAPPED" -eq 0 ]] && warn "0 reads mapped for ${SAMPLE} — check reference or input files"

ok "Mapping complete → ${BAM_SORTED}"