#!/usr/bin/env bash
# =============================================================================
# 03_call_consensus.sh — Step 03: bcftools variant call + consensus (SLURM array)
# =============================================================================

set -euo pipefail

REFERENCE=""; INPUT_DIR=""; OUTPUT_DIR=""; SAMPLE_LIST=""
PATTERN="R1R2_001"; THREADS=8
MIN_DEPTH=3; MIN_AF=0.7; MIN_BQ=20; CONDA_ENV=""
KEEP_INTERMEDIATES=false

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":r:i:o:s:p:t:d:f:m:ke:" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        p) PATTERN="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) MIN_DEPTH="$OPTARG" ;;
        f) MIN_AF="$OPTARG" ;;
        m) MIN_BQ="$OPTARG" ;;
        k) KEEP_INTERMEDIATES=true ;;
        e) CONDA_ENV="$OPTARG" ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$REFERENCE" ]]   && err "Reference required (-r)"
[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"
[[ -z "$SAMPLE_LIST" ]] && err "Sample list required (-s)"
[[ ! -f "$REFERENCE" ]] && err "Reference not found: $REFERENCE"

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
SAMPLE=$(sed -n "${TASK_ID}p" "$SAMPLE_LIST")
[[ -z "$SAMPLE" ]] && err "No sample at line ${TASK_ID} in ${SAMPLE_LIST}"

log "Array task ${TASK_ID} → Sample: ${BOLD}${SAMPLE}${RESET}"

BAM_FINAL="${OUTPUT_DIR}/bam/${SAMPLE}_final.bam"
VCF_RAW="${OUTPUT_DIR}/vcf/${SAMPLE}_raw.vcf.gz"
VCF_FILT="${OUTPUT_DIR}/vcf/${SAMPLE}_filtered.vcf.gz"
CONSENSUS="${OUTPUT_DIR}/consensus/${SAMPLE}_mito.fa"

[[ ! -f "$BAM_FINAL" ]] && err "Final BAM not found: $BAM_FINAL — did step 02 complete?"

mkdir -p "${OUTPUT_DIR}/vcf" "${OUTPUT_DIR}/consensus"

# ── Activate conda environment ────────────────────────────────────────────────
[[ -z "$CONDA_ENV" ]] && err "Conda environment path required (-e)"
set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"
set -u
src=$PWD
log "Conda env: ${BOLD}${CONDA_ENV}${RESET}"

# ── Warn if low coverage ──────────────────────────────────────────────────────
MEAN_DEPTH=$(samtools depth "$BAM_FINAL" \
    | awk '{sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}')
if (( $(echo "$MEAN_DEPTH < $MIN_DEPTH" | bc -l) )); then
    warn "Mean depth (${MEAN_DEPTH}×) below threshold (${MIN_DEPTH}×) — consensus may contain many Ns"
fi

# ── Variant calling ───────────────────────────────────────────────────────────
log "Calling variants (bcftools mpileup | call)..."
bcftools mpileup \
    -f "$REFERENCE" \
    -Q "$MIN_BQ" \
    -d 10000 \
    --annotate FORMAT/AD,FORMAT/DP \
    "$BAM_FINAL" \
    | bcftools call \
        -m \
        --ploidy 1 \
        -Oz \
        -o "$VCF_RAW"

tabix -p vcf "$VCF_RAW"

TOTAL_VARS=$(bcftools view "$VCF_RAW" | grep -vc "^#" || true)
log "Raw variants called: ${TOTAL_VARS}"

# ── Filter variants ───────────────────────────────────────────────────────────
log "Filtering variants (QUAL >= 20, DP >= ${MIN_DEPTH})..."
bcftools filter \
    -s 'LowQual' \
    -e "QUAL<20 || FORMAT/DP[0]<${MIN_DEPTH}" \
    -Oz \
    -o "$VCF_FILT" \
    "$VCF_RAW"

tabix -p vcf "$VCF_FILT"

PASS_VARS=$(bcftools view -f PASS "$VCF_FILT" | grep -vc "^#" || true)
log "PASS variants after filtering: ${PASS_VARS}"

# ── Generate consensus ────────────────────────────────────────────────────────
log "Generating consensus FASTA..."
bcftools consensus \
    -f "$REFERENCE" \
    --iupac-codes \
    -H 1 \
    "$VCF_FILT" \
    | sed "s/^>.*/>${SAMPLE}_mito/" \
    > "$CONSENSUS"

# ── Consensus stats ───────────────────────────────────────────────────────────
CONS_LEN=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | wc -c)
N_COUNT=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | tr -cd 'Nn' | wc -c)
PCT_N=$(python3 -c "print(round(100*${N_COUNT}/${CONS_LEN},2))" 2>/dev/null || echo "NA")
log "Consensus length: ${CONS_LEN} bp | N count: ${N_COUNT} (${PCT_N}%)"

if (( $(echo "$PCT_N > 20" | bc -l) )) 2>/dev/null; then
    warn "High N content (${PCT_N}%) — check coverage for this sample"
fi

# ── Cleanup ───────────────────────────────────────────────────────────────────
if [[ "$KEEP_INTERMEDIATES" == false ]]; then
    rm -f "$VCF_RAW" "${VCF_RAW}.tbi"
    log "Removed raw VCF"
fi

ok "Consensus complete → ${CONSENSUS}"