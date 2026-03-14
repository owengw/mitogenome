#!/usr/bin/env bash
# =============================================================================
# 04_qc_summary.sh — Step 04: Aggregate QC summary (single SLURM job)
# =============================================================================
# Collects stats from all per-sample outputs and writes a TSV summary table.
# Also concatenates all consensus FASTAs into one multi-FASTA.
# =============================================================================

set -euo pipefail

REFERENCE=""; INPUT_DIR=""; OUTPUT_DIR=""; SAMPLE_LIST=""
THREADS=8; DEPTH_WARN=10

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":r:i:o:s:p:t:w:ke:" opt; do
    case $opt in
        r) REFERENCE="$OPTARG" ;;
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        s) SAMPLE_LIST="$OPTARG" ;;
        p) ;;   # accepted but unused here
        t) THREADS="$OPTARG" ;;
        w) DEPTH_WARN="$OPTARG" ;;
        k) ;;   # accepted but unused here
        e) CONDA_ENV="$OPTARG" ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]]  && err "Output directory required (-o)"
[[ -z "$SAMPLE_LIST" ]] && err "Sample list required (-s)"
[[ ! -f "$SAMPLE_LIST" ]] && err "Sample list not found: $SAMPLE_LIST"

SUMMARY="${OUTPUT_DIR}/qc/pipeline_summary.tsv"
COMBINED="${OUTPUT_DIR}/all_samples_mito.fa"

# ── Activate conda environment ────────────────────────────────────────────────
[[ -z "$CONDA_ENV" ]] && err "Conda environment path required (-e)"
set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u
src=$PWD
log "Conda env: ${BOLD}${CONDA_ENV}${RESET}"

# ── Header ────────────────────────────────────────────────────────────────────
echo -e "sample\tmapped_reads\tmean_depth\tmedian_depth\tcoverage_breadth_pct\tnum_variants_pass\tconsensus_length_bp\tN_count\tpct_N\tstatus" \
    > "$SUMMARY"

PASS=0; WARN=0; FAIL=0

# ── Per-sample ────────────────────────────────────────────────────────────────
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    BAM_FINAL="${OUTPUT_DIR}/bam/${SAMPLE}_final.bam"
    VCF_FILT="${OUTPUT_DIR}/vcf/${SAMPLE}_filtered.vcf.gz"
    CONSENSUS="${OUTPUT_DIR}/consensus/${SAMPLE}_mito.fa"
    COV_FILE="${OUTPUT_DIR}/qc/${SAMPLE}_coverage.txt"

    log "Summarising: ${SAMPLE}"

    # Mapped reads
    if [[ -f "$BAM_FINAL" ]]; then
        MAPPED_READS=$(samtools view -c -F 4 "$BAM_FINAL" 2>/dev/null || echo "NA")
    else
        MAPPED_READS="NA"
        warn "BAM not found for ${SAMPLE}"
    fi

    # Depth stats
    if [[ -f "$BAM_FINAL" ]]; then
        DEPTH_DATA=$(samtools depth "$BAM_FINAL" 2>/dev/null || true)
        if [[ -n "$DEPTH_DATA" ]]; then
            MEAN_DEPTH=$(echo "$DEPTH_DATA" \
                | awk '{sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n; else print "0"}')
            MEDIAN_DEPTH=$(echo "$DEPTH_DATA" \
                | awk '{print $3}' | sort -n \
                | awk 'BEGIN{n=0} {a[n++]=$1}
                       END{if(n%2==1) print a[int(n/2)];
                           else printf "%.1f",(a[n/2-1]+a[n/2])/2}')
        else
            MEAN_DEPTH="0"; MEDIAN_DEPTH="0"
        fi
    else
        MEAN_DEPTH="NA"; MEDIAN_DEPTH="NA"
    fi

    # Coverage breadth from pre-computed file or recompute
    if [[ -f "$COV_FILE" ]]; then
        COV_PCT=$(awk 'NR>1 {print $6}' "$COV_FILE" | head -1 || echo "NA")
    elif [[ -f "$BAM_FINAL" ]]; then
        COV_PCT=$(samtools coverage "$BAM_FINAL" 2>/dev/null \
            | awk 'NR>1 {print $6}' | head -1 || echo "NA")
    else
        COV_PCT="NA"
    fi

    # PASS variants
    if [[ -f "$VCF_FILT" ]]; then
        NUM_VARIANTS=$(bcftools view -f PASS "$VCF_FILT" 2>/dev/null \
            | grep -vc "^#" || echo "0")
    else
        NUM_VARIANTS="NA"
    fi

    # Consensus stats
    if [[ -f "$CONSENSUS" ]]; then
        CONS_LEN=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | wc -c)
        N_COUNT=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | tr -cd 'Nn' | wc -c)
        if [[ "$CONS_LEN" -gt 0 ]]; then
            PCT_N=$(python3 -c "print(round(100*${N_COUNT}/${CONS_LEN},2))")
        else
            PCT_N="NA"
        fi
    else
        CONS_LEN="NA"; N_COUNT="NA"; PCT_N="NA"
    fi

    # ── Status ────────────────────────────────────────────────────────────────
    STATUS="PASS"
    if [[ "$CONS_LEN" == "NA" || "$CONS_LEN" == "0" ]]; then
        STATUS="FAIL:no_consensus"
        (( FAIL++ )) || true
    elif [[ "$MEAN_DEPTH" != "NA" ]] && (( $(echo "$MEAN_DEPTH < $DEPTH_WARN" | bc -l) )); then
        STATUS="WARN:low_depth"
        (( WARN++ )) || true
    elif [[ "$PCT_N" != "NA" ]] && (( $(echo "$PCT_N > 20" | bc -l) )); then
        STATUS="WARN:high_N"
        (( WARN++ )) || true
    else
        (( PASS++ )) || true
    fi

    echo -e "${SAMPLE}\t${MAPPED_READS}\t${MEAN_DEPTH}\t${MEDIAN_DEPTH}\t${COV_PCT}\t${NUM_VARIANTS}\t${CONS_LEN}\t${N_COUNT}\t${PCT_N}\t${STATUS}" \
        >> "$SUMMARY"

done < "$SAMPLE_LIST"

# ── Combine consensus FASTAs ──────────────────────────────────────────────────
log "Combining consensus FASTAs..."
> "$COMBINED"
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue
    CONSENSUS="${OUTPUT_DIR}/consensus/${SAMPLE}_mito.fa"
    if [[ -f "$CONSENSUS" ]]; then
        cat "$CONSENSUS" >> "$COMBINED"
    else
        warn "No consensus for ${SAMPLE} — skipping from combined FASTA"
    fi
done < "$SAMPLE_LIST"

# ── Print summary ─────────────────────────────────────────────────────────────
TOTAL=$(wc -l < "$SAMPLE_LIST")
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  QC SUMMARY  (${TOTAL} samples)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "  ${GREEN}PASS:${RESET}  ${PASS}"
echo -e "  ${YELLOW}WARN:${RESET}  ${WARN}"
echo -e "  ${RED}FAIL:${RESET}  ${FAIL}"
echo ""
echo -e "  Full table → ${SUMMARY}"
echo -e "  Combined FASTA → ${COMBINED}"
echo ""
column -t -s $'\t' "$SUMMARY"
echo ""

ok "QC summary complete"