#!/usr/bin/env bash
# =============================================================================
# m11_report.sh — Step 10: Generate HTML pipeline report
# =============================================================================
# Runs mito_report.py to produce two HTML reports:
#   - mito_report_all.html      (all populations, including small ones)
#   - mito_report_filtered.html (populations with n >= MIN_POP_SIZE only)
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""
CONDA_ENV=""
SPECIES="Unknown species"
MIN_POP_SIZE=5
# SCRIPT_DIR passed via -S at submission time from mito_pipeline.sh, which
# resolves the real path before sbatch copies the script to the SLURM spool.
# Falls back to dirname of this script for manual invocation.
SCRIPT_DIR_OVERRIDE=""

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[1;33m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:s:p:S:t:r:i:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        s) SPECIES="$OPTARG" ;;
        p) MIN_POP_SIZE="$OPTARG" ;;
        S) SCRIPT_DIR_OVERRIDE="$OPTARG" ;;
        t|r|i|q|d|f|m|k) : ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$CONDA_ENV"  ]] && err "Conda environment required (-e)"

# Use -S path if provided (SLURM-safe); otherwise resolve from script location
if [[ -n "$SCRIPT_DIR_OVERRIDE" ]]; then
    SCRIPT_DIR="$SCRIPT_DIR_OVERRIDE"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

REPORT_SCRIPT="${SCRIPT_DIR}/mito_report.py"
[[ ! -f "$REPORT_SCRIPT" ]] && err "mito_report.py not found: ${REPORT_SCRIPT}"

set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u

log "Generating reports for: ${SPECIES}"
log "Min population size filter: ${MIN_POP_SIZE}"
log "mito_report.py: ${REPORT_SCRIPT}"

# Full report (all populations)
log "Generating full report (all populations)..."
python3 "$REPORT_SCRIPT" \
    --outdir        "$OUTPUT_DIR" \
    --species       "$SPECIES" \
    --min-pop-size  0 \
    --output        "${OUTPUT_DIR}/mito_report_all.html" \
    && ok "Full report → ${OUTPUT_DIR}/mito_report_all.html" \
    || warn "Full report generation failed"

# Filtered report (populations >= MIN_POP_SIZE)
log "Generating filtered report (n >= ${MIN_POP_SIZE} per population)..."
python3 "$REPORT_SCRIPT" \
    --outdir        "$OUTPUT_DIR" \
    --species       "$SPECIES" \
    --min-pop-size  "$MIN_POP_SIZE" \
    --output        "${OUTPUT_DIR}/mito_report_filtered.html" \
    && ok "Filtered report → ${OUTPUT_DIR}/mito_report_filtered.html" \
    || warn "Filtered report generation failed"

ok "Report generation complete"
echo ""
echo -e "  Full report:     ${CYAN}${OUTPUT_DIR}/mito_report_all.html${RESET}"
echo -e "  Filtered report: ${CYAN}${OUTPUT_DIR}/mito_report_filtered.html${RESET}"