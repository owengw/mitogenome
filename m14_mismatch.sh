#!/usr/bin/env bash
# =============================================================================
# m14_mismatch.sh — Mismatch distribution analysis per population
#
# Uses D-loop aligned sequences from m13.
# Computes observed pairwise mismatch distributions per population,
# fits Rogers & Harpending (1992) sudden expansion model,
# outputs JSON summary for inclusion in m11 HTML report.
# =============================================================================
#SBATCH --job-name=m14_mismatch
#SBATCH --output=/mnt/parscratch/users/bi4og/genome/mito_out/mismatch/m14_%j.log
#SBATCH --error=/mnt/parscratch/users/bi4og/genome/mito_out/mismatch/m14_%j.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00

set -euo pipefail

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PIPELINE_ENV="/users/bi4og/conda_envs/mito_pipeline"
DLOOP_ALIGNED="/mnt/parscratch/users/bi4og/genome/mito_out/dloop/dloop_aligned.fa"
OUTDIR="/mnt/parscratch/users/bi4og/genome/mito_out/mismatch"
SCRIPT_DIR="/mnt/parscratch/users/bi4og/scripts"

POP_DEFS=(
    "SI:SI:#3b82f6"
    "NS:NS:#22c55e"
    "CAI:CAI:#f97316"
    "SB:SB:#a855f7"
)
# =============================================================================

export PATH="${PIPELINE_ENV}/bin:${PATH}"
PYTHON="${PIPELINE_ENV}/bin/python3"

GREEN='\033[0;32m'; CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log() { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
err() { echo -e "\033[0;31m[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

mkdir -p "$OUTDIR"

[[ ! -f "$DLOOP_ALIGNED" ]] && err "D-loop alignment not found: ${DLOOP_ALIGNED}. Run m13 first."

log "Running mismatch analysis..."

"$PYTHON" "${SCRIPT_DIR}/compute_mismatch.py" \
    --aligned  "$DLOOP_ALIGNED" \
    --outdir   "$OUTDIR" \
    --pop-defs "${POP_DEFS[@]}"

ok "Mismatch analysis complete → ${OUTDIR}/mismatch_summary.json"