#!/usr/bin/env bash
# =============================================================================
# m13_dloop.sh — D-loop extraction, NJ tree, and haplotype network
#
# Steps:
#   1. Extract D-loop from each sample's assembled mitogenome using MitoZ GFF
#   2. Align with MAFFT
#   3. Build unrooted NJ tree (IQ-TREE)
#   4. Generate TCS haplotype network (Python/networkx) with population colours
#   5. Write all figures + HTML report (picked up by m11_report)
# =============================================================================
#SBATCH --job-name=m13_dloop
#SBATCH --output=/mnt/parscratch/users/bi4og/genome/mito_out/dloop/m13_%j.log
#SBATCH --error=/mnt/parscratch/users/bi4og/genome/mito_out/dloop/m13_%j.err
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

set -euo pipefail

# =============================================================================
# USER CONFIGURATION
# =============================================================================
PIPELINE_ENV="/users/bi4og/conda_envs/mito_pipeline"
OUTDIR="/mnt/parscratch/users/bi4og/genome/mito_out/dloop"
SCRIPT_DIR="/mnt/parscratch/users/bi4og/scripts"
MITO_OUT="/mnt/parscratch/users/bi4og/genome/mito_out"

# MitoZ annotation output directory (m08 output)
ANNOT_DIR="${MITO_OUT}/annotation"

# Assembled mitogenome FASTA directory (one .fa per sample from m04/m05)
ASSEMBLY_DIR="${MITO_OUT}/aligned"

# Population assignments: pattern matched against sample name
# Format: "PATTERN:POPULATION_CODE:HEX_COLOUR"
POP_DEFS=(
    "SI:SI:#3b82f6"    # blue
    "NS:NS:#22c55e"    # green
    "CAI:CAI:#f97316"  # orange
    "SB:SB:#a855f7"    # purple
)

THREADS=${SLURM_CPUS_PER_TASK:-4}
# =============================================================================

DLOOP_EXTRACT="${SCRIPT_DIR}/extract_dloop.py"
NETWORK_SCRIPT="${SCRIPT_DIR}/make_dloop_report.py"

export PATH="${PIPELINE_ENV}/bin:${PATH}"
export CONDA_PREFIX="${PIPELINE_ENV}"

PYTHON="${PIPELINE_ENV}/bin/python3"
MAFFT="${PIPELINE_ENV}/bin/mafft"
IQTREE="${PIPELINE_ENV}/bin/iqtree"

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

mkdir -p "$OUTDIR"

[[ ! -x "$MAFFT"  ]] && err "mafft not found: ${MAFFT}"
[[ ! -x "$IQTREE" ]] && err "iqtree not found: ${IQTREE}"

# ── Step 1: Extract D-loop from each sample ───────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 1: Extract D-loop sequences from MitoZ GFF${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

DLOOP_RAW="${OUTDIR}/dloop_all_raw.fa"
> "$DLOOP_RAW"   # empty/create

"$PYTHON" "$DLOOP_EXTRACT" \
    --annot-dir  "$ANNOT_DIR" \
    --assembly-dir "$ASSEMBLY_DIR" \
    --outfile    "$DLOOP_RAW" \
    --pop-defs   "${POP_DEFS[@]}"

N_DLOOP=$(grep -c "^>" "$DLOOP_RAW" 2>/dev/null || echo 0)
[[ "$N_DLOOP" -eq 0 ]] && err "No D-loop sequences extracted — check annotation directory"
ok "Extracted ${N_DLOOP} D-loop sequences"

# ── Step 2: Align ─────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 2: Align D-loop sequences with MAFFT${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

DLOOP_ALIGNED="${OUTDIR}/dloop_aligned.fa"
log "Running MAFFT..."
"$MAFFT" --auto --thread "$THREADS" --reorder "$DLOOP_RAW" > "$DLOOP_ALIGNED"

ALN_LEN=$(grep -v "^>" "$DLOOP_ALIGNED" | head -1 | tr -d '\n' | wc -c)
ok "Aligned: ${N_DLOOP} sequences, ${ALN_LEN} bp"

# ── Step 3: Unrooted NJ tree ─────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 3: Unrooted NJ tree (IQ-TREE)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

log "Running IQ-TREE NJ tree on D-loop alignment..."
"$IQTREE" \
    -s "$DLOOP_ALIGNED" \
    -m GTR+G \
    -B 1000 \
    -T "$THREADS" \
    --prefix "${OUTDIR}/dloop_nj" \
    --redo

ok "D-loop NJ tree → ${OUTDIR}/dloop_nj.treefile"

# ── Step 4: Haplotype network + HTML report ───────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 4: Haplotype network + HTML report${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

DLOOP_REPORT="${OUTDIR}/dloop_report.html"
[[ ! -f "$NETWORK_SCRIPT" ]] && err "make_dloop_report.py not found: ${NETWORK_SCRIPT}"

"$PYTHON" "$NETWORK_SCRIPT" \
    --aligned   "$DLOOP_ALIGNED" \
    --treefile  "${OUTDIR}/dloop_nj.treefile" \
    --output    "$DLOOP_REPORT" \
    --pop-defs  "${POP_DEFS[@]}"

ok "D-loop HTML report → ${DLOOP_REPORT}"

echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  OUTPUT SUMMARY${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "  Raw D-loop FASTA:    ${DLOOP_RAW}"
echo -e "  Aligned D-loop:      ${DLOOP_ALIGNED}"
echo -e "  NJ tree:             ${OUTDIR}/dloop_nj.treefile"
echo -e "  HTML report:         ${DLOOP_REPORT}"
echo ""
ok "m13_dloop complete"