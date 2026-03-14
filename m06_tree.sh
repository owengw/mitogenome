#!/usr/bin/env bash
# =============================================================================
# m06_tree.sh — Step 06: IQ-TREE phylogenetic tree
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""; CONDA_ENV=""; THREADS=8; SPECIES="Unknown species"

RED='[0;31m'; GREEN='[0;32m'; CYAN='[0;36m'; BOLD='[1m'; RESET='[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:n:t:r:i:s:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        n) SPECIES="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        r|i|s|q|d|f|m|k) : ;;  # silently ignore unused common args
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$CONDA_ENV" ]]  && err "Conda environment required (-e)"

set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
set -u

ALIGNED="${OUTPUT_DIR}/aligned/all_samples_mito_aligned.fa"
TREE_DIR="${OUTPUT_DIR}/tree"

[[ ! -f "$ALIGNED" ]] && err "Aligned FASTA not found: $ALIGNED"
mkdir -p "$TREE_DIR"

log "Running IQ-TREE (species: ${SPECIES}, model: GTR+G, bootstrap: 1000, threads: ${THREADS})..."

iqtree \
    -s "$ALIGNED" \
    --prefix "${TREE_DIR}/mito_tree" \
    -m GTR+G \
    -B 1000 \
    -T "$THREADS" \
    --redo

ok "Tree complete → ${TREE_DIR}/mito_tree.treefile"