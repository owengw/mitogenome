#!/usr/bin/env bash
# =============================================================================
# m05_align.sh — Step 05: MAFFT alignment + majority-rule consensus
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""; CONDA_ENV=""; THREADS=8

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:t:r:i:s:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
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

FASTA="${OUTPUT_DIR}/all_samples_mito.fa"
ALIGNED="${OUTPUT_DIR}/aligned/all_samples_mito_aligned.fa"
CONSENSUS="${OUTPUT_DIR}/aligned/consensus_mito.fa"

[[ ! -f "$FASTA" ]] && err "Input FASTA not found: $FASTA"
mkdir -p "${OUTPUT_DIR}/aligned"

N_SEQS=$(grep -c '>' "$FASTA")
log "Aligning ${N_SEQS} sequences with MAFFT (${THREADS} threads)..."

mafft \
    --auto \
    --thread "$THREADS" \
    --reorder \
    "$FASTA" \
    > "$ALIGNED"

ok "Alignment complete → ${ALIGNED}"

log "Generating majority-rule consensus..."

python3 - "$ALIGNED" "$CONSENSUS" <<'EOF'
import sys
from collections import Counter

aligned_path = sys.argv[1]
out_path = sys.argv[2]

seqs = []
current = []
with open(aligned_path) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
            current = []
        else:
            current.append(line)
    if current:
        seqs.append("".join(current))

length = len(seqs[0])
consensus = []
for i in range(length):
    col = [s[i] for s in seqs if i < len(s)]
    counts = Counter(b for b in col if b not in "-N")
    consensus.append(counts.most_common(1)[0][0] if counts else "N")

seq = "".join(consensus)
with open(out_path, "w") as f:
    f.write(">Plestiodon_capito_consensus_mitogenome\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")

print(f"Consensus: {len(seq)} bp")
EOF

ok "Consensus sequence → ${CONSENSUS}"