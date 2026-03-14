#!/usr/bin/env bash
# =============================================================================
# m08_annotate.sh — Step 07: Extract per-gene alignments from individual mitogenomes
# =============================================================================
# Uses GBF coordinates from MitoZ annotation of the consensus to extract
# each CDS from all individual mitogenomes, then aligns with MAFFT.
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""
ANNOT_ENV=""
GENETIC_CODE=2
CLADE="Chordata"

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[1;33m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:C:c:t:r:i:s:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) ANNOT_ENV="$OPTARG" ;;
        C) GENETIC_CODE="$OPTARG" ;;
        c) CLADE="$OPTARG" ;;
        t|r|i|s|q|d|f|m|k) : ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$ANNOT_ENV"  ]] && err "Annotation conda environment required (-e)"

set +u
source ~/.bash_profile
conda activate "$ANNOT_ENV"
set -u

GBF="${OUTPUT_DIR}/annotation_test/test.consensus_mito.fa.result/test_consensus_mito.fa_mitoscaf.fa.gbf"
CONSENSUS_DIR="${OUTPUT_DIR}/consensus"
DNDS_DIR="${OUTPUT_DIR}/dnds"
ALN_DIR="${DNDS_DIR}/alignments"

[[ ! -f "$GBF" ]]          && err "GBF annotation file not found: $GBF"
[[ ! -d "$CONSENSUS_DIR" ]] && err "Consensus directory not found: $CONSENSUS_DIR"
mkdir -p "$ALN_DIR"

log "Genetic code: ${GENETIC_CODE} | MitoZ clade: ${CLADE}"
log "Parsing gene coordinates from GBF..."

python3 - "$GBF" "$CONSENSUS_DIR" "$ALN_DIR" "$GENETIC_CODE" <<'EOF'
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

gbf_path     = sys.argv[1]
cons_dir     = sys.argv[2]
aln_dir      = sys.argv[3]
genetic_code = int(sys.argv[4])

# ── Parse CDS coordinates from GBF ───────────────────────────────────────────
genes = {}
with open(gbf_path) as f:
    current_gene = None
    for line in f:
        m_cds  = re.match(r'\s+CDS\s+(?:complement\()?(\d+)\.\.(\d+)\)?', line)
        m_gene = re.match(r'\s+/gene="([^"]+)"', line)
        if m_cds:
            start  = int(m_cds.group(1)) - 1  # 0-based
            end    = int(m_cds.group(2))
            strand = -1 if 'complement' in line else 1
            current_gene = {'start': start, 'end': end, 'strand': strand}
        elif m_gene and current_gene is not None:
            gene_name = m_gene.group(1)
            if gene_name not in genes:
                genes[gene_name] = current_gene
            current_gene = None

print(f"Genes parsed: {list(genes.keys())}")

# ── Load all individual consensus FASTAs ─────────────────────────────────────
fa_files = sorted([f for f in os.listdir(cons_dir) if f.endswith('_mito.fa')])
print(f"Individual FASTAs found: {len(fa_files)}")

seqs = {}
for fa in fa_files:
    sample = fa.replace('_mito.fa', '')
    rec = next(SeqIO.parse(os.path.join(cons_dir, fa), 'fasta'))
    seqs[sample] = str(rec.seq)

# ── Extract per-gene sequences ────────────────────────────────────────────────
for gene, coords in genes.items():
    start  = coords['start']
    end    = coords['end']
    strand = coords['strand']

    records = []
    for sample, full_seq in seqs.items():
        subseq = full_seq[start:end]
        if strand == -1:
            subseq = str(Seq(subseq).reverse_complement())
        records.append(SeqRecord(Seq(subseq), id=sample, description=''))

    out_fa = os.path.join(aln_dir, f'{gene}_unaligned.fa')
    SeqIO.write(records, out_fa, 'fasta')
    print(f"  {gene}: {end - start} bp, {len(records)} sequences written")

print("Gene extraction complete.")
EOF

ok "Gene sequences extracted → ${ALN_DIR}"

# ── Align each gene with MAFFT ────────────────────────────────────────────────
log "Aligning each gene with MAFFT..."
for UNALIGNED in "${ALN_DIR}"/*_unaligned.fa; do
    GENE=$(basename "$UNALIGNED" _unaligned.fa)
    ALIGNED="${ALN_DIR}/${GENE}_aligned.fa"
    log "  Aligning ${GENE}..."
    mafft --auto --thread 8 --reorder "$UNALIGNED" > "$ALIGNED" 2>/dev/null
    rm -f "$UNALIGNED"
    ok "  ${GENE} aligned → ${ALIGNED}"
done

ok "All gene alignments complete → ${ALN_DIR}"