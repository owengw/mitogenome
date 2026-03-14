#!/usr/bin/env bash
# =============================================================================
# run_cytb_phylo.sh — Download CYTB sequences, combine with P. longirostris
#                     consensus, align with MAFFT, and build IQ-TREE phylogeny
# =============================================================================
#SBATCH --job-name=cytb_phylo
#SBATCH --output=/mnt/parscratch/users/bi4og/genome/mito_out/cytb_phylo/cytb_phylo_%j.log
#SBATCH --error=/mnt/parscratch/users/bi4og/genome/mito_out/cytb_phylo/cytb_phylo_%j.err
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

set -euo pipefail

# =============================================================================
# USER CONFIGURATION — edit these before submitting
# =============================================================================
EMAIL="bi4og@sheffield.ac.uk"
CONDA_ENV="/users/bi4og/conda_envs/mito_pipeline"
OUTDIR="/mnt/parscratch/users/bi4og/genome/mito_out/cytb_phylo"
SCRIPT_DIR="/mnt/parscratch/users/bi4og/scripts"

# Path to P. longirostris CYTB from m09 — leave blank to auto-detect
LONGIROSTRIS_CYTB=""

# Analysis parameters
MAX_PER_SPECIES=3
THREADS=${SLURM_CPUS_PER_TASK:-8}   # uses SLURM allocation, falls back to 8
BOOTSTRAP=1000
MIN_LEN=900
MAX_LEN=1200
SKIP_FETCH=false   # set true to skip GenBank download and reuse existing sequences
# ── Optional advanced features ───────────────────────────────────────────────

# Expand sampling across all skinks (Scincidae)
GLOBAL_SCINCIDAE_SAMPLING=false
MAX_GLOBAL_SCINCIDAE=200

# Reproducible subsampling of longirostris haplotypes
MAX_LONGIROSTRIS_HAPS=30
HAP_RANDOM_SEED=42

# Remove duplicate sequences across dataset
DEDUP_SEQUENCES=true

# Faster IQ-TREE bootstrap
USE_ULTRAFAST_BOOTSTRAP=true

# Generate publication-ready circular tree
GENERATE_CIRCULAR_TREE=true
# =============================================================================

FETCH_SCRIPT="${SCRIPT_DIR}/fetch_cytb.py"

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

# ── Help ──────────────────────────────────────────────────────────────────────
usage() {
cat << EOF
${BOLD}run_cytb_phylo.sh${RESET} — CYTB phylogenetic placement of Plestiodon longirostris

${BOLD}USAGE${RESET}
  # Submit to SLURM (recommended) — edit USER CONFIGURATION block first
  sbatch run_cytb_phylo.sh

  # Run interactively with argument overrides
  bash run_cytb_phylo.sh -e <conda_env> -o <outdir> --email <email> [OPTIONS]

${BOLD}OPTIONS (override USER CONFIGURATION block)${RESET}
  -e  Full path to conda environment (needs Biopython, MAFFT, IQ-TREE)
  -o  Output directory
  -l  Path to P. longirostris CYTB FASTA from m09
      (default: auto-detected from <outdir>/../dnds/alignments/CYTB_aligned.fa)
  -t  Threads for IQ-TREE (default: 8)
  -b  Bootstrap replicates (default: 1000)
  --email             Email address for NCBI Entrez
  --max-per-species   Max GenBank sequences per species (default: 3)
  --skip-fetch        Skip GenBank download, reuse existing sequences in outdir
  --fetch-script      Path to fetch_cytb.py (default: SCRIPT_DIR/fetch_cytb.py)

${BOLD}EXAMPLES${RESET}
  # Standard sbatch submission
  sbatch /mnt/parscratch/users/bi4og/scripts/run_cytb_phylo.sh

  # Re-run alignment + tree only (sequences already downloaded)
  sbatch /mnt/parscratch/users/bi4og/scripts/run_cytb_phylo.sh --skip-fetch

EOF
exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --email)           EMAIL="$2";           shift 2 ;;
        --max-per-species) MAX_PER_SPECIES="$2"; shift 2 ;;
        --skip-fetch)      SKIP_FETCH=true;      shift ;;
        --fetch-script)    FETCH_SCRIPT="$2";    shift 2 ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

while getopts ":e:o:l:t:b:h" opt; do
    case $opt in
        e) CONDA_ENV="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        l) LONGIROSTRIS_CYTB="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        b) BOOTSTRAP="$OPTARG" ;;
        h) usage ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
[[ -z "$CONDA_ENV" ]] && err "Conda environment required (-e)"
[[ -z "$CONDA_ENV" ]] && err "Conda environment not set — edit USER CONFIGURATION block"
[[ -z "$OUTDIR" ]]    && err "Output directory not set — edit USER CONFIGURATION block"
[[ -z "$EMAIL" && "$SKIP_FETCH" == false ]] && err "Email not set — edit USER CONFIGURATION block"
[[ ! -d "$CONDA_ENV" ]] && err "Conda environment not found: $CONDA_ENV"
[[ "$SKIP_FETCH" == false && ! -f "$FETCH_SCRIPT" ]] && \
    err "fetch_cytb.py not found: ${FETCH_SCRIPT} — check SCRIPT_DIR in USER CONFIGURATION"

# ── Activate conda (prepend bin dir directly for SLURM compatibility) ─────────
# conda activate is unreliable in SLURM as it requires shell init hooks.
# Prepending the env bin dir achieves the same effect reliably.
export PATH="${CONDA_ENV}/bin:${PATH}"
export CONDA_PREFIX="${CONDA_ENV}"
log "Conda env: ${BOLD}${CONDA_ENV}${RESET}"

# ── Tool paths (explicit, never rely on PATH alone) ───────────────────────────
MAFFT="${CONDA_ENV}/bin/mafft"
IQTREE="${CONDA_ENV}/bin/iqtree"
PYTHON="${CONDA_ENV}/bin/python3"
[[ ! -x "$MAFFT"  ]] && err "mafft not found: ${MAFFT}"
[[ ! -x "$IQTREE" ]] && err "iqtree not found: ${IQTREE}"
"$PYTHON" -c "from Bio import Entrez" 2>/dev/null || \
    err "Biopython not found in ${CONDA_ENV}"

mkdir -p "$OUTDIR"

# ── Key file paths ────────────────────────────────────────────────────────────
GENBANK_FA="${OUTDIR}/all_cytb_outgroups.fa"
COMBINED_FA="${OUTDIR}/all_cytb_combined.fa"
ALIGNED_FA="${OUTDIR}/all_cytb_aligned.fa"
TREE_PREFIX="${OUTDIR}/cytb_phylo"

# ── Auto-detect P. longirostris CYTB if not specified ─────────────────────────
if [[ -z "$LONGIROSTRIS_CYTB" ]]; then
    # Try standard m09 output location relative to outdir
    CANDIDATE="${OUTDIR}/../dnds/alignments/CYTB_aligned.fa"
    CANDIDATE="$(realpath "$CANDIDATE" 2>/dev/null || echo "$CANDIDATE")"
    if [[ -f "$CANDIDATE" ]]; then
        LONGIROSTRIS_CYTB="$CANDIDATE"
        log "Auto-detected P. longirostris CYTB: ${LONGIROSTRIS_CYTB}"
    else
        warn "P. longirostris CYTB not found at default location (${CANDIDATE})"
        warn "Proceeding without — tree will only contain GenBank sequences"
        warn "Specify with -l to include your consensus sequences"
    fi
fi

# ── Step 1: Fetch GenBank sequences ───────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 1: Download CYTB sequences from GenBank${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

if [[ "$SKIP_FETCH" == true ]]; then
    log "Skipping fetch (--skip-fetch)"
    [[ ! -f "$GENBANK_FA" ]] && err "Expected ${GENBANK_FA} not found — cannot skip fetch"
    N_GENBANK=$(grep -c "^>" "$GENBANK_FA" || true)
    log "Using existing sequences: ${N_GENBANK} sequences in ${GENBANK_FA}"
else
    log "Fetching sequences (this may take a few minutes)..."
    "$PYTHON" "$FETCH_SCRIPT" \
        --email "$EMAIL" \
        --outdir "$OUTDIR" \
        --max-per-species "$MAX_PER_SPECIES" \
        --min-len "$MIN_LEN" \
        --max-len "$MAX_LEN"
    [[ ! -f "$GENBANK_FA" ]] && err "Fetch failed — ${GENBANK_FA} not created"
    N_GENBANK=$(grep -c "^>" "$GENBANK_FA" || true)
    ok "Downloaded ${N_GENBANK} sequences → ${GENBANK_FA}"
fi

# ── Step 1b: Optional global Scincidae sampling ──────────────────────────────
if [[ "$GLOBAL_SCINCIDAE_SAMPLING" == true ]]; then

echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 1b: Expanding dataset across Scincidae${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

GLOBAL_SCINCIDAE_FA="${OUTDIR}/scincidae_global_sample.fa"

"$PYTHON" <<PYEOF
from Bio import Entrez, SeqIO
import random

Entrez.email="${EMAIL}"

query="Scincidae[Organism] AND cytb[Gene]"

handle=Entrez.esearch(db="nucleotide", term=query, retmax=2000)
ids=Entrez.read(handle)["IdList"]

handle=Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")

records=[r for r in SeqIO.parse(handle,"fasta") if ${MIN_LEN} <= len(r.seq) <= ${MAX_LEN}]

random.seed(${HAP_RANDOM_SEED})
records=random.sample(records, min(len(records), ${MAX_GLOBAL_SCINCIDAE}))

SeqIO.write(records,"${GLOBAL_SCINCIDAE_FA}","fasta")

print("Added",len(records),"Scincidae CYTB sequences")
PYEOF

cat "$GLOBAL_SCINCIDAE_FA" >> "$GENBANK_FA"

fi

# ── Step 2: Combine with P. longirostris consensus ────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 2: Combine with P. longirostris consensus CYTB${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

cp "$GENBANK_FA" "$COMBINED_FA"

if [[ -n "$LONGIROSTRIS_CYTB" && -f "$LONGIROSTRIS_CYTB" ]]; then
    # Rename sequences to make them identifiable in the tree
    # Input headers are like ">CYTB_SI75_mito" — prefix with species
    "$PYTHON" - "$LONGIROSTRIS_CYTB" "$COMBINED_FA" << 'PYEOF'
import sys
from Bio import SeqIO

cytb_fa   = sys.argv[1]
combined  = sys.argv[2]

records = list(SeqIO.parse(cytb_fa, "fasta"))
print(f"  Adding {len(records)} P. longirostris CYTB sequences")

MAX_HAPS = 30   # limit number of P. longirostris haplotypes used in tree

# Deduplicate identical sequences, keeping one per unique haplotype
seen = {}
haplotypes = []
for rec in records:
    seq_str = str(rec.seq).upper()
    if seq_str not in seen:
        seen[seq_str] = rec
        haplotypes.append(rec)

print(f"  Unique haplotypes: {len(haplotypes)}")

if len(haplotypes) > MAX_HAPS:
    print(f"  Subsampling to {MAX_HAPS} haplotypes")
    haplotypes = haplotypes[:MAX_HAPS]

# Rename for clarity in tree output
for i, rec in enumerate(haplotypes, 1):
    rec.id          = f"Plestiodon_longirostris_hap{i}"
    rec.name        = rec.id
    rec.description = ""

with open(combined, "a") as fh:
    SeqIO.write(haplotypes, fh, "fasta")

print(f"  Written to {combined}")
PYEOF
    N_COMBINED=$(grep -c "^>" "$COMBINED_FA" || true)
    ok "Combined FASTA: ${N_COMBINED} sequences → ${COMBINED_FA}"
else
    warn "No P. longirostris CYTB added — proceeding with GenBank sequences only"
    N_COMBINED=$(grep -c "^>" "$COMBINED_FA" || true)
fi

#log "Running BLAST sanity check on sampled sequences..."

#BLASTN="${CONDA_ENV}/bin/blastn"
#BLAST_DB="nt"

#SAMPLE_FA="${OUTDIR}/blast_check.fa"

#seqkit sample -n 5 "$COMBINED_FA" > "$SAMPLE_FA"

#"$BLASTN" \
#    -query "$SAMPLE_FA" \
#    -db "$BLAST_DB" \
#    -max_target_seqs 5 \
#    -outfmt "6 qseqid sscinames pident length evalue" \
#    > "${OUTDIR}/blast_check.tsv"

#log "BLAST check results:"
#cat "${OUTDIR}/blast_check.tsv"

# ── remove duplicate sequences ───────────────────────────────

if [[ "$DEDUP_SEQUENCES" == true ]]; then

echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  EXTRA STEP: Removing duplicate sequences${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

DEDUP_TMP="${OUTDIR}/combined_dedup_tmp.fa"

seqkit rmdup -s "$COMBINED_FA" > "$DEDUP_TMP"

mv "$DEDUP_TMP" "$COMBINED_FA"

N_DEDUP=$(grep -c "^>" "$COMBINED_FA" || true)

log "After deduplication: ${N_DEDUP} sequences"

fi

# ── Step 3: Align with MAFFT ──────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 3: Align with MAFFT${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

log "Running MAFFT (--auto, ${THREADS} threads)..."
"$MAFFT" \
    --auto \
    --thread "$THREADS" \
    --reorder \
    "$COMBINED_FA" \
    > "$ALIGNED_FA"

[[ ! -s "$ALIGNED_FA" ]] && err "MAFFT produced an empty alignment: $ALIGNED_FA"

N_ALIGNED=$(grep -c "^>" "$ALIGNED_FA" || true)

ALN_LEN=$(awk '
/^>/ {next}
{print length($0); exit}
' "$ALIGNED_FA" 2>/dev/null || echo 0)
ok "Alignment complete: ${N_ALIGNED} sequences, ${ALN_LEN} bp → ${ALIGNED_FA}"

TRIMAL="${CONDA_ENV}/bin/trimal"

TRIMMED_FA="${OUTDIR}/all_cytb_aligned_trimmed.fa"

log "Trimming alignment (trimAl automated1)..."
"$TRIMAL" -in "$ALIGNED_FA" -out "$TRIMMED_FA" -automated1

ALIGNED_FA="$TRIMMED_FA"

# ── Step 4: IQ-TREE ───────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 4: Build phylogeny with IQ-TREE${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

# ── choose bootstrap strategy ────────────────────────────────

if [[ "$USE_ULTRAFAST_BOOTSTRAP" == true ]]; then
BOOTSTRAP_ARGS="-bb ${BOOTSTRAP} -alrt 1000"
else
BOOTSTRAP_ARGS="-B ${BOOTSTRAP}"
fi

log "Running IQ-TREE (GTR+G, ${BOOTSTRAP} ultrafast bootstraps, ${BOOTSTRAP} SH-aLRT, ${THREADS} threads)..."
OUTGROUP_SEQ=$(grep "^>Eumeces_schneiderii" "$ALIGNED_FA" | head -n1 | tr -d ">")

"$IQTREE" \
    -s "$ALIGNED_FA" \
    -m MFP \
    -B "$BOOTSTRAP" \
    -alrt "$BOOTSTRAP" \
    -T "$THREADS" \
    -o "$OUTGROUP_SEQ,Scincella_lateralis_AY217806.1" \
    --prefix "$TREE_PREFIX" \
    --redo

ok "IQ-TREE complete → ${TREE_PREFIX}.treefile"

# ── Publication circular tree ────────────────────────────────────

if [[ "$GENERATE_CIRCULAR_TREE" == true ]]; then

CIRC_SCRIPT="${SCRIPT_DIR}/make_circular_tree.py"
CIRC_TREE="${OUTDIR}/cytb_tree_circular.svg"

log "Generating circular phylogenetic tree figure"

"$PYTHON" "$CIRC_SCRIPT" \
"${TREE_PREFIX}.treefile" \
"$CIRC_TREE"

fi

# ── Step 5: Generate HTML report ──────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 5: Generate HTML tree report${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

# ── Dataset summary table ────────────────────────────────────────

SUMMARY_TSV="${OUTDIR}/dataset_summary.tsv"

echo -e "metric\tvalue" > "$SUMMARY_TSV"
echo -e "genbank_sequences\t$(grep -c '^>' $GENBANK_FA)" >> "$SUMMARY_TSV"
echo -e "combined_sequences\t$(grep -c '^>' $COMBINED_FA)" >> "$SUMMARY_TSV"
echo -e "alignment_length\t${ALN_LEN}" >> "$SUMMARY_TSV"

log "Dataset summary written → ${SUMMARY_TSV}"

HTML_REPORT="${OUTDIR}/cytb_phylo_report.html"
REPORT_SCRIPT="${SCRIPT_DIR}/make_cytb_report.py"
[[ ! -f "$REPORT_SCRIPT" ]] && err "make_cytb_report.py not found: ${REPORT_SCRIPT} — deploy alongside this script"

log "Building HTML report → ${HTML_REPORT}"
"$PYTHON" "$REPORT_SCRIPT" "${TREE_PREFIX}.treefile" "$HTML_REPORT"
ok "HTML report → ${HTML_REPORT}"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  OUTPUT SUMMARY${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "  GenBank sequences:   ${GENBANK_FA}"
echo -e "  Combined FASTA:      ${COMBINED_FA}"
echo -e "  Aligned FASTA:       ${ALIGNED_FA}"
echo -e "  Tree file:           ${TREE_PREFIX}.treefile"
echo -e "  HTML report:         ${HTML_REPORT}"
echo -e "  Fetch summary:       ${OUTDIR}/fetch_summary.tsv"
echo ""
ok "Phylogenetic analysis complete"