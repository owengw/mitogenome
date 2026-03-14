#!/usr/bin/env bash
# =============================================================================
# m12_multigene.sh — Multigene phylogeny of Plestiodon longirostris
#
# Steps:
#   1. Install MrBayes into mito_beast conda env (if not present)
#   2. Fetch 12S, ND4, CYTB from GenBank via fetch_multigene.py
#   3. Add one P. longirostris representative per population to each gene
#   4. Align each gene with MAFFT
#   5. Run modeltest-ng on each gene alignment
#   6. Concatenate alignments + write partitioned MrBayes nexus
#   7. Run MrBayes (2 runs × 4 chains, 2M generations)
#   8. Run IQ-TREE ML tree on concatenated alignment
#   9. Run IQ-TREE NJ tree on CYTB alignment
#  10. Generate HTML report (figures collected by m11_report)
# =============================================================================
#SBATCH --job-name=m12_multigene
#SBATCH --output=/mnt/parscratch/users/bi4og/genome/mito_out/multigene/m12_%j.log
#SBATCH --error=/mnt/parscratch/users/bi4og/genome/mito_out/multigene/m12_%j.err
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

set -euo pipefail

# =============================================================================
# USER CONFIGURATION
# =============================================================================
EMAIL="bi4og@sheffield.ac.uk"
PIPELINE_ENV="/users/bi4og/conda_envs/mito_pipeline"
BEAST_ENV="/users/bi4og/conda_envs/mito_beast"
OUTDIR="/mnt/parscratch/users/bi4og/genome/mito_out/multigene"
SCRIPT_DIR="/mnt/parscratch/users/bi4og/scripts"
MITO_OUT="/mnt/parscratch/users/bi4og/genome/mito_out"

# Population map: sample prefix -> population code
# Adjust prefixes to match your actual sample naming
POP_MAP=(
    "SI:SI"   # Santo Island — match samples with SI in name
    "NS:NS"   # North Shore
    "CAI:CAI" # Castle Harbour Island
    "SB:SB"   # Somerset Bridge
)

MAX_PER_SPECIES=3
THREADS=${SLURM_CPUS_PER_TASK:-8}
BOOTSTRAP=1000
MRBAYES_NGEN=2000000     # MCMC generations
MRBAYES_SAMPLEFREQ=500   # sample every N generations
MRBAYES_BURNIN=25        # % burnin to discard
SKIP_FETCH=false
SKIP_MRBAYES=false
# =============================================================================

FETCH_SCRIPT="${SCRIPT_DIR}/fetch_multigene.py"
REPORT_SCRIPT="${SCRIPT_DIR}/make_multigene_report.py"
GENES=(12S ND4 CYTB)

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

# ── Parse args ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-fetch)   SKIP_FETCH=true;   shift ;;
        --skip-mrbayes) SKIP_MRBAYES=true; shift ;;
        --email)        EMAIL="$2";        shift 2 ;;
        --outdir)       OUTDIR="$2";       shift 2 ;;
        *) shift ;;
    esac
done

mkdir -p "$OUTDIR"

# ── Activate envs (PATH prepend for SLURM compatibility) ──────────────────────
export PATH="${PIPELINE_ENV}/bin:${PATH}"
export CONDA_PREFIX="${PIPELINE_ENV}"

PYTHON="${PIPELINE_ENV}/bin/python3"
MAFFT="${PIPELINE_ENV}/bin/mafft"
IQTREE="${PIPELINE_ENV}/bin/iqtree"
MODELTEST="${PIPELINE_ENV}/bin/modeltest-ng"

[[ ! -x "$MAFFT"  ]] && err "mafft not found: ${MAFFT}"
[[ ! -x "$IQTREE" ]] && err "iqtree not found: ${IQTREE}"
"$PYTHON" -c "from Bio import Entrez" 2>/dev/null || err "Biopython not found"

# ── Step 0: Install MrBayes if needed ─────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 0: Check / install MrBayes${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

MRBAYES="${BEAST_ENV}/bin/mb"
if [[ ! -x "$MRBAYES" ]]; then
    log "MrBayes not found — installing into ${BEAST_ENV}..."
    "${BEAST_ENV}/bin/conda" install -y -c bioconda mrbayes --prefix "$BEAST_ENV" \
        || "${BEAST_ENV}/bin/mamba" install -y -c bioconda mrbayes --prefix "$BEAST_ENV" \
        || err "Failed to install MrBayes. Try manually: conda install -c bioconda mrbayes"
    ok "MrBayes installed → ${MRBAYES}"
else
    ok "MrBayes found: ${MRBAYES}"
fi

# Also check/install modeltest-ng
if [[ ! -x "$MODELTEST" ]]; then
    log "modeltest-ng not found — installing into ${PIPELINE_ENV}..."
    "${PIPELINE_ENV}/bin/conda" install -y -c bioconda modeltest-ng --prefix "$PIPELINE_ENV" \
        || warn "modeltest-ng install failed — will use GTR+G for all partitions"
    MODELTEST="${PIPELINE_ENV}/bin/modeltest-ng"
fi

# ── Step 1: Fetch sequences ───────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 1: Fetch 12S, ND4, CYTB from GenBank${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

if [[ "$SKIP_FETCH" == true ]]; then
    log "Skipping fetch (--skip-fetch)"
    for gene in "${GENES[@]}"; do
        [[ ! -f "${OUTDIR}/all_${gene}_genbank.fa" ]] && \
            err "Expected ${OUTDIR}/all_${gene}_genbank.fa not found"
    done
else
    "$PYTHON" "$FETCH_SCRIPT" \
        --email "$EMAIL" \
        --outdir "$OUTDIR" \
        --max-per-species "$MAX_PER_SPECIES"

    for gene in "${GENES[@]}"; do
        n=$(grep -c "^>" "${OUTDIR}/all_${gene}_genbank.fa" 2>/dev/null || echo 0)
        ok "${gene}: ${n} GenBank sequences"
    done
fi

# ── Step 2: Add P. longirostris representatives ───────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 2: Add P. longirostris (one rep per population)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

# Source of per-gene longirostris sequences from m09 dnds alignments
# m09 produces aligned FASTAs per gene; we use those directly.
DNDS_DIR="${MITO_OUT}/dnds/alignments"

for gene in "${GENES[@]}"; do
    COMBINED="${OUTDIR}/all_${gene}_combined.fa"
    cp "${OUTDIR}/all_${gene}_genbank.fa" "$COMBINED"

    LONGI_FA="${DNDS_DIR}/${gene}_aligned.fa"
    if [[ ! -f "$LONGI_FA" ]]; then
        warn "${gene}: no longirostris alignment found at ${LONGI_FA} — skipping"
        continue
    fi

    # Python: pick first sequence per population as representative
    "$PYTHON" - "$LONGI_FA" "$COMBINED" "${POP_MAP[@]}" << 'PYEOF'
import sys
from Bio import SeqIO

longi_fa  = sys.argv[1]
combined  = sys.argv[2]
pop_map_args = sys.argv[3:]  # "SI:SI NS:NS ..." etc.

# Build prefix->pop dict
pop_map = {}
for entry in pop_map_args:
    prefix, pop = entry.split(":")
    pop_map[prefix] = pop

records = list(SeqIO.parse(longi_fa, "fasta"))
# Remove gap-only columns by ungapping each sequence
for rec in records:
    rec.seq = rec.seq.ungap("-")

# Pick first representative per population
seen_pops = {}
for rec in records:
    sid = rec.id
    for prefix, pop in pop_map.items():
        if prefix.upper() in sid.upper() and pop not in seen_pops:
            seen_pops[pop] = rec
            break

reps = list(seen_pops.values())
for pop, rec in seen_pops.items():
    rec.id          = f"Plestiodon_longirostris_{pop}"
    rec.name        = rec.id
    rec.description = ""

print(f"  Adding {len(reps)} P. longirostris representatives: {list(seen_pops.keys())}")
with open(combined, "a") as fh:
    SeqIO.write(reps, fh, "fasta")
PYEOF

    n=$(grep -c "^>" "$COMBINED")
    ok "${gene}: combined FASTA has ${n} sequences"
done

# ── Step 3: Align each gene with MAFFT ───────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 3: Align with MAFFT${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

declare -A ALN_LENS
for gene in "${GENES[@]}"; do
    COMBINED="${OUTDIR}/all_${gene}_combined.fa"
    ALIGNED="${OUTDIR}/all_${gene}_aligned.fa"

    log "Aligning ${gene}..."
    "$MAFFT" --auto --thread "$THREADS" --reorder "$COMBINED" > "$ALIGNED"

    ALN_LEN=$(grep -v "^>" "$ALIGNED" | head -1 | tr -d '\n' | wc -c)
    ALN_LENS[$gene]=$ALN_LEN
    N=$(grep -c "^>" "$ALIGNED")
    ok "${gene}: ${N} sequences, ${ALN_LEN} bp"
done

# ── Step 4: modeltest-ng per gene ─────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 4: Model selection (modeltest-ng)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

declare -A BEST_MODELS
for gene in "${GENES[@]}"; do
    ALIGNED="${OUTDIR}/all_${gene}_aligned.fa"
    MT_PREFIX="${OUTDIR}/modeltest_${gene}"

    if [[ -x "$MODELTEST" ]]; then
        log "Running modeltest-ng for ${gene}..."
        "$MODELTEST" \
            -i "$ALIGNED" \
            -o "$MT_PREFIX" \
            -t ml \
            --force \
            --threads "$THREADS" 2>/dev/null || true

        # Parse best BIC model from modeltest-ng output
        BEST=""
        if [[ -f "${MT_PREFIX}.out" ]]; then
            BEST=$(grep "BIC" "${MT_PREFIX}.out" | grep "^>" | head -1 | \
                   awk '{print $2}' 2>/dev/null || true)
        fi
        # Fallback if parsing fails
        if [[ -z "$BEST" ]]; then
            BEST="GTR+G"
            warn "${gene}: could not parse modeltest-ng output — using GTR+G"
        fi
    else
        BEST="GTR+G"
        warn "${gene}: modeltest-ng not available — using GTR+G"
    fi

    BEST_MODELS[$gene]="$BEST"
    ok "${gene}: best model = ${BEST}"
done

# ── Step 5: Concatenate alignments + write MrBayes nexus ─────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 5: Concatenate alignments + write MrBayes nexus${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

CONCAT_FA="${OUTDIR}/all_genes_concat.fa"
NEXUS="${OUTDIR}/mrbayes_input.nex"

"$PYTHON" - \
    "$OUTDIR" \
    "${GENES[@]}" \
    "${BEST_MODELS[12S]}" \
    "${BEST_MODELS[ND4]}" \
    "${BEST_MODELS[CYTB]}" \
    "$CONCAT_FA" \
    "$NEXUS" \
    "$MRBAYES_NGEN" \
    "$MRBAYES_SAMPLEFREQ" \
    "$MRBAYES_BURNIN" << 'PYEOF'
import sys
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

outdir     = sys.argv[1]
genes      = sys.argv[2:5]          # 12S ND4 CYTB
models     = sys.argv[5:8]          # model per gene
concat_fa  = sys.argv[8]
nexus_out  = sys.argv[9]
ngen       = int(sys.argv[10])
samplefreq = int(sys.argv[11])
burnin     = int(sys.argv[12])

# Load each alignment; keyed by sequence ID
aln = {}
gene_lens = {}
for gene, model in zip(genes, models):
    aligned_fa = f"{outdir}/all_{gene}_aligned.fa"
    records = {r.id: str(r.seq) for r in SeqIO.parse(aligned_fa, "fasta")}
    aln[gene]       = records
    gene_lens[gene] = len(next(iter(records.values())))

# Find taxa present in ALL three genes
taxa_per_gene = [set(aln[g].keys()) for g in genes]
shared_taxa   = taxa_per_gene[0]
for t in taxa_per_gene[1:]:
    shared_taxa = shared_taxa & t

# Also include taxa present in at least one gene, padding missing with Ns
all_taxa = set()
for g in genes:
    all_taxa |= set(aln[g].keys())

print(f"  Shared in all 3 genes: {len(shared_taxa)} taxa")
print(f"  Total unique taxa:     {len(all_taxa)} taxa")
print(f"  Gene lengths: " + ", ".join(f"{g}={gene_lens[g]}bp" for g in genes))

# Build concatenated sequences (pad with Ns where gene missing)
concat_seqs = OrderedDict()
for taxon in sorted(all_taxa):
    seq = ""
    for gene in genes:
        if taxon in aln[gene]:
            seq += aln[gene][taxon]
        else:
            seq += "N" * gene_lens[gene]
    concat_seqs[taxon] = seq

total_len = sum(gene_lens[g] for g in genes)

# Write concatenated FASTA
records = [SeqRecord(Seq.Seq(seq), id=tid, name=tid, description="")
           for tid, seq in concat_seqs.items()]
with open(concat_fa, "w") as fh:
    SeqIO.write(records, fh, "fasta")
print(f"  Concatenated FASTA: {len(records)} taxa, {total_len} bp → {concat_fa}")

# ── Write MrBayes nexus ──────────────────────────────────────────────────────
# Partition positions
pos = 1
partitions = []
for gene in genes:
    end = pos + gene_lens[gene] - 1
    partitions.append((gene, pos, end))
    pos = end + 1

def mrbayes_model_block(gene, model, part_num):
    """
    Convert modeltest-ng model string to MrBayes lset/prset commands.
    Handles common models: GTR, HKY, SYM, TrN, TVM, TIM etc.
    Falls back to GTR+G if unrecognised.
    """
    m = model.upper()
    lines = []

    # Rate variation
    if "+G4" in m or "+G" in m:
        lines.append(f"  lset applyto=({part_num}) rates=gamma ngammacat=4;")
    elif "+I+G" in m:
        lines.append(f"  lset applyto=({part_num}) rates=invgamma ngammacat=4;")
    elif "+I" in m:
        lines.append(f"  lset applyto=({part_num}) rates=propinv;")
    else:
        lines.append(f"  lset applyto=({part_num}) rates=gamma ngammacat=4;")

    # Substitution model
    if "GTR" in m:
        lines.append(f"  lset applyto=({part_num}) nst=6;")
    elif "HKY" in m or "TrN" in m:
        lines.append(f"  lset applyto=({part_num}) nst=2;")
    elif "SYM" in m or "TVM" in m or "TIM" in m or "TVMef" in m:
        lines.append(f"  lset applyto=({part_num}) nst=6;")
    else:
        lines.append(f"  lset applyto=({part_num}) nst=6;")

    return "\n".join(lines)

nsamples = ngen // samplefreq
burnin_n = int(nsamples * burnin / 100)

with open(nexus_out, "w") as fh:
    # Data block
    fh.write("#NEXUS\n\n")
    fh.write("Begin data;\n")
    fh.write(f"  Dimensions ntax={len(concat_seqs)} nchar={total_len};\n")
    fh.write("  Format datatype=dna missing=N gap=-;\n")
    fh.write("  Matrix\n")
    for tid, seq in concat_seqs.items():
        fh.write(f"    {tid:<50} {seq}\n")
    fh.write("  ;\nEnd;\n\n")

    # Sets block — define partitions
    fh.write("Begin sets;\n")
    for gene, start, end in partitions:
        fh.write(f"  charset {gene} = {start}-{end};\n")
    fh.write("End;\n\n")

    # MrBayes block
    fh.write("Begin mrbayes;\n")
    fh.write(f"  partition genes = {len(genes)}: " +
             ", ".join(g for g, _, _ in partitions) + ";\n")
    fh.write("  set partition=genes;\n\n")

    for i, (gene, model) in enumerate(zip(genes, models), 1):
        fh.write(f"  [ Gene {i}: {gene} — best model: {model} ]\n")
        fh.write(mrbayes_model_block(gene, model, i) + "\n\n")

    fh.write("  [ Unlink substitution parameters across partitions ]\n")
    fh.write("  unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);\n\n")

    fh.write("  [ MCMC settings ]\n")
    fh.write(f"  mcmcp ngen={ngen} samplefreq={samplefreq} "
             f"nruns=2 nchains=4 temp=0.2\n")
    fh.write(f"         filename={outdir}/mrbayes_output\n")
    fh.write(f"         printfreq=10000 diagnfreq=10000;\n")
    fh.write(f"  mcmc;\n\n")

    fh.write(f"  [ Summarise — discard {burnin}% burnin ]\n")
    fh.write(f"  sumt filename={outdir}/mrbayes_output burnin={burnin_n} "
             f"contype=allcompat;\n")
    fh.write(f"  sump filename={outdir}/mrbayes_output burnin={burnin_n};\n")
    fh.write("End;\n")

print(f"  MrBayes nexus → {nexus_out}")
print(f"  Partitions: " + " | ".join(f"{g} {s}-{e}" for g,s,e in partitions))
PYEOF

ok "Concatenated alignment and MrBayes nexus written"

# ── Step 6: Run MrBayes ───────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 6: MrBayes partitioned Bayesian analysis${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

if [[ "$SKIP_MRBAYES" == true ]]; then
    warn "Skipping MrBayes (--skip-mrbayes)"
else
    log "Running MrBayes (${MRBAYES_NGEN} generations, 2 runs × 4 chains)..."
    log "This typically takes 2–6 hours. Monitor: tail -f ${OUTDIR}/mrbayes_output.run1.p"

    # MrBayes must be run from the output directory
    cd "$OUTDIR"
    "$MRBAYES" mrbayes_input.nex
    cd -

    # Check convergence
    log "Checking convergence (average std deviation of split frequencies)..."
    if [[ -f "${OUTDIR}/mrbayes_output.run1.p" ]]; then
        ASDSF=$(grep "Average standard deviation" "${OUTDIR}/mrbayes_output.run1.t" \
                2>/dev/null | tail -1 | awk '{print $NF}' || echo "unknown")
        if [[ "$ASDSF" != "unknown" ]]; then
            ok "ASDSF = ${ASDSF}  (target <0.01)"
            if (( $(echo "$ASDSF > 0.05" | bc -l 2>/dev/null || echo 0) )); then
                warn "ASDSF > 0.05 — consider increasing ngen or checking mixing"
            fi
        fi
    fi

    CONTREE="${OUTDIR}/mrbayes_output.con.tre"
    [[ -f "$CONTREE" ]] && ok "Consensus tree → ${CONTREE}" || \
        warn "Consensus tree not found — check MrBayes log"
fi

# ── Step 7: IQ-TREE ML tree (concatenated) ────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 7: IQ-TREE ML tree (concatenated alignment)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

log "Running IQ-TREE on concatenated alignment..."
"$IQTREE" \
    -s "${OUTDIR}/all_genes_concat.fa" \
    -m GTR+G \
    -B "$BOOTSTRAP" \
    -T "$THREADS" \
    --prefix "${OUTDIR}/iqtree_concat" \
    --redo

ok "IQ-TREE ML tree → ${OUTDIR}/iqtree_concat.treefile"

# ── Step 8: IQ-TREE NJ tree (CYTB only) ──────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 8: NJ tree (CYTB alignment)${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

log "Running IQ-TREE NJ tree on CYTB..."
"$IQTREE" \
    -s "${OUTDIR}/all_CYTB_aligned.fa" \
    -m GTR+G \
    --tree-fix \
    -T "$THREADS" \
    --prefix "${OUTDIR}/nj_cytb" \
    --redo

ok "NJ tree → ${OUTDIR}/nj_cytb.treefile"

# ── Step 9: Generate HTML report ─────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  STEP 9: Generate HTML tree report${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"

[[ ! -f "$REPORT_SCRIPT" ]] && \
    err "make_multigene_report.py not found: ${REPORT_SCRIPT}"

"$PYTHON" "$REPORT_SCRIPT" \
    --outdir "$OUTDIR" \
    --mrbayes-tree  "${OUTDIR}/mrbayes_output.con.tre" \
    --ml-tree       "${OUTDIR}/iqtree_concat.treefile" \
    --nj-tree       "${OUTDIR}/nj_cytb.treefile" \
    --output        "${OUTDIR}/multigene_report.html"

ok "HTML report → ${OUTDIR}/multigene_report.html"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
echo -e "${BOLD}  OUTPUT SUMMARY${RESET}"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
for gene in "${GENES[@]}"; do
    echo -e "  ${gene} aligned:        ${OUTDIR}/all_${gene}_aligned.fa"
done
echo -e "  Concat alignment:    ${OUTDIR}/all_genes_concat.fa"
echo -e "  MrBayes nexus:       ${OUTDIR}/mrbayes_input.nex"
echo -e "  MrBayes tree:        ${OUTDIR}/mrbayes_output.con.tre"
echo -e "  IQ-TREE ML tree:     ${OUTDIR}/iqtree_concat.treefile"
echo -e "  NJ tree (CYTB):      ${OUTDIR}/nj_cytb.treefile"
echo -e "  HTML report:         ${OUTDIR}/multigene_report.html"
echo -e "  Fetch summary:       ${OUTDIR}/fetch_summary.tsv"
echo ""
ok "m12_multigene complete"