#!/usr/bin/env bash
# =============================================================================
# m09_dnds.sh — Step 09: dN/dS analysis with PAML (yn00 + codeml)
# =============================================================================
# Runs yn00 (pairwise dN/dS) and codeml branch model on each gene alignment.
# Uses IQ-TREE treefile from m06 with population labels as branch annotations.
# Populations with fewer than MIN_POP_SIZE samples are skipped for branch model.
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""
BEAST_ENV=""
GENETIC_CODE=2
SPECIES="Unknown species"
MIN_POP_SIZE=5

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[1;33m'; BOLD='\033[1m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:C:n:p:t:r:i:s:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) BEAST_ENV="$OPTARG" ;;
        C) GENETIC_CODE="$OPTARG" ;;
        n) SPECIES="$OPTARG" ;;
        p) MIN_POP_SIZE="$OPTARG" ;;
        t|r|i|s|q|d|f|m|k) : ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$BEAST_ENV"  ]] && err "BEAST2 conda environment required (-e)"

set +u
source ~/.bash_profile
conda activate "$BEAST_ENV"
set -u

ALN_DIR="${OUTPUT_DIR}/dnds/alignments"
DNDS_DIR="${OUTPUT_DIR}/dnds"
YN00_DIR="${DNDS_DIR}/yn00"
CODEML_DIR="${DNDS_DIR}/codeml"
TREEFILE="${OUTPUT_DIR}/tree/mito_tree.treefile"

[[ ! -d "$ALN_DIR"  ]] && err "Alignment directory not found: $ALN_DIR"
[[ ! -f "$TREEFILE" ]] && err "IQ-TREE treefile not found: $TREEFILE"

mkdir -p "$YN00_DIR" "$CODEML_DIR"

declare -A CODE_MAP=([1]=0 [2]=1 [3]=2 [4]=3 [5]=4 [6]=5 [9]=8 [13]=12 [14]=13)
PAML_ICODE="${CODE_MAP[$GENETIC_CODE]:-1}"
log "Species: ${SPECIES}"
log "Genetic code: NCBI=${GENETIC_CODE} → PAML icode=${PAML_ICODE}"
log "Minimum population size for branch model: ${MIN_POP_SIZE}"

# ── Step 1: Prepare population-labelled trees + size check ───────────────────
log "Preparing population-labelled trees for codeml branch model..."

python3 - "$TREEFILE" "$CODEML_DIR" "$MIN_POP_SIZE" <<'PYEOF'
import sys, os, re, json

treefile     = sys.argv[1]
codeml_dir   = sys.argv[2]
min_pop_size = int(sys.argv[3])

def get_pop(name):
    m = re.search(r'-([A-Z]+)\d', name)
    return m.group(1) if m else "Unknown"

with open(treefile) as f:
    tree_str = f.read().strip()

pop_counts = {}
for name in re.findall(r'[A-Za-z0-9_\-]+(?=[:,()\s])', tree_str):
    if '-' in name:
        pop = get_pop(name)
        if pop != "Unknown":
            pop_counts[pop] = pop_counts.get(pop, 0) + 1

print(f"Population sample counts: {pop_counts}")

included = {p: n for p, n in pop_counts.items() if n >= min_pop_size}
excluded = {p: n for p, n in pop_counts.items() if n < min_pop_size}

if excluded:
    print(f"WARNING: Skipping branch model for populations with n < {min_pop_size}: {excluded}")
print(f"Included for branch model: {included}")

with open(os.path.join(codeml_dir, "pop_sizes.json"), 'w') as f:
    json.dump({'min_pop_size': min_pop_size, 'all_pops': pop_counts,
               'included_pops': included, 'excluded_pops': excluded}, f, indent=2)

for focal_pop in sorted(included.keys()):
    labelled = re.sub(
        r'[A-Za-z0-9_\-]+-[A-Z]+\d+[A-Za-z0-9_\-]*',
        lambda m: m.group(0) + ' #1' if get_pop(m.group(0)) == focal_pop else m.group(0),
        tree_str
    )
    out_tree = os.path.join(codeml_dir, f"tree_{focal_pop}.nwk")
    with open(out_tree, 'w') as f:
        f.write(labelled + '\n')
    print(f"  Written: tree_{focal_pop}.nwk  (n={included[focal_pop]})")

with open(os.path.join(codeml_dir, "tree_null.nwk"), 'w') as f:
    f.write(tree_str + '\n')
print("  Written: tree_null.nwk")
PYEOF

ok "Population trees written → ${CODEML_DIR}"

# ── Step 2: Convert alignments to PHYLIP ─────────────────────────────────────
log "Converting alignments to PHYLIP format..."

python3 - "$ALN_DIR" <<'PYEOF'
import sys, os
from Bio import AlignIO

aln_dir = sys.argv[1]
for fa in sorted(os.listdir(aln_dir)):
    if not fa.endswith('_aligned.fa'):
        continue
    gene     = fa.replace('_aligned.fa', '')
    in_path  = os.path.join(aln_dir, fa)
    out_path = os.path.join(aln_dir, f'{gene}.phy')
    aln = AlignIO.read(in_path, 'fasta')
    for rec in aln:
        rec.id = rec.id[:30]; rec.description = ''
    AlignIO.write(aln, out_path, 'phylip-relaxed')
    print(f"  {gene}: {len(aln)} sequences, {aln.get_alignment_length()} bp → {gene}.phy")
print("PHYLIP conversion complete.")
PYEOF

ok "PHYLIP files written → ${ALN_DIR}"

# ── Step 3: Run yn00 (all samples, no population filter needed) ───────────────
log "Running yn00 (pairwise dN/dS) for each gene..."

for PHY in "${ALN_DIR}"/*.phy; do
    GENE=$(basename "$PHY" .phy)
    WORKDIR="${YN00_DIR}/${GENE}"
    mkdir -p "$WORKDIR"
    cp "$PHY" "${WORKDIR}/seqfile.phy"

    cat > "${WORKDIR}/yn00.ctl" <<CTLEOF
seqfile   = seqfile.phy
outfile   = yn00_results.txt
verbose   = 0
icode     = ${PAML_ICODE}
weighting = 0
CTLEOF

    log "  Running yn00 for ${GENE}..."
    (cd "$WORKDIR" && yn00 yn00.ctl > yn00.log 2>&1) \
        && ok "  ${GENE} done" \
        || warn "  yn00 failed for ${GENE} — check ${WORKDIR}/yn00.log"
done

ok "yn00 complete → ${YN00_DIR}"

# ── Step 4: Run codeml (included populations only) ────────────────────────────
log "Running codeml branch model (populations with n >= ${MIN_POP_SIZE})..."

for PHY in "${ALN_DIR}"/*.phy; do
    GENE=$(basename "$PHY" .phy)
    log "  codeml for gene: ${GENE}"

    NULL_DIR="${CODEML_DIR}/${GENE}/null"
    mkdir -p "$NULL_DIR"
    cp "$PHY" "${NULL_DIR}/seqfile.phy"
    cp "${CODEML_DIR}/tree_null.nwk" "${NULL_DIR}/tree.nwk"

    cat > "${NULL_DIR}/codeml.ctl" <<CTLEOF
seqfile   = seqfile.phy
treefile  = tree.nwk
outfile   = codeml_null.txt
noisy     = 0
verbose   = 0
runmode   = 0
seqtype   = 1
codonFreq = 2
model     = 0
NSsites   = 0
icode     = ${PAML_ICODE}
fix_kappa = 0
kappa     = 2
fix_omega = 0
omega     = 0.5
CTLEOF

    (cd "$NULL_DIR" && codeml codeml.ctl > codeml.log 2>&1) \
        || warn "  codeml null failed for ${GENE}"

    # Only tree files for included populations will exist — loop over what's there
    for POPTREE in "${CODEML_DIR}"/tree_*.nwk; do
        POP=$(basename "$POPTREE" .nwk | sed 's/tree_//')
        [[ "$POP" == "null" ]] && continue

        POP_DIR="${CODEML_DIR}/${GENE}/${POP}"
        mkdir -p "$POP_DIR"
        cp "$PHY" "${POP_DIR}/seqfile.phy"
        cp "$POPTREE" "${POP_DIR}/tree.nwk"

        cat > "${POP_DIR}/codeml.ctl" <<CTLEOF
seqfile   = seqfile.phy
treefile  = tree.nwk
outfile   = codeml_branch.txt
noisy     = 0
verbose   = 0
runmode   = 0
seqtype   = 1
codonFreq = 2
model     = 2
NSsites   = 0
icode     = ${PAML_ICODE}
fix_kappa = 0
kappa     = 2
fix_omega = 0
omega     = 0.5
CTLEOF

        (cd "$POP_DIR" && codeml codeml.ctl > codeml.log 2>&1) \
            || warn "  codeml branch failed for ${GENE}/${POP}"
    done

    ok "  ${GENE} codeml done"
done

ok "codeml complete → ${CODEML_DIR}"

# ── Step 5: Parse results ─────────────────────────────────────────────────────
log "Parsing results into summary tables..."

python3 - "$YN00_DIR" "$CODEML_DIR" "$DNDS_DIR" <<'PYEOF'
import sys, os, re, csv, json

yn00_dir   = sys.argv[1]
codeml_dir = sys.argv[2]
out_dir    = sys.argv[3]

pop_size_file = os.path.join(codeml_dir, "pop_sizes.json")
pop_sizes, excluded_pops = {}, {}
if os.path.exists(pop_size_file):
    with open(pop_size_file) as f:
        d = json.load(f)
        pop_sizes     = d.get('all_pops', {})
        excluded_pops = d.get('excluded_pops', {})

# yn00
yn00_rows = []
for gene in sorted(os.listdir(yn00_dir)):
    result_file = os.path.join(yn00_dir, gene, 'yn00_results.txt')
    if not os.path.exists(result_file): continue
    with open(result_file) as f: content = f.read()
    dn_vals, ds_vals, omega_vals = [], [], []
    for line in content.split('\n'):
        parts = line.split()
        if len(parts) >= 7:
            try:
                ds = float(parts[4]); dn = float(parts[5]); omega = float(parts[6])
                if 0 <= omega <= 99 and ds >= 0:
                    dn_vals.append(dn); ds_vals.append(ds); omega_vals.append(omega)
            except (ValueError, IndexError): pass
    if omega_vals:
        yn00_rows.append({'gene': gene, 'n_pairs': len(omega_vals),
            'mean_dN': f"{sum(dn_vals)/len(dn_vals):.6f}",
            'mean_dS': f"{sum(ds_vals)/len(ds_vals):.6f}",
            'mean_dNdS': f"{sum(omega_vals)/len(omega_vals):.6f}",
            'min_dNdS': f"{min(omega_vals):.6f}",
            'max_dNdS': f"{max(omega_vals):.6f}"})

with open(os.path.join(out_dir, 'yn00_summary.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['gene','n_pairs','mean_dN','mean_dS',
                                       'mean_dNdS','min_dNdS','max_dNdS'], delimiter='\t')
    w.writeheader(); w.writerows(yn00_rows)
print(f"yn00 summary written.")

# codeml
def parse_lnl(path):
    if not os.path.exists(path): return None
    with open(path) as f:
        for line in f:
            m = re.search(r'lnL\s*\([^)]+\)\s*=\s*([-\d.]+)', line)
            if m: return float(m.group(1))
    return None

def parse_omega(path):
    if not os.path.exists(path): return None, None
    with open(path) as f:
        for line in f:
            m = re.search(r'w \(dN/dS\) for branches:\s*([\d.]+)\s+([\d.]+)', line)
            if m: return float(m.group(1)), float(m.group(2))
    return None, None

codeml_rows = []
for gene in sorted(os.listdir(codeml_dir)):
    gene_dir = os.path.join(codeml_dir, gene)
    if not os.path.isdir(gene_dir): continue
    null_lnl = parse_lnl(os.path.join(gene_dir, 'null', 'codeml_null.txt'))
    if null_lnl is None: continue
    for pop_dir in sorted(os.listdir(gene_dir)):
        if pop_dir == 'null': continue
        pop_path = os.path.join(gene_dir, pop_dir)
        if not os.path.isdir(pop_path): continue
        alt_lnl    = parse_lnl(os.path.join(pop_path, 'codeml_branch.txt'))
        bg_w, fg_w = parse_omega(os.path.join(pop_path, 'codeml_branch.txt'))
        if alt_lnl is None: continue
        lrt = 2 * (alt_lnl - null_lnl)
        codeml_rows.append({
            'gene': gene, 'population': pop_dir,
            'n_samples': pop_sizes.get(pop_dir, 'NA'),
            'lnL_null': f"{null_lnl:.4f}", 'lnL_branch': f"{alt_lnl:.4f}",
            'LRT': f"{lrt:.4f}", 'p<0.05': "YES" if lrt > 3.841 else "no",
            'omega_background': f"{bg_w:.4f}" if bg_w is not None else 'NA',
            'omega_foreground': f"{fg_w:.4f}" if fg_w is not None else 'NA'})

if excluded_pops:
    with open(os.path.join(out_dir, 'codeml_excluded_pops.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['population','n_samples','reason'], delimiter='\t')
        w.writeheader()
        for pop, n in sorted(excluded_pops.items()):
            w.writerow({'population': pop, 'n_samples': n,
                        'reason': 'below min_pop_size threshold'})
    print(f"Excluded populations written.")

with open(os.path.join(out_dir, 'codeml_branch_summary.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['gene','population','n_samples','lnL_null',
                                       'lnL_branch','LRT','p<0.05',
                                       'omega_background','omega_foreground'], delimiter='\t')
    w.writeheader(); w.writerows(codeml_rows)
print(f"codeml summary written.")
PYEOF

ok "Results parsed → ${DNDS_DIR}/yn00_summary.tsv and ${DNDS_DIR}/codeml_branch_summary.tsv"