#!/usr/bin/env bash
# =============================================================================
# m07_popgen.sh — Step 07: Population genetics statistics
# =============================================================================
# Computes nucleotide diversity (π), Tajima's D, haplotype diversity (Hd),
# and Fu's Fs from the whole-genome aligned FASTA, per population and combined.
# Also produces a sliding window analysis of π and Tajima's D.
# Input: aligned/all_samples_mito_aligned.fa + population codes in sample names.
# Output: popgen/ directory with TSV summaries and sliding window data.
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""
ENV_DIR=""
MIN_POP_SIZE=5
WINDOW_SIZE=200
STEP_SIZE=50
POP_REGEX='-([A-Z]+)\d'

RED='[0;31m'; GREEN='[0;32m'; CYAN='[0;36m'; YELLOW='[1;33m'; RESET='[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:p:w:s:R:" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) ENV_DIR="$OPTARG" ;;
        p) MIN_POP_SIZE="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        s) STEP_SIZE="$OPTARG" ;;
        R) POP_REGEX="$OPTARG" ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$ENV_DIR"    ]] && err "Conda environment required (-e)"

set +u
source ~/.bash_profile
conda activate "$ENV_DIR"
set -u

FASTA="${OUTPUT_DIR}/aligned/all_samples_mito_aligned.fa"
POPGEN_DIR="${OUTPUT_DIR}/popgen"

[[ ! -f "$FASTA" ]] && err "Aligned FASTA not found: $FASTA"

mkdir -p "$POPGEN_DIR"

log "Running population genetics analysis..."
log "  Input: $FASTA"
log "  Min population size: ${MIN_POP_SIZE}"
log "  Sliding window: ${WINDOW_SIZE} bp, step ${STEP_SIZE} bp"

python3 - "$FASTA" "$POPGEN_DIR" "$MIN_POP_SIZE" "$WINDOW_SIZE" "$STEP_SIZE" "$POP_REGEX" <<'PYEOF'
import sys, os, re, csv, math, itertools
from collections import defaultdict

fasta_path   = sys.argv[1]
out_dir      = sys.argv[2]
min_pop_size = int(sys.argv[3])
window_size  = int(sys.argv[4])
step_size    = int(sys.argv[5])
pop_regex    = sys.argv[6]

# ── 1. Load FASTA ─────────────────────────────────────────────────────────────
def load_fasta(path):
    seqs = {}
    name = None
    buf  = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf).lower()
                name = line[1:].split()[0]
                buf  = []
            else:
                buf.append(line)
    if name:
        seqs[name] = ''.join(buf).lower()
    return seqs

print("Loading alignment...")
seqs = load_fasta(fasta_path)
print(f"  {len(seqs)} sequences loaded")

# Sample names have _mito suffix — strip it for pop parsing
def get_pop(name):
    clean = name.replace('_mito', '')
    m = re.search(pop_regex, clean)
    return m.group(1) if m else 'Unknown'

# ── 2. Group by population ────────────────────────────────────────────────────
pop_seqs = defaultdict(dict)
for name, seq in seqs.items():
    pop = get_pop(name)
    pop_seqs[pop][name] = seq

pop_counts = {p: len(s) for p, s in pop_seqs.items()}
print(f"  Population counts: {dict(sorted(pop_counts.items()))}")

included = {p: s for p, s in pop_seqs.items()
            if len(s) >= min_pop_size and p != 'Unknown'}
excluded = {p: n for p, n in pop_counts.items()
            if n < min_pop_size or p == 'Unknown'}

if excluded:
    print(f"  Skipping (n < {min_pop_size}): {excluded}")

# Add 'all' group
included['all'] = seqs

# ── 3. Core statistics functions ──────────────────────────────────────────────

def get_seq_matrix(seq_dict):
    """Return list of sequences as list of strings, gap-stripped columns removed."""
    names = list(seq_dict.keys())
    seqs  = [seq_dict[n] for n in names]
    if not seqs:
        return []
    # Remove positions where ALL seqs have gap or N
    L = len(seqs[0])
    keep = []
    for i in range(L):
        col = [s[i] for s in seqs]
        if not all(c in '-n' for c in col):
            keep.append(i)
    return [''.join(s[i] for i in keep) for s in seqs]

def count_segregating_sites(seqs):
    """Number of sites with >= 2 distinct non-gap bases."""
    if len(seqs) < 2:
        return 0, len(seqs[0]) if seqs else 0
    L = len(seqs[0])
    S = 0
    L_used = 0
    for i in range(L):
        col = [s[i] for s in seqs if s[i] not in '-n']
        if len(col) < 2:
            continue
        L_used += 1
        bases = set(col)
        if len(bases) > 1:
            S += 1
    return S, L_used

def nucleotide_diversity(seqs):
    """π = mean pairwise differences per site."""
    n = len(seqs)
    if n < 2:
        return None
    L = len(seqs[0])
    total_diff = 0
    total_sites = 0
    n_pairs = 0
    for i, j in itertools.combinations(range(n), 2):
        s1, s2 = seqs[i], seqs[j]
        diffs = 0
        sites = 0
        for k in range(L):
            b1, b2 = s1[k], s2[k]
            if b1 in '-n' or b2 in '-n':
                continue
            sites += 1
            if b1 != b2:
                diffs += 1
        if sites > 0:
            total_diff  += diffs / sites
            total_sites += 1
        n_pairs += 1
    if n_pairs == 0:
        return None
    return total_diff / n_pairs

def tajimas_d(seqs):
    """Tajima's D statistic."""
    n = len(seqs)
    if n < 4:
        return None
    S, L = count_segregating_sites(seqs)
    if S == 0:
        return None

    # Normalise π to per-site
    pi = nucleotide_diversity(seqs)
    if pi is None:
        return None

    # Watterson's theta
    a1 = sum(1.0/i for i in range(1, n))
    theta_w = S / (a1 * L) if L > 0 else 0

    if theta_w == 0:
        return None

    # Tajima's D variance
    a2 = sum(1.0/i**2 for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - 1.0/a1
    c2 = b2 - (n + 2)/(a1 * n) + a2/a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    var_d = e1 * S + e2 * S * (S - 1)
    if var_d <= 0:
        return None

    D = (pi - theta_w) / math.sqrt(var_d)
    return D

def haplotype_diversity(seqs):
    """Haplotype diversity Hd = n/(n-1) * (1 - sum(pi^2))."""
    n = len(seqs)
    if n < 2:
        return None
    # Count unique haplotypes (ignoring gaps/N)
    def clean(s):
        return ''.join(c if c not in '-n' else '?' for c in s)
    haps = defaultdict(int)
    for s in seqs:
        haps[clean(s)] += 1
    freq = [count/n for count in haps.values()]
    hd = (n / (n - 1)) * (1 - sum(f**2 for f in freq))
    return hd, len(haps)

def fus_fs(seqs):
    """Fu's Fs statistic (approximate — uses observed haplotype count and π)."""
    n = len(seqs)
    if n < 4:
        return None
    result = haplotype_diversity(seqs)
    if result is None:
        return None
    hd, k = result
    pi = nucleotide_diversity(seqs)
    if pi is None or pi == 0:
        return None

    # Theta estimate from π
    theta = pi
    # Expected number of haplotypes under neutrality
    # Using the approximation: E[k] ≈ theta * H_n where H_n = sum(1/i)
    Hn = sum(1.0/i for i in range(1, n))
    # Fu's Fs: Fs = ln(S_obs / (1 - S_obs)) where S_obs = P(k <= k_obs | theta)
    # Simplified approximation using log-odds of observing <= k haplotypes
    # Full implementation uses the exact probability — here we use the standard
    # approximation: Fs = ln(p / (1-p)), negative values = excess haplotypes
    # We compute the probability via the Ewens sampling formula approximation
    # log(prob of k or fewer haplotypes given theta)
    # For practical purposes use the Depaulis & Veuille (1998) formula:
    # Fs = ln[ P(K <= k_obs | theta_hat) / (1 - P(K <= k_obs | theta_hat)) ]
    # theta_hat from pi
    theta_hat = theta * n

    # Compute P(K <= k | theta) via stirling numbers of the first kind (recursive)
    # for small n; for large n use normal approximation
    def ewens_prob_k_or_fewer(n, k_obs, theta):
        """P(K <= k_obs) under Ewens sampling formula."""
        # Unnormalised probabilities via unsigned Stirling numbers s(n,k)
        # s(n,k) = s(n-1,k-1) + (n-1)*s(n-1,k)
        if n > 100:
            # Normal approximation: K ~ N(mu, sigma^2)
            mu    = sum(theta/(theta + i - 1) for i in range(1, n+1))
            var   = sum(theta*(i-1)/(theta + i - 1)**2 for i in range(1, n+1))
            if var <= 0:
                return 0.5
            z = (k_obs - mu) / math.sqrt(var)
            # Standard normal CDF approximation
            def norm_cdf(x):
                t = 1/(1 + 0.2316419*abs(x))
                poly = t*(0.319381530 + t*(-0.356563782 +
                       t*(1.781477937 + t*(-1.821255978 + t*1.330274429))))
                p = 1 - (1/math.sqrt(2*math.pi)) * math.exp(-x**2/2) * poly
                return p if x >= 0 else 1 - p
            return norm_cdf(z)
        # Exact via Stirling numbers
        # s[i][j] = unsigned Stirling number of first kind
        s = [[0.0]*(n+1) for _ in range(n+1)]
        s[0][0] = 1.0
        for i in range(1, n+1):
            for j in range(1, i+1):
                s[i][j] = s[i-1][j-1] + (i-1)*s[i-1][j]
        # Rising factorial theta^(n) = theta*(theta+1)*...*(theta+n-1)
        rising = 1.0
        for i in range(n):
            rising *= (theta + i)
        probs = []
        for k in range(1, n+1):
            probs.append(s[n][k] * (theta**k) / rising)
        cumulative = sum(probs[:k_obs])
        return min(cumulative, 1.0 - 1e-10)

    p = ewens_prob_k_or_fewer(n, k, theta_hat)
    p = max(p, 1e-10)
    p = min(p, 1 - 1e-10)
    Fs = math.log(p / (1 - p))
    return Fs

# ── 4. Compute stats per population ──────────────────────────────────────────
print("Computing population statistics...")
summary_rows = []

for pop_label in sorted(included.keys()):
    pop_dict = included[pop_label]
    n = len(pop_dict)
    print(f"  {pop_label}: n={n}")

    seqs_clean = get_seq_matrix(pop_dict)
    if not seqs_clean:
        warn(f"  No usable sites for {pop_label}")
        continue

    L = len(seqs_clean[0])
    S, L_used = count_segregating_sites(seqs_clean)
    pi    = nucleotide_diversity(seqs_clean)
    td    = tajimas_d(seqs_clean)
    hd_r  = haplotype_diversity(seqs_clean)
    hd    = hd_r[0] if hd_r else None
    n_hap = hd_r[1] if hd_r else None
    fs    = fus_fs(seqs_clean)

    summary_rows.append({
        'population':       pop_label,
        'n':                n,
        'alignment_length': L,
        'sites_used':       L_used,
        'seg_sites':        S,
        'pi':               f"{pi:.6f}"    if pi  is not None else 'NA',
        'tajimas_d':        f"{td:.4f}"    if td  is not None else 'NA',
        'haplotypes':       n_hap          if n_hap is not None else 'NA',
        'hap_diversity':    f"{hd:.4f}"    if hd  is not None else 'NA',
        'fus_fs':           f"{fs:.4f}"    if fs  is not None else 'NA',
    })
    print(f"    π={pi:.6f}  D={td:.4f}  Hd={hd:.4f}  k={n_hap}  Fs={fs:.4f}"
          if all(x is not None for x in [pi, td, hd, n_hap, fs])
          else f"    Some stats unavailable (n too small?)")

# Write summary
summary_path = os.path.join(out_dir, 'popgen_summary.tsv')
with open(summary_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=[
        'population','n','alignment_length','sites_used','seg_sites',
        'pi','tajimas_d','haplotypes','hap_diversity','fus_fs'], delimiter='	')
    w.writeheader()
    w.writerows(summary_rows)
print(f"Summary written → {summary_path}")

# ── 5. Sliding window (all samples + per population) ─────────────────────────
print(f"Running sliding window (w={window_size}, step={step_size})...")

def sliding_window_stats(seq_dict, window_size, step_size):
    seqs_raw = list(seq_dict.values())
    if len(seqs_raw) < 4:
        return []
    L = len(seqs_raw[0])
    rows = []
    for start in range(0, L - window_size + 1, step_size):
        end = start + window_size
        window_seqs = [s[start:end] for s in seqs_raw]
        # Remove all-gap columns within window
        keep = [i for i in range(window_size)
                if not all(s[i] in '-n' for s in window_seqs)]
        if len(keep) < 10:
            continue
        w_seqs = [''.join(s[i] for i in keep) for s in window_seqs]
        pi = nucleotide_diversity(w_seqs)
        td = tajimas_d(w_seqs)
        rows.append({
            'start':     start + 1,
            'end':       end,
            'mid':       (start + end) // 2,
            'sites':     len(keep),
            'pi':        f"{pi:.6f}" if pi is not None else 'NA',
            'tajimas_d': f"{td:.4f}" if td is not None else 'NA',
        })
    return rows

window_rows = {}
for pop_label in sorted(included.keys()):
    print(f"  Sliding window: {pop_label}...")
    rows = sliding_window_stats(included[pop_label], window_size, step_size)
    window_rows[pop_label] = rows
    out_path = os.path.join(out_dir, f'window_{pop_label}.tsv')
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['start','end','mid','sites','pi','tajimas_d'],
                           delimiter='	')
        w.writeheader()
        w.writerows(rows)
    print(f"    {len(rows)} windows → {out_path}")

# ── 6. Write pop size JSON for report ────────────────────────────────────────
import json
pop_json = {
    'all_pops':      pop_counts,
    'included_pops': {p: len(s) for p, s in included.items() if p != 'all'},
    'excluded_pops': excluded,
    'min_pop_size':  min_pop_size,
}
with open(os.path.join(out_dir, 'popgen_pop_summary.json'), 'w') as f:
    json.dump(pop_json, f, indent=2)

print("Popgen analysis complete.")
PYEOF

ok "Population genetics complete → ${POPGEN_DIR}"
