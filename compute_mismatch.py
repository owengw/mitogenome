#!/usr/bin/env python3
# =============================================================================
# compute_mismatch.py — Mismatch distribution analysis
#
# For each population:
#   1. Compute observed pairwise difference distribution
#   2. Fit Rogers & Harpending (1992) sudden expansion model
#   3. Calculate raggedness index (r) and SSD statistic
#   4. Output JSON summary + per-population TSV data
#
# The JSON is read by mito_report.py to render mismatch plots in the HTML report.
#
# Rogers & Harpending model: expected mismatch under sudden population expansion
#   f(x) = (theta_1 / (theta_0 + theta_1)) * (theta_0 / (theta_0 + theta_1))^x
#   Parameters: theta_0 (ancestral), theta_1 (modern), tau (expansion time)
# =============================================================================

import argparse
import json
import os
import sys
from collections import defaultdict, Counter
from itertools import combinations
import math

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--aligned",  required=True)
parser.add_argument("--outdir",   required=True)
parser.add_argument("--pop-defs", nargs="+", default=[])
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# Parse population definitions
pop_patterns = []
pop_colours  = {}
for pd in args.pop_defs:
    parts = pd.split(":")
    if len(parts) >= 2:
        pattern = parts[0]
        pop     = parts[1]
        colour  = parts[2] if len(parts) > 2 else "#94a3b8"
        pop_patterns.append((pattern, pop))
        pop_colours[pop] = colour

def get_pop(sample_id):
    for pattern, pop in pop_patterns:
        if pattern.upper() in sample_id.upper():
            return pop
    return "UNK"

# Load sequences
records = list(SeqIO.parse(args.aligned, "fasta"))
print(f"Loaded {len(records)} sequences")

# Group by population
pop_seqs = defaultdict(list)
for rec in records:
    pop = get_pop(rec.id)
    pop_seqs[pop].append(str(rec.seq).upper())

# Also compute for all combined
all_seqs = [str(r.seq).upper() for r in records]
pop_seqs["all"] = all_seqs
pop_colours["all"] = "#e2e8f0"

def pairwise_differences(seqs):
    """Count pairwise differences for all sequence pairs."""
    diffs = []
    for s1, s2 in combinations(seqs, 2):
        d = sum(1 for a, b in zip(s1, s2)
                if a != b and a not in "-N" and b not in "-N")
        diffs.append(d)
    return diffs

def raggedness(observed_freq, expected_freq):
    """
    Harpending's raggedness index r.
    r = sum((f_i - f_{i+1})^2) over all i
    Small r = smooth (consistent with expansion)
    Large r = ragged (consistent with stable/declining population)
    """
    obs = list(observed_freq)
    if len(obs) < 2:
        return 0.0
    r = sum((obs[i] - obs[i+1])**2 for i in range(len(obs)-1))
    return r

def ssd_statistic(observed_freq, expected_freq):
    """
    Sum of squared deviations between observed and expected.
    """
    n = min(len(observed_freq), len(expected_freq))
    ssd = sum((observed_freq[i] - expected_freq[i])**2 for i in range(n))
    return ssd

def rogers_harpending_expected(max_diff, theta0, theta1, tau, n_pairs):
    """
    Rogers & Harpending (1992) sudden expansion model.
    Returns expected frequency distribution.

    Under sudden expansion:
      mismatch_distribution(x) approximated by Poisson with mean tau
      for large theta_1 (modern pop size >> ancestral):
      f(x) ~ (tau^x * exp(-tau)) / x!

    Full model accounts for theta_0 and theta_1, but for most purposes
    the Poisson approximation (valid when theta_1 >> theta_0) is used.
    """
    expected = []
    for x in range(max_diff + 1):
        # Poisson approximation (Rogers & Harpending eq. 7 simplified)
        try:
            prob = (math.exp(-tau) * (tau ** x)) / math.factorial(x)
        except (OverflowError, ValueError):
            prob = 0.0
        expected.append(prob * n_pairs)
    return expected

def fit_expansion_model(diffs):
    """
    Fit sudden expansion model by minimising SSD over tau.
    Simple grid search over tau values.
    Returns (tau, theta0, theta1, expected_dist, ssd, r).
    """
    if not diffs:
        return None

    n_pairs  = len(diffs)
    max_diff = max(diffs) if diffs else 0
    counts   = Counter(diffs)
    observed = [counts.get(i, 0) for i in range(max_diff + 1)]
    obs_freq = [c / n_pairs for c in observed]

    # Mean of observed distribution ≈ tau under expansion model
    mean_diff = sum(d for d in diffs) / n_pairs if n_pairs > 0 else 0

    # Grid search for best tau
    best_tau = mean_diff
    best_ssd = float("inf")
    for tau_try in [x * 0.1 for x in range(1, int(max_diff * 30) + 1)]:
        exp = rogers_harpending_expected(max_diff, 0.01, 100, tau_try, n_pairs)
        exp_freq = [e / n_pairs for e in exp]
        s = ssd_statistic(obs_freq, exp_freq)
        if s < best_ssd:
            best_ssd = s
            best_tau = tau_try

    expected = rogers_harpending_expected(max_diff, 0.01, 100, best_tau, n_pairs)
    exp_freq = [e / n_pairs for e in expected]

    r = raggedness(obs_freq, exp_freq)

    return {
        "tau":       round(best_tau, 3),
        "theta0":    0.01,
        "theta1":    100,
        "ssd":       round(best_ssd, 6),
        "raggedness": round(r, 6),
        "n_pairs":   n_pairs,
        "mean_diff": round(mean_diff, 3),
        "max_diff":  max_diff,
        "observed":  observed,
        "observed_freq": [round(f, 6) for f in obs_freq],
        "expected":  [round(e, 2) for e in expected],
        "expected_freq": [round(f, 6) for f in exp_freq],
        "x_values":  list(range(max_diff + 1)),
    }

# ── Run analysis per population ───────────────────────────────────────────────
summary = {}

for pop in sorted(pop_seqs.keys()):
    seqs = pop_seqs[pop]
    n    = len(seqs)
    print(f"\n  Population: {pop} (n={n})")

    if n < 4:
        print(f"    Skipping — too few sequences (n<4)")
        summary[pop] = {
            "n": n, "skipped": True,
            "reason": "too few sequences (n<4)",
            "colour": pop_colours.get(pop, "#94a3b8")
        }
        continue

    diffs = pairwise_differences(seqs)
    result = fit_expansion_model(diffs)

    if result is None:
        summary[pop] = {"n": n, "skipped": True,
                         "colour": pop_colours.get(pop, "#94a3b8")}
        continue

    result["n"] = n
    result["population"] = pop
    result["colour"] = pop_colours.get(pop, "#94a3b8")

    print(f"    n pairs:     {result['n_pairs']}")
    print(f"    mean diff:   {result['mean_diff']:.3f}")
    print(f"    tau (fit):   {result['tau']:.3f}")
    print(f"    SSD:         {result['ssd']:.6f}")
    print(f"    raggedness:  {result['raggedness']:.6f}")
    if result['raggedness'] < 0.02:
        print(f"    → Consistent with sudden expansion (r < 0.02)")
    else:
        print(f"    → Raggedness suggests stable/structured population (r ≥ 0.02)")

    summary[pop] = result

    # Write per-population TSV
    tsv_path = os.path.join(args.outdir, f"mismatch_{pop}.tsv")
    with open(tsv_path, "w") as fh:
        fh.write("x\tobserved\tobserved_freq\texpected\texpected_freq\n")
        for i, x in enumerate(result["x_values"]):
            fh.write(f"{x}\t{result['observed'][i]}\t"
                     f"{result['observed_freq'][i]:.6f}\t"
                     f"{result['expected'][i]:.2f}\t"
                     f"{result['expected_freq'][i]:.6f}\n")

# ── Write JSON summary ────────────────────────────────────────────────────────
json_path = os.path.join(args.outdir, "mismatch_summary.json")
with open(json_path, "w") as fh:
    json.dump(summary, fh, indent=2)

print(f"\n{'='*60}")
print(f"  Mismatch summary JSON → {json_path}")
print(f"  Populations analysed: {len([p for p in summary if not summary[p].get('skipped')])}")