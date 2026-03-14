#!/usr/bin/env bash
# =============================================================================
# m10_bsp.sh — Step 10: Bayesian Skyline Plot (BEAST2 v2.6.7)
# =============================================================================
# Runs BEAST2 BSP analysis on CYTB alignment:
#   - Per-population runs (n >= MIN_POP_SIZE) + 1 combined run
#   - Runs execute IN PARALLEL (background jobs, wait at end)
# Clock: uncorrelated lognormal relaxed clock
# Substitution: HKY+G
# Clock rate prior: uniform(0.005, 0.02) subs/site/Myr
# Coalescent: Bayesian Skyline
# Chain: 300M (doubled from 150M for ESS > 200)
# =============================================================================

set -euo pipefail

OUTPUT_DIR=""; CONDA_ENV=""
BEAST_ENV="";  BEAST_BIN=""
MIN_POP_SIZE=5
CLOCK_RATE_MIN=0.005
CLOCK_RATE_MAX=0.02
SPECIES="Unknown species"
POP_REGEX='-([A-Z]+)\d'

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[0;33m'; RESET='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${RESET}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔  $*${RESET}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠  $*${RESET}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘  $*${RESET}" >&2; exit 1; }

while getopts ":o:e:b:p:l:u:n:R:t:r:i:s:q:d:f:m:k" opt; do
    case $opt in
        o) OUTPUT_DIR="$OPTARG" ;;
        e) CONDA_ENV="$OPTARG" ;;
        b) BEAST_BIN="$OPTARG" ;;
        p) MIN_POP_SIZE="$OPTARG" ;;
        l) CLOCK_RATE_MIN="$OPTARG" ;;
        u) CLOCK_RATE_MAX="$OPTARG" ;;
        n) SPECIES="$OPTARG" ;;
        R) POP_REGEX="$OPTARG" ;;
        r|i|s|t|q|d|f|m|k) : ;;
        :) err "Option -$OPTARG requires an argument." ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

[[ -z "$OUTPUT_DIR" ]] && err "Output directory required (-o)"
[[ -z "$CONDA_ENV"  ]] && err "Conda environment required (-e)"

# BEAST_ENV: conda env root (for treeannotator etc); BEAST_BIN: beast binary
[[ -z "$BEAST_ENV"  ]] && BEAST_ENV="$CONDA_ENV"
[[ -z "$BEAST_BIN"  ]] && BEAST_BIN="${BEAST_ENV}/bin/beast"
# If BEAST_BIN was set directly, derive BEAST_ENV from it (strip /bin/beast)
[[ "$BEAST_ENV" == "$CONDA_ENV" && -n "$BEAST_BIN" ]] && \
    BEAST_ENV="$(dirname "$(dirname "$BEAST_BIN")")"
[[ ! -f "$BEAST_BIN" ]] && err "BEAST2 binary not found: ${BEAST_BIN}"

set +u
source ~/.bash_profile
conda activate "$CONDA_ENV"
export PATH="${BEAST_ENV}/bin:${PATH}"
set -u

ALN_DIR="${OUTPUT_DIR}/dnds/alignments"
BSP_DIR="${OUTPUT_DIR}/bsp"
CYTB_ALN="${ALN_DIR}/CYTB_aligned.fa"
META="${OUTPUT_DIR}/bsp_pop_summary.json"

log "Species: ${SPECIES}"
log "BEAST2 binary: ${BEAST_BIN}"
log "Clock rate prior: Uniform(${CLOCK_RATE_MIN}, ${CLOCK_RATE_MAX}) subs/site/Myr"
log "Minimum population size for BSP: ${MIN_POP_SIZE}"

[[ ! -f "$CYTB_ALN" ]] && err "CYTB alignment not found: $CYTB_ALN"

mkdir -p "$BSP_DIR"

# ── Step 1: Split CYTB alignment by population ───────────────────────────────
log "Splitting CYTB alignment by population..."

python3 - "$CYTB_ALN" "$BSP_DIR" "$MIN_POP_SIZE" "$POP_REGEX" <<'PYEOF'
import sys, os, re, json
from Bio import SeqIO

cytb_path    = sys.argv[1]
bsp_dir      = sys.argv[2]
min_pop_size = int(sys.argv[3])
pop_regex    = sys.argv[4]

records = list(SeqIO.parse(cytb_path, 'fasta'))
print(f"Total sequences: {len(records)}")

pop_records = {}
for rec in records:
    m = re.search(pop_regex, rec.id)
    if m:
        pop_records.setdefault(m.group(1), []).append(rec)

included_pops = {}
excluded_pops = {}
for pop, recs in sorted(pop_records.items()):
    n = len(recs)
    if n >= min_pop_size:
        included_pops[pop] = n
        out_fa = os.path.join(bsp_dir, f'CYTB_{pop}.fa')
        SeqIO.write(recs, out_fa, 'fasta')
        print(f"  {pop}: n={n} \u2192 {out_fa}")
    else:
        excluded_pops[pop] = n
        print(f"  {pop}: n={n} \u2014 SKIPPED (below min_pop_size={min_pop_size})")

all_out = os.path.join(bsp_dir, 'CYTB_all.fa')
SeqIO.write(records, all_out, 'fasta')
print(f"  ALL: {len(records)} sequences \u2192 {all_out}")

with open(os.path.join(bsp_dir, 'bsp_pop_summary.json'), 'w') as jf:
    json.dump({'min_pop_size': min_pop_size,
               'all_pops':      {p: len(r) for p, r in pop_records.items()},
               'included_pops': included_pops,
               'excluded_pops': excluded_pops}, jf, indent=2)

print(f"Included: {sorted(included_pops.keys())}")
if excluded_pops:
    print(f"Excluded: {excluded_pops}")
PYEOF

ok "CYTB alignments split \u2192 ${BSP_DIR}"

# ── Step 2: Generate BEAST2 XML for each dataset ──────────────────────────────
# Remove any stale XMLs so they are always regenerated fresh
rm -f "${BSP_DIR}"/BSP_*.xml
log "Generating BEAST2 XML files..."

python3 - "$BSP_DIR" "$CLOCK_RATE_MIN" "$CLOCK_RATE_MAX" "$SPECIES" <<'PYEOF'
import sys
import os
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

bsp_dir        = sys.argv[1]
clock_rate_min = float(sys.argv[2])
clock_rate_max = float(sys.argv[3])
species        = sys.argv[4] if len(sys.argv) > 4 else "Unknown species"

def make_beast2_xml(fasta_path, run_name, out_xml, clock_rate_min, clock_rate_max, species):
    """Generate BEAST2 XML for BSP with:
    - HKY+G substitution model
    - Uncorrelated lognormal relaxed clock
    - Uniform clock rate prior parameterised from args
    - Bayesian Skyline coalescent (10 groups)
    - 300M chain, sampled every 15000
    """
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    n_seq   = len(records)
    aln_len = len(records[0].seq)

    # Bayesian Skyline groups: min(10, n_seq-1)
    n_groups = min(10, max(2, n_seq - 1))

    # Format sequence data block
    seq_lines = '\n'.join(
        f'                    <sequence id="seq_{i}_{r.id}" taxon="{r.id}" '
        f'totalcount="4" value="{str(r.seq)}"/>'
        for i, r in enumerate(records)
    )

    taxa_lines = '\n'.join(
        f'            <taxon id="{r.id}" spec="Taxon"/>'
        for r in records
    )

    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate="Standard" beautistatus="" namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.6">

    <!-- BEAST2 BSP XML for {run_name} -->
    <!-- Species: {species} -->
    <!-- {n_seq} sequences, {aln_len} bp CYTB -->
    <!-- Relaxed clock, HKY+G, Bayesian Skyline coalescent -->
    <!-- Clock rate prior: Uniform(0.005, 0.02) subs/site/Myr -->

    <data id="{run_name}" dataType="nucleotide" name="alignment">
{seq_lines}
    </data>

    <map name="Uniform">beast.math.distributions.Uniform</map>
    <map name="Exponential">beast.math.distributions.Exponential</map>
    <map name="LogNormal">beast.math.distributions.LogNormalDistributionModel</map>
    <map name="Normal">beast.math.distributions.Normal</map>
    <map name="Beta">beast.math.distributions.Beta</map>
    <map name="Gamma">beast.math.distributions.Gamma</map>
    <map name="LaplaceDistribution">beast.math.distributions.LaplaceDistribution</map>
    <map name="prior">beast.math.distributions.Prior</map>
    <map name="InverseGamma">beast.math.distributions.InverseGamma</map>
    <map name="OneOnX">beast.math.distributions.OneOnX</map>

    <run id="mcmc" spec="MCMC" chainLength="300000000">
        <state id="state" storeEvery="15000">

            <!-- Tree -->
            <tree id="Tree.t:{run_name}" name="stateNode">
                <taxonset id="TaxonSet.{run_name}" spec="TaxonSet">
{taxa_lines}
                </taxonset>
            </tree>

            <!-- Bayesian Skyline population sizes -->
            <parameter id="bPopSizes.t:{run_name}" spec="parameter.RealParameter"
                       dimension="{n_groups}" lower="0.0" name="stateNode" value="1.0"/>
            <parameter id="bGroupSizes.t:{run_name}" spec="parameter.IntegerParameter"
                       dimension="{n_groups}" lower="1" name="stateNode" value="1"/>

            <!-- HKY parameters -->
            <parameter id="kappa.s:{run_name}" lower="0.0" name="stateNode" value="2.0"/>
            <parameter id="freqParameter.s:{run_name}" dimension="4" lower="0.0"
                       name="stateNode" upper="1.0" value="0.25"/>

            <!-- Gamma shape -->
            <parameter id="gammaShape.s:{run_name}" name="stateNode" value="0.5"/>

            <!-- Clock rate: uniform prior 0.005-0.02 subs/site/Myr -->
            <parameter id="ucldMean.c:{run_name}" lower="{clock_rate_min}" upper="{clock_rate_max}"
                       name="stateNode" value="0.01"/>
            <parameter id="ucldStdev.c:{run_name}" lower="0.0" name="stateNode"
                       upper="5.0" value="0.5"/>
            <parameter id="rateCategories.c:{run_name}" spec="parameter.IntegerParameter"
                       dimension="{2*n_seq - 2}" lower="1" name="stateNode" value="1"/>

        </state>

        <!-- Initialiser -->
        <init id="RandomTree.t:{run_name}"
              spec="beast.evolution.tree.RandomTree"
              estimate="false"
              initial="@Tree.t:{run_name}"
              taxa="@{run_name}">
            <populationModel id="ConstantPopulation.t:{run_name}"
                             spec="beast.evolution.tree.coalescent.ConstantPopulation">
                <parameter id="randomPopSize.t:{run_name}" name="popSize" value="1.0"/>
            </populationModel>
        </init>

        <!-- Posterior -->
        <distribution id="posterior" spec="util.CompoundDistribution">

            <distribution id="prior" spec="util.CompoundDistribution">

                <!-- Bayesian Skyline prior -->
                <distribution id="BayesianSkyline.t:{run_name}"
                              spec="beast.evolution.tree.coalescent.BayesianSkyline"
                              groupSizes="@bGroupSizes.t:{run_name}"
                              popSizes="@bPopSizes.t:{run_name}">
                    <treeIntervals id="TreeIntervals.t:{run_name}"
                                   spec="beast.evolution.tree.coalescent.TreeIntervals"
                                   tree="@Tree.t:{run_name}"/>
                </distribution>

                <!-- Clock rate prior: Uniform(0.005, 0.02) -->
                <prior id="ucldMeanPrior.c:{run_name}" name="distribution"
                       x="@ucldMean.c:{run_name}">
                    <Uniform lower="{clock_rate_min}" upper="{clock_rate_max}" name="distr"/>
                </prior>

                <!-- ucldStdev prior: Exponential(mean=0.3333) -->
                <prior id="ucldStdevPrior.c:{run_name}" name="distribution"
                       x="@ucldStdev.c:{run_name}">
                    <Exponential mean="0.3333" name="distr"/>
                </prior>

                <!-- kappa prior -->
                <prior id="KappaPrior.s:{run_name}" name="distribution"
                       x="@kappa.s:{run_name}">
                    <LogNormal M="1.0" S="1.25" name="distr"/>
                </prior>

                <!-- Gamma shape prior -->
                <prior id="GammaShapePrior.s:{run_name}" name="distribution"
                       x="@gammaShape.s:{run_name}">
                    <Exponential mean="0.5" name="distr"/>
                </prior>

                <!-- Population size prior: 1/x -->
                <prior id="popSizesPrior.t:{run_name}" name="distribution"
                       x="@bPopSizes.t:{run_name}">
                    <OneOnX name="distr"/>
                </prior>

            </distribution>

            <!-- Likelihood -->
            <distribution id="likelihood" spec="util.CompoundDistribution"
                          useThreads="true">
                <distribution id="treeLikelihood.{run_name}"
                              spec="ThreadedTreeLikelihood"
                              data="@{run_name}" tree="@Tree.t:{run_name}">

                    <!-- HKY+G site model -->
                    <siteModel id="SiteModel.s:{run_name}" spec="SiteModel"
                               gammaCategoryCount="4">
                        <parameter id="mutationRate.s:{run_name}" estimate="false"
                                   name="mutationRate" value="1.0"/>
                        <parameter id="proportionInvariant.s:{run_name}" estimate="false"
                                   lower="0.0" name="proportionInvariant"
                                   upper="1.0" value="0.0"/>
                        <substModel id="hky.s:{run_name}" spec="HKY"
                                    kappa="@kappa.s:{run_name}">
                            <frequencies id="estimatedFreqs.s:{run_name}"
                                         spec="Frequencies"
                                         frequencies="@freqParameter.s:{run_name}"/>
                        </substModel>
                        <shape idref="gammaShape.s:{run_name}"/>
                    </siteModel>

                    <!-- Relaxed clock -->
                    <branchRateModel id="RelaxedClock.c:{run_name}"
                                     spec="beast.evolution.branchratemodel.UCRelaxedClockModel"
                                     rateCategories="@rateCategories.c:{run_name}"
                                     clock.rate="@ucldMean.c:{run_name}"
                                     tree="@Tree.t:{run_name}">
                        <LogNormal id="LogNormalDistributionModel.c:{run_name}"
                                   name="distr" M="0.0"
                                   S="@ucldStdev.c:{run_name}"
                                   meanInRealSpace="false"/>
                    </branchRateModel>

                </distribution>
            </distribution>
        </distribution>

        <!-- ── Operators ──────────────────────────────────────────────── -->

        <!-- Tree operators -->
        <operator id="CoalescentExchangeNarrow.t:{run_name}"
                  spec="Exchange" tree="@Tree.t:{run_name}" weight="15.0"/>
        <operator id="CoalescentExchangeWide.t:{run_name}"
                  spec="Exchange" isNarrow="false" tree="@Tree.t:{run_name}"
                  weight="3.0"/>
        <operator id="CoalescentTreeScaler.t:{run_name}"
                  spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:{run_name}"
                  weight="3.0"/>
        <operator id="CoalescentTreeRootScaler.t:{run_name}"
                  spec="ScaleOperator" rootOnly="true" scaleFactor="0.5"
                  tree="@Tree.t:{run_name}" weight="3.0"/>
        <operator id="CoalescentUniformOperator.t:{run_name}"
                  spec="Uniform" tree="@Tree.t:{run_name}" weight="30.0"/>
        <operator id="CoalescentSubtreeSlide.t:{run_name}"
                  spec="SubtreeSlide" tree="@Tree.t:{run_name}" weight="15.0"/>
        <operator id="CoalescentWilsonBalding.t:{run_name}"
                  spec="WilsonBalding" tree="@Tree.t:{run_name}" weight="3.0"/>

        <!-- BSP operators -->
        <operator id="popSizesScaler.t:{run_name}"
                  spec="ScaleOperator" parameter="@bPopSizes.t:{run_name}"
                  scaleFactor="0.75" weight="15.0"/>
        <operator id="groupSizesDelta.t:{run_name}"
                  spec="DeltaExchangeOperator"
                  intparameter="@bGroupSizes.t:{run_name}"
                  delta="1" weight="5.0"/>

        <!-- HKY operators -->
        <operator id="KappaScaler.s:{run_name}"
                  spec="ScaleOperator" parameter="@kappa.s:{run_name}"
                  scaleFactor="0.5" weight="0.1"/>
        <operator id="FrequenciesExchanger.s:{run_name}"
                  spec="DeltaExchangeOperator" delta="0.01"
                  parameter="@freqParameter.s:{run_name}" weight="0.1"/>
        <operator id="gammaShapeScaler.s:{run_name}"
                  spec="ScaleOperator" parameter="@gammaShape.s:{run_name}"
                  scaleFactor="0.5" weight="0.1"/>

        <!-- Relaxed clock operators -->
        <operator id="ucldMeanScaler.c:{run_name}"
                  spec="ScaleOperator" parameter="@ucldMean.c:{run_name}"
                  scaleFactor="0.5" weight="1.0"/>
        <operator id="ucldStdevScaler.c:{run_name}"
                  spec="ScaleOperator" parameter="@ucldStdev.c:{run_name}"
                  scaleFactor="0.5" weight="3.0"/>
        <operator id="rateCategoriesRandomWalk.c:{run_name}"
                  spec="IntRandomWalkOperator"
                  parameter="@rateCategories.c:{run_name}"
                  windowSize="1" weight="10.0"/>
        <operator id="rateCategoriesSwapOperator.c:{run_name}"
                  spec="SwapOperator"
                  intparameter="@rateCategories.c:{run_name}" weight="10.0"/>
        <operator id="rateCategoriesUniformOperator.c:{run_name}"
                  spec="UniformOperator"
                  parameter="@rateCategories.c:{run_name}" weight="10.0"/>

        <!-- ── Loggers ─────────────────────────────────────────────────── -->

        <logger id="tracelog" fileName="{run_name}.log" logEvery="15000"
                model="@posterior" sanitiseHeaders="true" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.{run_name}"/>
            <log idref="BayesianSkyline.t:{run_name}"/>
            <log idref="ucldMean.c:{run_name}"/>
            <log idref="ucldStdev.c:{run_name}"/>
            <log idref="kappa.s:{run_name}"/>
            <log idref="gammaShape.s:{run_name}"/>
            <log idref="freqParameter.s:{run_name}"/>
            <log idref="bPopSizes.t:{run_name}"/>
            <log idref="bGroupSizes.t:{run_name}"/>
            <log id="TreeHeight.t:{run_name}" spec="beast.evolution.tree.TreeHeightLogger"
                 tree="@Tree.t:{run_name}"/>
            <log idref="TreeHeight.t:{run_name}"/>
        </logger>

        <logger id="treelog" fileName="{run_name}.trees" logEvery="15000"
                mode="tree">
            <log id="treeWithMetaDataLogger.t:{run_name}"
                 spec="beast.evolution.tree.TreeWithMetaDataLogger"
                 tree="@Tree.t:{run_name}"/>
        </logger>

        <logger id="screenlog" logEvery="100000">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log id="ESS.0" spec="util.ESS" arg="@posterior" name="log"/>
        </logger>

    </run>
</beast>
"""
    with open(out_xml, 'w') as f:
        f.write(xml)
    print(f"  Written: {out_xml} ({n_seq} sequences, {aln_len} bp, {n_groups} BSP groups)")

# Generate XML for each population + combined
for fa in sorted(os.listdir(bsp_dir)):
    if not fa.endswith('.fa'):
        continue
    label    = fa.replace('CYTB_', '').replace('.fa', '')
    run_name = f"BSP_{label}"
    out_xml  = os.path.join(bsp_dir, f"{run_name}.xml")
    make_beast2_xml(os.path.join(bsp_dir, fa), run_name, out_xml, clock_rate_min, clock_rate_max, species)

print("All XML files generated.")
PYEOF

ok "BEAST2 XML files generated → ${BSP_DIR}"

# ── Step 3: Run BEAST2 for each XML (PARALLEL) ────────────────────────────────
log "Running BEAST2..."

BEAST_PIDS=()
BEAST_RUNS=()

for XML in "${BSP_DIR}"/BSP_*.xml; do
    RUN=$(basename "$XML" .xml)
    LOG="${BSP_DIR}/${RUN}.log"
    STDERR_LOG="${BSP_DIR}/${RUN}.beast.stderr"

    log "  Launching BEAST2: ${RUN}..."
    (cd "${BSP_DIR}" && "${BEAST_BIN}" -threads 4 -overwrite "${RUN}.xml" \
        > "${STDERR_LOG}" 2>&1) &
    BEAST_PIDS+=($!)
    BEAST_RUNS+=("$RUN")
done

# Wait for all background jobs, report individually
ALL_OK=true
for i in "${!BEAST_PIDS[@]}"; do
    PID="${BEAST_PIDS[$i]}"
    RUN="${BEAST_RUNS[$i]}"
    if wait "$PID"; then
        ok "  ${RUN} complete"
    else
        warn "  ${RUN} BEAST2 FAILED — check ${BSP_DIR}/${RUN}.beast.stderr"
        ALL_OK=false
    fi
done

$ALL_OK && ok "All BEAST2 runs complete" || warn "One or more BEAST2 runs failed — check stderr logs above"

# ── Step 4: Check convergence with Tracer (headless ESS check) ───────────────
log "Checking ESS from BEAST2 log files..."

python3 - "$BSP_DIR" <<'EOF'
import sys
import os
import csv
import math

bsp_dir = sys.argv[1]

def compute_ess(values):
    """Compute ESS using autocorrelation method."""
    n = len(values)
    if n < 10:
        return 0
    mean = sum(values) / n
    var  = sum((x - mean)**2 for x in values) / n
    if var == 0:
        return 0
    # Autocorrelation up to lag 100 or n//2
    max_lag = min(100, n // 2)
    rho_sum = 0.0
    for lag in range(1, max_lag + 1):
        rho = sum((values[i] - mean) * (values[i + lag] - mean)
                  for i in range(n - lag)) / ((n - lag) * var)
        if rho < 0.05:
            break
        rho_sum += rho
    ess = n / (1 + 2 * rho_sum)
    return ess

ess_results = []
for fname in sorted(os.listdir(bsp_dir)):
    if not (fname.startswith("BSP_") and fname.endswith(".log")):
        continue
    run_dir  = fname[:-4]          # strip .log
    log_file = os.path.join(bsp_dir, fname)
    # Read log, skip comments
    headers = None
    rows    = []
    with open(log_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if headers is None:
                headers = parts
            else:
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    pass
    if not headers or not rows:
        continue
    # Discard first 10% as burnin
    burnin = len(rows) // 10
    rows   = rows[burnin:]
    # Compute ESS for posterior (column index 1 typically)
    post_idx = headers.index('posterior') if 'posterior' in headers else 1
    post_vals = [r[post_idx] for r in rows]
    ess = compute_ess(post_vals)
    status = "OK" if ess >= 200 else "LOW — consider longer chain"
    ess_results.append({'run': run_dir, 'ESS_posterior': f"{ess:.0f}", 'status': status})
    print(f"  {run_dir}: ESS(posterior) = {ess:.0f} [{status}]")

# Write ESS summary
out_tsv = os.path.join(bsp_dir, 'ess_summary.tsv')
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['run', 'ESS_posterior', 'status'], delimiter='\t')
    w.writeheader()
    w.writerows(ess_results)
print(f"ESS summary written: {out_tsv}")
EOF

ok "ESS check complete → ${BSP_DIR}/ess_summary.tsv"

# ── Step 5: Run TreeAnnotator for each run ────────────────────────────────────
log "Running TreeAnnotator (burnin 10%)..."

for TREES in "${BSP_DIR}"/BSP_*.trees; do
    [[ -f "$TREES" ]] || continue
    RUN=$(basename "$TREES" .trees)
    MCC="${BSP_DIR}/${RUN}_MCC.tree"

    log "  TreeAnnotator: ${RUN} (burnin 10%)..."
    "${BEAST_ENV}/bin/treeannotator" -burnin 10 -heights mean \
        "$TREES" "$MCC" > "${BSP_DIR}/${RUN}.treeannotator.log" 2>&1 \
        && ok "  ${RUN} MCC tree written" \
        || warn "  TreeAnnotator failed for ${RUN}"
done

ok "TreeAnnotator complete"

# ── Step 6: Extract BSP skyline data for plotting ─────────────────────────────
log "Extracting BSP skyline data for R plotting..."

python3 - "$BSP_DIR" "$SPECIES" <<'EOF'
import sys
import os
import csv
import re

bsp_dir = sys.argv[1]
species = sys.argv[2] if len(sys.argv) > 2 else "Unknown species"

# Write an R script to plot all BSPs using the boa/coda approach
# User runs this manually in R after the pipeline completes
r_script = """
# =============================================================================
# BSP plotting script — generated by m09_bsp.sh
# Requires: ggplot2, reshape2, tracerer or beastio
# Run in R after BEAST2 completes
# =============================================================================

library(ggplot2)

# If you have the 'beastio' or 'phylodyn' package:
# install.packages('devtools')
# devtools::install_github('laduplessis/beastio')

# Manual approach using BEAST2 log files:
read_bsp_log <- function(log_file, burnin_frac = 0.1) {
  dat <- read.table(log_file, header=TRUE, comment.char='#', sep='\\t')
  burnin <- floor(nrow(dat) * burnin_frac)
  dat <- dat[(burnin+1):nrow(dat), ]
  return(dat)
}

bsp_dir <- "{bsp_dir}"
species_name <- "{species_name}"
runs <- list.dirs(bsp_dir, recursive=FALSE, full.names=FALSE)
runs <- runs[grepl("^BSP_", runs)]

plots <- list()
for (run in runs) {
  log_file <- file.path(bsp_dir, paste0(run, ".log"))
  if (!file.exists(log_file)) next
  dat  <- read_bsp_log(log_file)
  pop_cols <- grep("^bPopSizes", names(dat), value=TRUE)
  if (length(pop_cols) == 0) next
  pop_means <- colMeans(dat[, pop_cols, drop=FALSE])
  pop_lower <- apply(dat[, pop_cols, drop=FALSE], 2, quantile, 0.025)
  pop_upper <- apply(dat[, pop_cols, drop=FALSE], 2, quantile, 0.975)
  df <- data.frame(
    group  = seq_along(pop_means),
    mean   = pop_means,
    lower  = pop_lower,
    upper  = pop_upper,
    run    = gsub("BSP_", "", run)
  )
  plots[[run]] <- df
}

all_df <- do.call(rbind, plots)

p <- ggplot(all_df, aes(x=group, y=mean, ymin=lower, ymax=upper, colour=run, fill=run)) +
  geom_ribbon(alpha=0.2, colour=NA) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~run, scales='free_y') +
  labs(x="Coalescent interval", y="Effective population size (log scale)",
       title=paste0("Bayesian Skyline Plot \u2014 ", species_name, " CYTB")) +
  theme_bw() +
  theme(legend.position='none')

ggsave(file.path(bsp_dir, "BSP_all_populations.pdf"), p, width=14, height=10)
message("BSP plot saved: ", file.path(bsp_dir, "BSP_all_populations.pdf"))
"""

r_script = r_script.replace("{bsp_dir}", bsp_dir)
r_script = r_script.replace("{species_name}", species)
r_out = os.path.join(bsp_dir, "plot_bsp.R")
with open(r_out, "w") as f:
    f.write(r_script)
print(f"R plotting script written: {r_out}")
EOF

ok "R plotting script written → ${BSP_DIR}/plot_bsp.R"
log "Run manually in R: source('${BSP_DIR}/plot_bsp.R')"