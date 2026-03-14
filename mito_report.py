#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime
from io import BytesIO
import base64

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import networkx as nx

from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate mitogenome HTML report")
parser.add_argument("--outdir",       required=True)
parser.add_argument("--species",      default="Unknown species")
parser.add_argument("--output",       default="mito_report.html")
parser.add_argument("--min-pop-size", type=int, default=0,
                    help="Minimum samples per population for per-pop analyses. "
                         "0 = include all (default). Excluded pops shown as greyed "
                         "out in tables and omitted from heatmap/BSP figures.")
parser.add_argument("--genetic-code", type=int, default=2,
                    help="NCBI genetic code used (for methods text)")
parser.add_argument("--clade",        default="Chordata",
                    help="MitoZ clade used (for methods text)")
parser.add_argument("--clock-rate-min", type=float, default=0.005,
                    help="BEAST2 clock rate prior lower bound (for methods text)")
parser.add_argument("--clock-rate-max", type=float, default=0.02,
                    help="BEAST2 clock rate prior upper bound (for methods text)")
args = parser.parse_args()

OUTDIR        = args.outdir
SPECIES       = args.species
OUTHTML       = args.output
MIN_POP_SIZE  = args.min_pop_size
GENETIC_CODE  = args.genetic_code
CLADE         = args.clade
CLOCK_MIN     = args.clock_rate_min
CLOCK_MAX     = args.clock_rate_max

# Determine if this is a filtered report
IS_FILTERED   = MIN_POP_SIZE > 0

# ── File paths ─────────────────────────────────────────────────────────────────
SUMMARY_TSV        = os.path.join(OUTDIR, "qc",     "pipeline_summary.tsv")
CONSENSUS_FA       = os.path.join(OUTDIR, "aligned", "consensus_mito.fa")
ALL_MITO_FA        = os.path.join(OUTDIR, "aligned", "all_samples_mito_aligned.fa")
TREEFILE           = os.path.join(OUTDIR, "tree",    "mito_tree.treefile")
YN00_TSV           = os.path.join(OUTDIR, "dnds",   "yn00_summary.tsv")
CODEML_TSV         = os.path.join(OUTDIR, "dnds",   "codeml_branch_summary.tsv")
CODEML_EXCL_TSV    = os.path.join(OUTDIR, "dnds",   "codeml_excluded_pops.tsv")
BSP_ESS_TSV        = os.path.join(OUTDIR, "bsp",    "ess_summary.tsv")
BSP_POP_JSON       = os.path.join(OUTDIR, "bsp",    "bsp_pop_summary.json")
CODEML_POP_JSON    = os.path.join(OUTDIR, "dnds",   "codeml", "pop_sizes.json")
ANNOTATION_SUMMARY = os.path.join(OUTDIR, "annotation_test",
                     "test.consensus_mito.fa.result", "summary.txt")

# ── Helpers ───────────────────────────────────────────────────────────────────
def median(lst):
    if not lst: return 0
    lst = sorted(lst)
    n = len(lst)
    return lst[n//2] if n % 2 else (lst[n//2-1] + lst[n//2]) / 2

def safe_float(x):
    try: return float(x)
    except: return None

def read_fasta_seq(path):
    seq = ""
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                if not line.startswith(">"): seq += line.strip()
    return seq

def encode_fig(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def get_pop(name):
    m = re.search(r'-([A-Z]+)\d', name)
    return m.group(1) if m else "Unknown"

def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='	'))

def load_json(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)

# ── Load pop size metadata ─────────────────────────────────────────────────────
bsp_pop_summary    = load_json(BSP_POP_JSON)
codeml_pop_summary = load_json(CODEML_POP_JSON)

# Determine which pops were excluded by each analysis
bsp_excluded    = set(bsp_pop_summary.get('excluded_pops', {}).keys())
codeml_excluded = set(codeml_pop_summary.get('excluded_pops', {}).keys())
all_pop_sizes   = bsp_pop_summary.get('all_pops', codeml_pop_summary.get('all_pops', {}))

# For the filtered report, also apply the --min-pop-size filter at report level
# (in case report is run with a different threshold than the pipeline)
if IS_FILTERED:
    report_excluded = {p for p, n in all_pop_sizes.items() if n < MIN_POP_SIZE}
else:
    report_excluded = set()

# ── Load QC summary ───────────────────────────────────────────────────────────
samples = []
with open(SUMMARY_TSV) as f:
    header = f.readline().strip().split("	")
    for line in f:
        parts = line.strip().split("	")
        if len(parts) < 2 or parts[0] in ("", "0"): continue
        samples.append(dict(zip(header, parts)))

n_total = len(samples)
n_pass  = sum(1 for s in samples if "PASS" in s.get("status",""))
n_warn  = sum(1 for s in samples if "WARN" in s.get("status",""))
n_fail  = sum(1 for s in samples if "FAIL" in s.get("status",""))

depths, breadths, mapped_reads = [], [], []
for s in samples:
    d = safe_float(s.get("mean_depth"))
    b = safe_float(s.get("coverage_breadth_pct"))
    m = safe_float(s.get("mapped_reads"))
    if d is not None: depths.append(d)
    if b is not None: breadths.append(b)
    if m is not None: mapped_reads.append(m)

# ── Consensus ─────────────────────────────────────────────────────────────────
consensus   = read_fasta_seq(CONSENSUS_FA)
cons_len    = len(consensus)
base_counts = Counter(consensus.upper())
gc = (base_counts.get("G",0) + base_counts.get("C",0)) / cons_len * 100 if cons_len else 0

# ── Load alignment + population codes ─────────────────────────────────────────
print("Loading alignment...")
records = list(SeqIO.parse(ALL_MITO_FA, "fasta"))

seqs_by_pop = defaultdict(list)
for rec in records:
    pop = get_pop(rec.id)
    seqs_by_pop[pop].append(str(rec.seq))

# All pops and filtered pops
ALL_POP_ORDER      = sorted(seqs_by_pop.keys())
FILTERED_POP_ORDER = sorted(p for p in ALL_POP_ORDER if p not in report_excluded)
POP_ORDER          = FILTERED_POP_ORDER if IS_FILTERED else ALL_POP_ORDER
N_POPS             = len(POP_ORDER)

print(f"  All populations:      {ALL_POP_ORDER}")
if IS_FILTERED:
    print(f"  Filtered populations: {FILTERED_POP_ORDER}")
    print(f"  Excluded:             {sorted(report_excluded)}")

cmap = matplotlib.colormaps.get_cmap("tab10").resampled(max(len(ALL_POP_ORDER), 1))
pop_colors_map = {pop: cmap(i) for i, pop in enumerate(ALL_POP_ORDER)}

# ── Load annotation summary ───────────────────────────────────────────────────
print("Loading annotation data...")
annot_genes   = []
annot_stats   = {}
annot_related = ""
if os.path.exists(ANNOTATION_SUMMARY):
    with open(ANNOTATION_SUMMARY) as f:
        content = f.read()
    for line in content.split('
'):
        parts = line.split()
        if len(parts) >= 9 and re.match(r'\S+_mito', parts[0]):
            annot_genes.append({
                'seq_id':    parts[0],
                'start':     parts[1],
                'end':       parts[2],
                'length_bp': parts[3],
                'direction': parts[4],
                'type':      parts[5],
                'gene':      parts[6],
                'product':   ' '.join(parts[7:-1]),
                'freq':      parts[-1],
            })
    for line in content.split('
'):
        m = re.match(r'(Protein coding|tRNA|rRNA) genes totally found:\s+(\d+)', line)
        if m:
            annot_stats[m.group(1)] = int(m.group(2))
    m2 = re.search(r'Closely_related_species\s*
\S+\s+\d+\s+\S+\s+(.+)', content)
    if m2:
        annot_related = m2.group(1).strip()

# ── Load dN/dS results ────────────────────────────────────────────────────────
print("Loading dN/dS results...")
yn00_rows   = read_tsv(YN00_TSV)
codeml_rows_all = read_tsv(CODEML_TSV)

# Filter codeml rows for filtered report
if IS_FILTERED:
    codeml_rows = [r for r in codeml_rows_all
                   if r.get('population','') not in report_excluded]
else:
    codeml_rows = codeml_rows_all

# Collect excluded pop info for display
excl_rows = read_tsv(CODEML_EXCL_TSV)
# Merge pipeline exclusions with report-level exclusions
excl_display = {}
for r in excl_rows:
    excl_display[r['population']] = r['n_samples']
if IS_FILTERED:
    for p in report_excluded:
        if p not in excl_display:
            excl_display[p] = all_pop_sizes.get(p, '?')

# ── Load BSP ESS summary ──────────────────────────────────────────────────────
print("Loading BSP ESS data...")
bsp_ess_rows_all = read_tsv(BSP_ESS_TSV)

if IS_FILTERED:
    bsp_ess_rows = [r for r in bsp_ess_rows_all
                    if r.get('run','').replace('BSP_','') not in report_excluded
                    or r.get('run','') == 'BSP_all']
else:
    bsp_ess_rows = bsp_ess_rows_all

# ── Figure 1: Depth histogram ─────────────────────────────────────────────────
print("Generating figures...")
fig1, ax1 = plt.subplots(figsize=(7, 4))
ax1.hist(depths, bins=30, color="#2196F3", edgecolor="white", linewidth=0.5)
ax1.axvline(median(depths), color="#F44336", linestyle="--", linewidth=1.5,
            label=f"Median: {median(depths):.0f}x")
ax1.set_xlabel("Mean depth (x)", fontsize=11)
ax1.set_ylabel("Number of samples", fontsize=11)
ax1.set_title("Coverage Depth Distribution", fontsize=12, fontweight="bold")
ax1.legend()
ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)
depth_hist_b64 = encode_fig(fig1)

# ── Figure 2: Breadth histogram ───────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.hist(breadths, bins=30, color="#4CAF50", edgecolor="white", linewidth=0.5)
ax2.axvline(median(breadths), color="#F44336", linestyle="--", linewidth=1.5,
            label=f"Median: {median(breadths):.1f}%")
ax2.set_xlabel("Coverage breadth (%)", fontsize=11)
ax2.set_ylabel("Number of samples", fontsize=11)
ax2.set_title("Genome Coverage Breadth", fontsize=12, fontweight="bold")
ax2.legend()
ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
breadth_hist_b64 = encode_fig(fig2)

# ── Figure 3: Scatter mapped reads vs depth ───────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(7, 4))
status_colors = []
for s in samples:
    st = s.get("status","")
    if "PASS" in st:   status_colors.append("#4CAF50")
    elif "WARN" in st: status_colors.append("#FF9800")
    else:              status_colors.append("#F44336")
ax3.scatter(mapped_reads, depths, c=status_colors, alpha=0.6, s=30, edgecolors="none")
ax3.set_xlabel("Mapped reads", fontsize=11)
ax3.set_ylabel("Mean depth (x)", fontsize=11)
ax3.set_title("Mapped Reads vs Coverage Depth", fontsize=12, fontweight="bold")
ax3.legend(handles=[
    mpatches.Patch(color="#4CAF50", label="PASS"),
    mpatches.Patch(color="#FF9800", label="WARN"),
    mpatches.Patch(color="#F44336", label="FAIL"),
])
ax3.spines["top"].set_visible(False); ax3.spines["right"].set_visible(False)
scatter_b64 = encode_fig(fig3)

# ── Figure 4: Status pie ──────────────────────────────────────────────────────
fig4, ax4 = plt.subplots(figsize=(5, 4))
pie_vals   = [x for x in [n_pass, n_warn, n_fail] if x > 0]
pie_labels = [f"{l}
({v})" for l, v in zip(["PASS","WARN","FAIL"],
              [n_pass, n_warn, n_fail]) if v > 0]
pie_colors = [c for c, v in zip(["#4CAF50","#FF9800","#F44336"],
              [n_pass, n_warn, n_fail]) if v > 0]
ax4.pie(pie_vals, labels=pie_labels, colors=pie_colors,
        startangle=90, textprops={"fontsize": 10})
ax4.set_title("QC Status", fontsize=12, fontweight="bold")
pie_b64 = encode_fig(fig4)

# ── Figure 5: Base composition ────────────────────────────────────────────────
fig5, ax5 = plt.subplots(figsize=(5, 4))
bases     = ["A","T","G","C","N"]
base_pcts = [base_counts.get(b,0)/cons_len*100 for b in bases]
bar_cols  = ["#2196F3","#F44336","#4CAF50","#FF9800","#9E9E9E"]
bars = ax5.bar(bases, base_pcts, color=bar_cols, edgecolor="white")
for bar, pct in zip(bars, base_pcts):
    if pct > 0.5:
        ax5.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.3,
                 f"{pct:.1f}%", ha="center", va="bottom", fontsize=9)
ax5.set_xlabel("Base"); ax5.set_ylabel("Percentage (%)")
ax5.set_title("Consensus Base Composition", fontsize=12, fontweight="bold")
ax5.spines["top"].set_visible(False); ax5.spines["right"].set_visible(False)
basecomp_b64 = encode_fig(fig5)

# ── Figure 6: Population NJ tree ─────────────────────────────────────────────
print("Building population NJ tree...")
tree_b64 = None
try:
    pop_cons_records = []
    for pop in POP_ORDER:
        seqs   = seqs_by_pop[pop]
        length = len(seqs[0])
        cons_seq = []
        for i in range(length):
            col = [s[i] for s in seqs if i < len(s) and s[i] not in "-N"]
            cons_seq.append(Counter(col).most_common(1)[0][0] if col else "N")
        pop_cons_records.append(SeqRecord(Seq("".join(cons_seq)), id=pop, description=""))

    tmp_fa = os.path.join(OUTDIR, "tmp_pop_consensus.fa")
    SeqIO.write(pop_cons_records, tmp_fa, "fasta")
    alignment   = AlignIO.read(tmp_fa, "fasta")
    calculator  = DistanceCalculator("identity")
    dm          = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    pop_tree    = constructor.nj(dm)
    os.remove(tmp_fa)

    fig6, ax6 = plt.subplots(figsize=(10, 6))
    Phylo.draw(pop_tree, do_show=False, axes=ax6,
               label_func=lambda x: x.name if x.is_terminal() else "",
               label_colors={pop: pop_colors_map[pop] for pop in POP_ORDER})
    for text_obj in ax6.texts:
        label = text_obj.get_text().strip()
        if label in seqs_by_pop:
            text_obj.set_text(f"{label}  (n={len(seqs_by_pop[label])})")
    ax6.set_title(f"Population-Level Mitogenome NJ Tree
{SPECIES}",
                  fontsize=13, fontweight="bold")
    ax6.spines["top"].set_visible(False); ax6.spines["right"].set_visible(False)
    ax6.legend(handles=[
        mpatches.Patch(color=pop_colors_map[p],
                       label=f"{p}  (n={len(seqs_by_pop[p])})")
        for p in POP_ORDER
    ], title="Population", loc="lower right", fontsize=9)
    plt.subplots_adjust(left=0.05, right=0.75, top=0.92, bottom=0.05)
    tree_b64 = encode_fig(fig6)
    print("  Population tree done")
except Exception as e:
    print(f"  Tree failed: {e}")

# ── Figure 7: Population network ──────────────────────────────────────────────
print("Building population network...")
haplo_b64 = None
try:
    seq_length = len(records[0].seq)
    dist_between_pops = {}
    for i, p1 in enumerate(POP_ORDER):
        for j, p2 in enumerate(POP_ORDER):
            if i >= j: continue
            dists = []
            for s1 in seqs_by_pop[p1]:
                for s2 in seqs_by_pop[p2]:
                    d = sum(a != b and a not in "-N" and b not in "-N"
                            for a, b in zip(s1, s2)) / seq_length
                    dists.append(d)
            dist_between_pops[(p1, p2)] = np.mean(dists) if dists else 1.0

    G = nx.Graph()
    for pop in POP_ORDER:
        G.add_node(pop, size=len(seqs_by_pop[pop]))
    for (p1, p2), dist in dist_between_pops.items():
        G.add_edge(p1, p2, weight=dist)

    pos      = nx.spring_layout(G, weight="weight", seed=42, k=3)
    max_dist = max(dist_between_pops.values()) if dist_between_pops else 1
    edge_widths = [max(0.5, (1 - G[u][v]["weight"] / max_dist) * 6)
                   for u, v in G.edges()]

    fig7, ax7 = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, ax=ax7, width=edge_widths,
                           alpha=0.5, edge_color="#888888")
    nx.draw_networkx_nodes(G, pos, ax=ax7,
                           node_color=[pop_colors_map[p] for p in G.nodes()],
                           node_size=[len(seqs_by_pop[p])*60+200 for p in G.nodes()],
                           alpha=0.9)
    nx.draw_networkx_labels(G, pos, ax=ax7,
                            labels={p: f"{p}
(n={len(seqs_by_pop[p])})"
                                    for p in G.nodes()},
                            font_size=9, font_weight="bold")
    ax7.set_title(f"Population Network — {SPECIES}
"
                  "Node size = sample count, edge thickness = genetic similarity",
                  fontsize=12, fontweight="bold")
    ax7.axis("off")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.02)
    haplo_b64 = encode_fig(fig7)
    print("  Population network done")
except Exception as e:
    print(f"  Network failed: {e}")

# ── Figure 8: dN/dS per gene bar chart ───────────────────────────────────────
print("Generating dN/dS figure...")
dnds_fig_b64 = None
if yn00_rows:
    try:
        genes     = [r['gene'] for r in yn00_rows]
        mean_dnds = [float(r['mean_dNdS']) for r in yn00_rows]
        colors    = ['#F44336' if v > 1 else '#2196F3' for v in mean_dnds]

        fig8, ax8 = plt.subplots(figsize=(12, 5))
        bars = ax8.bar(genes, mean_dnds, color=colors, edgecolor='white')
        ax8.axhline(1.0, color='black', linestyle='--', linewidth=1,
                    label='dN/dS = 1 (neutral)')
        ax8.axhline(0.0, color='grey', linestyle='-', linewidth=0.5)
        for bar, val in zip(bars, mean_dnds):
            ax8.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                     f"{val:.3f}", ha='center', va='bottom', fontsize=8)
        ax8.set_xlabel("Gene", fontsize=11)
        ax8.set_ylabel("Mean pairwise dN/dS (yn00)", fontsize=11)
        ax8.set_title(f"Pairwise dN/dS per Mitochondrial Gene — {SPECIES}",
                      fontsize=12, fontweight="bold")
        ax8.legend()
        ax8.spines["top"].set_visible(False); ax8.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        dnds_fig_b64 = encode_fig(fig8)
        print("  dN/dS figure done")
    except Exception as e:
        print(f"  dN/dS figure failed: {e}")

# ── Figure 9: codeml heatmap ──────────────────────────────────────────────────
print("Generating codeml heatmap...")
codeml_fig_b64 = None
if codeml_rows:
    try:
        genes = sorted(set(r['gene'] for r in codeml_rows))
        pops  = sorted(set(r['population'] for r in codeml_rows))

        omega_mat = np.full((len(pops), len(genes)), np.nan)
        sig_mat   = np.zeros((len(pops), len(genes)), dtype=bool)
        for r in codeml_rows:
            gi = genes.index(r['gene'])
            pi = pops.index(r['population'])
            try:
                omega_mat[pi, gi] = float(r['omega_foreground'])
                sig_mat[pi, gi]   = r['p<0.05'] == 'YES'
            except (ValueError, KeyError):
                pass

        fig9, ax9 = plt.subplots(figsize=(max(10, len(genes)*1.2),
                                          max(5, len(pops)*0.8)))
        im = ax9.imshow(omega_mat, aspect='auto', cmap='RdYlBu_r', vmin=0, vmax=1)
        plt.colorbar(im, ax=ax9, label='Foreground ω (dN/dS)')
        ax9.set_xticks(range(len(genes)))
        ax9.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
        ax9.set_yticks(range(len(pops)))
        ax9.set_yticklabels(pops, fontsize=9)
        ax9.set_title(f"codeml Branch Model — Foreground ω per Population × Gene
"
                      f"* = significant (LRT p<0.05)", fontsize=11, fontweight="bold")
        for pi in range(len(pops)):
            for gi in range(len(genes)):
                if sig_mat[pi, gi]:
                    ax9.text(gi, pi, '*', ha='center', va='center',
                             fontsize=14, fontweight='bold', color='white')
        plt.tight_layout()
        codeml_fig_b64 = encode_fig(fig9)
        print("  codeml heatmap done")
    except Exception as e:
        print(f"  codeml heatmap failed: {e}")

# ── Figure 9: BSP population size plots ──────────────────────────────────────
print("Generating BSP figures...")
bsp_plots_b64 = {}
if bsp_pop_data:
    POP_COLOURS = {"SI":"#1976D2","NS":"#388E3C","CAI":"#F57C00",
                   "SB":"#7B1FA2","COI":"#D32F2F","SP":"#0097A7",
                   "DH":"#5D4037","all":"#455A64"}
    for pop_label, data in bsp_pop_data.items():
        pop_sizes = data["pop_sizes"]
        if not pop_sizes: continue
        n_groups = len(pop_sizes[0])
        import numpy as np
        ps_arr = np.array(pop_sizes)  # shape: (n_samples, n_groups)
        medians = np.median(ps_arr, axis=0)
        lo95    = np.percentile(ps_arr, 2.5,  axis=0)
        hi95    = np.percentile(ps_arr, 97.5, axis=0)
        groups  = np.arange(1, n_groups + 1)
        colour  = POP_COLOURS.get(pop_label, "#607D8B")
        try:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.fill_between(groups, lo95, hi95, alpha=0.25, color=colour)
            ax.plot(groups, medians, color=colour, linewidth=2)
            ax.set_yscale("log")
            ax.set_xlabel("Coalescent interval (recent → ancient)", fontsize=10)
            ax.set_ylabel("Effective population size (Ne, log scale)", fontsize=10)
            n_samples = bsp_pop_summary.get("included_pops", {}).get(pop_label,
                        bsp_pop_summary.get("all_pops", {}).get(pop_label, ""))
            title_pop = "All samples" if pop_label == "all" else f"Population {pop_label}"
            n_label   = f" (n={n_samples})" if n_samples else ""
            ax.set_title(f"BSP — {title_pop}{n_label}", fontsize=11)
            ax.tick_params(labelsize=9)
            fig.tight_layout()
            bsp_plots_b64[pop_label] = encode_fig(fig)
            plt.close(fig)
            print(f"  BSP plot done for {pop_label}")
        except Exception as e:
            print(f"  BSP plot failed for {pop_label}: {e}")

# ── Figure 10: BSP ESS bar chart ─────────────────────────────────────────────
print("Generating BSP ESS figure...")
bsp_ess_fig_b64 = None
if bsp_ess_rows:
    try:
        runs = [r['run'].replace('BSP_','') for r in bsp_ess_rows]
        ess  = [float(r['ESS_posterior']) for r in bsp_ess_rows]
        cols = ['#4CAF50' if e >= 200 else '#F44336' for e in ess]

        fig10, ax10 = plt.subplots(figsize=(10, 4))
        bars = ax10.bar(runs, ess, color=cols, edgecolor='white')
        ax10.axhline(200, color='black', linestyle='--', linewidth=1,
                     label='ESS = 200 (minimum threshold)')
        for bar, val in zip(bars, ess):
            ax10.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                      f"{val:.0f}", ha='center', va='bottom', fontsize=8)
        ax10.set_xlabel("Run", fontsize=11)
        ax10.set_ylabel("ESS (posterior)", fontsize=11)
        ax10.set_title("BEAST2 ESS — Convergence Check per BSP Run",
                       fontsize=12, fontweight="bold")
        ax10.legend()
        ax10.spines["top"].set_visible(False); ax10.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        bsp_ess_fig_b64 = encode_fig(fig10)
        print("  BSP ESS figure done")
    except Exception as e:
        print(f"  BSP ESS figure failed: {e}")

# ── Build excluded-population banner ─────────────────────────────────────────
excl_banner = ""
if excl_display:
    excl_items = ", ".join(f"<b>{p}</b> (n={n})" for p, n in sorted(excl_display.items()))
    if IS_FILTERED:
        excl_banner = f"""
        <div style="background:#FFF8E1;border-left:4px solid #FF9800;padding:12px 16px;
                    border-radius:4px;margin-bottom:16px;font-size:0.9em;">
          ⚠ <b>Population filter applied (n &ge; {MIN_POP_SIZE}).</b>
          The following populations were excluded from per-population analyses
          (codeml branch model, BSP): {excl_items}.
          They are still included in pairwise yn00 dN/dS and the combined BSP run.
          See <code>mito_report_all.html</code> for the unfiltered report.
        </div>"""
    else:
        excl_banner = f"""
        <div style="background:#E3F2FD;border-left:4px solid #1565C0;padding:12px 16px;
                    border-radius:4px;margin-bottom:16px;font-size:0.9em;">
          ℹ <b>Small populations detected.</b>
          The following populations have fewer samples than the pipeline filter threshold
          and were excluded from per-population analyses during pipeline execution: {excl_items}.
          See <code>mito_report_filtered.html</code> for the report with these excluded.
        </div>"""

# ── Table: excluded populations ───────────────────────────────────────────────
excl_table_rows = ""
for pop, n in sorted(excl_display.items()):
    excl_table_rows += f"""
    <tr>
        <td><b>{pop}</b></td>
        <td>{n}</td>
        <td>Below min_pop_size threshold — excluded from codeml branch model and BSP</td>
    </tr>"""

# ── QC table rows ─────────────────────────────────────────────────────────────
qc_rows = ""
for s in sorted(samples, key=lambda x: x.get("sample","")):
    status = s.get("status","")
    css = "pass" if "PASS" in status else ("warn" if "WARN" in status else "fail")
    qc_rows += f"""
    <tr class="{css}">
        <td>{s.get('sample','')}</td>
        <td>{int(float(s.get('mapped_reads',0))):,}</td>
        <td>{float(s.get('mean_depth',0)):.1f}x</td>
        <td>{float(s.get('coverage_breadth_pct',0)):.1f}%</td>
        <td>{s.get('num_variants_pass','')}</td>
        <td>{int(s.get('consensus_length_bp',0)):,} bp</td>
        <td>{s.get('N_count','0')}</td>
        <td><b>{status}</b></td>
    </tr>"""

# ── Annotation table rows ─────────────────────────────────────────────────────
annot_rows = ""
for g in annot_genes:
    type_col = "#2196F3" if g['type'] == 'CDS' else "#4CAF50"
    annot_rows += f"""
    <tr>
        <td><b>{g['gene']}</b></td>
        <td>{g['product']}</td>
        <td><span style="color:{type_col};font-weight:bold">{g['type']}</span></td>
        <td>{g['start']}–{g['end']}</td>
        <td>{g['length_bp']} bp</td>
        <td>{'&minus;' if g['direction'] == '-' else '+'}</td>
    </tr>"""

# ── yn00 table rows ───────────────────────────────────────────────────────────
yn00_table_rows = ""
for r in yn00_rows:
    try:
        omega = float(r['mean_dNdS'])
        col = "#F44336" if omega > 1 else ("#FF9800" if omega > 0.5 else "#4CAF50")
        yn00_table_rows += f"""
        <tr>
            <td><b>{r['gene']}</b></td>
            <td>{r['n_pairs']}</td>
            <td>{float(r['mean_dN']):.5f}</td>
            <td>{float(r['mean_dS']):.5f}</td>
            <td style="color:{col};font-weight:bold">{omega:.5f}</td>
            <td>{float(r['min_dNdS']):.5f}</td>
            <td>{float(r['max_dNdS']):.5f}</td>
        </tr>"""
    except (ValueError, KeyError):
        pass

# ── codeml table rows ─────────────────────────────────────────────────────────
codeml_table_rows = ""
for r in codeml_rows:
    sig   = r.get('p<0.05','no')
    css   = 'pass' if sig == 'YES' else ''
    badge = "<span style='color:#F44336;font-weight:bold'>YES ✱</span>"             if sig == 'YES' else "no"
    n_s   = r.get('n_samples','')
    codeml_table_rows += f"""
    <tr class="{css}">
        <td><b>{r.get('gene','')}</b></td>
        <td>{r.get('population','')} <span style="color:#999;font-size:0.85em">(n={n_s})</span></td>
        <td>{r.get('lnL_null','')}</td>
        <td>{r.get('lnL_branch','')}</td>
        <td>{r.get('LRT','')}</td>
        <td>{badge}</td>
        <td>{r.get('omega_background','')}</td>
        <td>{r.get('omega_foreground','')}</td>
    </tr>"""

# ── BSP ESS table rows ────────────────────────────────────────────────────────
bsp_ess_table_rows = ""
for r in bsp_ess_rows:
    css = 'pass' if r.get('status','').startswith('OK') else 'fail'
    bsp_ess_table_rows += f"""
    <tr class="{css}">
        <td>{r.get('run','').replace('BSP_','')}</td>
        <td>{r.get('ESS_posterior','')}</td>
        <td><b>{r.get('status','')}</b></td>
    </tr>"""

# ── Methods text (dynamic) ────────────────────────────────────────────────────
# Genetic code description
gc_descriptions = {
    1: "standard", 2: "vertebrate mitochondrial", 4: "mold/protozoan mitochondrial",
    5: "invertebrate mitochondrial", 9: "echinoderm mitochondrial",
    13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial"
}
gc_desc = gc_descriptions.get(GENETIC_CODE, f"code {GENETIC_CODE}")

# Clock rate description
clock_desc = f"U({CLOCK_MIN}, {CLOCK_MAX}) subs/site/Myr"

# Filtered report subtitle
report_subtitle = ""
if IS_FILTERED:
    report_subtitle = f' <span style="background:#FFF8E1;color:#E65100;padding:2px 8px;border-radius:10px;font-size:0.7em;font-weight:bold;">FILTERED (n≥{MIN_POP_SIZE})</span>'
else:
    report_subtitle = f' <span style="background:#E8F5E9;color:#2E7D32;padding:2px 8px;border-radius:10px;font-size:0.7em;font-weight:bold;">ALL POPULATIONS</span>'

# ── Write HTML ────────────────────────────────────────────────────────────────
print("Writing HTML...")
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Mitogenome Report — {SPECIES}{"  [filtered]" if IS_FILTERED else ""}</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 0;
          background: #f4f6f9; color: #222; }}
  .header {{ background: linear-gradient(135deg, #1A237E, #3949AB);
             color: white; padding: 40px; text-align: center; }}
  .header h1 {{ margin: 0; font-size: 2em; }}
  .header p  {{ margin: 8px 0 0; opacity: 0.85; font-size: 1.1em; }}
  .container {{ max-width: 1200px; margin: 30px auto; padding: 0 20px; }}
  .card {{ background: white; border-radius: 8px; padding: 24px;
           margin-bottom: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
  h2 {{ color: #1A237E; border-bottom: 2px solid #E8EAF6;
        padding-bottom: 8px; margin-top: 0; }}
  h3 {{ color: #37474F; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(4, 1fr);
                   gap: 16px; margin-bottom: 8px; }}
  .stat-box {{ background: #F5F7FF; border-radius: 6px; padding: 16px;
               text-align: center; border-left: 4px solid #3949AB; }}
  .stat-box.pass {{ border-left-color: #4CAF50; }}
  .stat-box.warn {{ border-left-color: #FF9800; }}
  .stat-box.fail {{ border-left-color: #F44336; }}
  .stat-box .val {{ font-size: 2em; font-weight: bold; color: #1A237E; }}
  .stat-box.pass .val {{ color: #2E7D32; }}
  .stat-box.warn .val {{ color: #E65100; }}
  .stat-box.fail .val {{ color: #B71C1C; }}
  .stat-box .lbl {{ font-size: 0.85em; color: #666; margin-top: 4px; }}
  .plots-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  img {{ max-width: 100%; border-radius: 4px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; }}
  th {{ background: #1A237E; color: white; padding: 8px 10px;
        text-align: left; position: sticky; top: 0; }}
  td {{ border-bottom: 1px solid #eee; padding: 6px 10px; }}
  .pass td {{ background: #F1F8E9; }}
  .warn td {{ background: #FFF8E1; }}
  .fail td {{ background: #FFEBEE; }}
  .table-wrap {{ max-height: 500px; overflow-y: auto; border-radius: 6px;
                 border: 1px solid #ddd; }}
  pre {{ background: #F5F5F5; padding: 16px; border-radius: 6px;
         font-size: 0.8em; overflow-x: auto; white-space: pre-wrap; }}
  .methods p {{ line-height: 1.7; text-align: justify; }}
  .footer {{ text-align: center; color: #999; font-size: 0.8em; padding: 20px; }}
  .muted {{ color: #999; font-style: italic; }}
  .nav {{ background: white; border-bottom: 1px solid #ddd; padding: 12px 20px;
          position: sticky; top: 0; z-index: 100;
          box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
  .nav a {{ color: #3949AB; text-decoration: none; margin-right: 20px;
            font-size: 0.9em; font-weight: 500; }}
  .nav a:hover {{ text-decoration: underline; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 10px;
            font-size: 0.8em; font-weight: bold; }}
  .badge-blue  {{ background:#E3F2FD; color:#1565C0; }}
  .badge-green {{ background:#E8F5E9; color:#2E7D32; }}
</style>
</head>
<body>

<div class="header">
  <h1>Mitogenome Extraction Pipeline{report_subtitle}</h1>
  <p><i>{SPECIES}</i> &nbsp;|&nbsp; {datetime.now().strftime("%d %B %Y, %H:%M")}
     &nbsp;|&nbsp; {len(POP_ORDER)} populations
     {"&nbsp;|&nbsp; min pop size: " + str(MIN_POP_SIZE) if IS_FILTERED else ""}
  </p>
</div>

<div class="nav">
  <a href="#summary">Summary</a>
  <a href="#coverage">Coverage</a>
  <a href="#qc-table">Per-Sample QC</a>
  <a href="#consensus">Consensus</a>
  <a href="#annotation">Annotation</a>
  <a href="#phylo">Phylogenetics</a>
  <a href="#dnds">dN/dS</a>
  <a href="#bsp">Bayesian Skyline</a>
  <a href="#methods">Methods</a>
  <a href="#outputs">Output Files</a>
</div>

<div class="container">

{excl_banner}

<!-- ── Run Summary ──────────────────────────────────────────────────────── -->
<div class="card" id="summary">
  <h2>Run Summary</h2>
  <div class="summary-grid">
    <div class="stat-box">
      <div class="val">{n_total}</div><div class="lbl">Total Samples</div></div>
    <div class="stat-box pass">
      <div class="val">{n_pass}</div>
      <div class="lbl">PASS ({n_pass/n_total*100:.0f}%)</div></div>
    <div class="stat-box warn">
      <div class="val">{n_warn}</div>
      <div class="lbl">WARN ({n_warn/n_total*100:.0f}%)</div></div>
    <div class="stat-box fail">
      <div class="val">{n_fail}</div>
      <div class="lbl">FAIL ({n_fail/n_total*100:.0f}%)</div></div>
  </div>
  <div class="summary-grid" style="margin-top:12px">
    <div class="stat-box">
      <div class="val">{median(depths):.0f}x</div>
      <div class="lbl">Median Depth</div></div>
    <div class="stat-box">
      <div class="val">{median(breadths):.1f}%</div>
      <div class="lbl">Median Breadth</div></div>
    <div class="stat-box">
      <div class="val">{cons_len:,}</div>
      <div class="lbl">Consensus Length (bp)</div></div>
    <div class="stat-box">
      <div class="val">{gc:.1f}%</div>
      <div class="lbl">GC Content</div></div>
  </div>
  <div class="summary-grid" style="margin-top:12px">
    <div class="stat-box">
      <div class="val">{annot_stats.get('Protein coding', '—')}</div>
      <div class="lbl">Protein-Coding Genes</div></div>
    <div class="stat-box">
      <div class="val">{annot_stats.get('rRNA', '—')}</div>
      <div class="lbl">rRNA Genes</div></div>
    <div class="stat-box">
      <div class="val">{annot_stats.get('tRNA', '—')}</div>
      <div class="lbl">tRNA Genes</div></div>
    <div class="stat-box">
      <div class="val">{len(yn00_rows)}</div>
      <div class="lbl">Genes with dN/dS</div></div>
  </div>
</div>

<!-- ── Coverage ─────────────────────────────────────────────────────────── -->
<div class="card" id="coverage">
  <h2>Coverage Statistics</h2>
  <div class="plots-grid">
    <div><h3>Depth Distribution</h3>
      <img src="data:image/png;base64,{depth_hist_b64}" /></div>
    <div><h3>Coverage Breadth</h3>
      <img src="data:image/png;base64,{breadth_hist_b64}" /></div>
    <div><h3>Mapped Reads vs Depth</h3>
      <img src="data:image/png;base64,{scatter_b64}" /></div>
    <div><h3>QC Status</h3>
      <img src="data:image/png;base64,{pie_b64}" /></div>
  </div>
</div>

<!-- ── Per-Sample QC ────────────────────────────────────────────────────── -->
<div class="card" id="qc-table">
  <h2>Per-Sample QC</h2>
  <div class="table-wrap">
  <table>
    <thead><tr>
      <th>Sample</th><th>Mapped Reads</th><th>Mean Depth</th>
      <th>Breadth</th><th>Variants</th><th>Cons. Length</th>
      <th>N Count</th><th>Status</th>
    </tr></thead>
    <tbody>{qc_rows}</tbody>
  </table>
  </div>
</div>

<!-- ── Consensus ────────────────────────────────────────────────────────── -->
<div class="card" id="consensus">
  <h2>Consensus Mitogenome</h2>
  <div class="plots-grid">
    <div><h3>Base Composition</h3>
      <img src="data:image/png;base64,{basecomp_b64}" /></div>
    <div><h3>Sequence Statistics</h3>
      <table>
        <tr><th>Property</th><th>Value</th></tr>
        <tr><td>Length</td><td>{cons_len:,} bp</td></tr>
        <tr><td>GC content</td><td>{gc:.2f}%</td></tr>
        <tr><td>AT content</td><td>{100-gc:.2f}%</td></tr>
        <tr><td>N count</td><td>{base_counts.get('N',0):,}</td></tr>
        <tr><td>A count</td><td>{base_counts.get('A',0):,}</td></tr>
        <tr><td>T count</td><td>{base_counts.get('T',0):,}</td></tr>
        <tr><td>G count</td><td>{base_counts.get('G',0):,}</td></tr>
        <tr><td>C count</td><td>{base_counts.get('C',0):,}</td></tr>
        {'<tr><td>Closest relative</td><td><i>' + annot_related + '</i></td></tr>'
         if annot_related else ''}
      </table>
    </div>
  </div>
  <h3>First 600 bp</h3>
  <pre>{consensus[:600]}</pre>
</div>

<!-- ── Annotation ───────────────────────────────────────────────────────── -->
<div class="card" id="annotation">
  <h2>Mitogenome Annotation</h2>
  <p>Annotation was performed using MitoZ v3.5 with the {CLADE} CDS database
     and {gc_desc} genetic code (translation table {GENETIC_CODE}).
     {f'Closest reference species: <i>{annot_related}</i>.' if annot_related else ''}
  </p>
  <div class="summary-grid" style="margin-bottom:16px">
    <div class="stat-box">
      <div class="val">{annot_stats.get('Protein coding','—')}</div>
      <div class="lbl">Protein-Coding Genes <span class="badge badge-blue">CDS</span></div></div>
    <div class="stat-box">
      <div class="val">{annot_stats.get('rRNA','—')}</div>
      <div class="lbl">rRNA Genes <span class="badge badge-green">rRNA</span></div></div>
    <div class="stat-box">
      <div class="val">{annot_stats.get('tRNA','—')}</div>
      <div class="lbl">tRNA Genes</div></div>
    <div class="stat-box">
      <div class="val">{sum(annot_stats.values()) if annot_stats else '—'}</div>
      <div class="lbl">Total Genes Found</div></div>
  </div>
  {f'''<div class="table-wrap"><table>
    <thead><tr>
      <th>Gene</th><th>Product</th><th>Type</th>
      <th>Position</th><th>Length</th><th>Strand</th>
    </tr></thead>
    <tbody>{annot_rows}</tbody>
  </table></div>''' if annot_rows else '<p class="muted">Annotation results not found.</p>'}
</div>

<!-- ── Phylogenetics ─────────────────────────────────────────────────────── -->
<div class="card" id="phylo">
  <h2>Population-Level Phylogenetic Tree</h2>
  <p>Neighbour-joining tree built from majority-rule consensus sequences per population,
     using identity distance. Tip labels show population code and sample count.
     {"Only populations with n &ge; " + str(MIN_POP_SIZE) + " are shown." if IS_FILTERED else ""}
  </p>
  {"<img src='data:image/png;base64," + tree_b64 + "' />" if tree_b64
   else "<p class='muted'>Tree could not be rendered.</p>"}
  <h2 style="margin-top:24px">Population Network</h2>
  <p>Each node represents a population. Node size is proportional to sample count.
     Edge thickness reflects genetic similarity (thicker = more similar).</p>
  {"<img src='data:image/png;base64," + haplo_b64 + "' />" if haplo_b64
   else "<p class='muted'>Network could not be generated.</p>"}
</div>

<!-- ── dN/dS ────────────────────────────────────────────────────────────── -->
<div class="card" id="dnds">
  <h2>Selection Analysis — dN/dS</h2>
  <p>Synonymous (dS) and non-synonymous (dN) substitution rates were estimated
     for each of the {len(annot_genes)} protein-coding mitochondrial genes.
     Pairwise dN/dS ratios were calculated using yn00 (Yang &amp; Nielsen 2000, PAML)
     across all samples. A branch model test was additionally performed using codeml
     (model=2) to identify population-specific deviations in ω (dN/dS),
     with likelihood ratio tests against a null model (model=0).
     {"Branch model analyses were restricted to populations with n &ge; " + str(MIN_POP_SIZE) + "." if IS_FILTERED else ""}
     Values &lt; 1 indicate purifying selection; values &gt; 1 indicate positive selection.</p>

  {excl_banner if excl_display else ""}

  <h3>Pairwise dN/dS (yn00)</h3>
  {f'<img src="data:image/png;base64,{dnds_fig_b64}" />' if dnds_fig_b64
   else '<p class="muted">yn00 results not yet available.</p>'}
  {f'''<div class="table-wrap" style="margin-top:16px"><table>
    <thead><tr>
      <th>Gene</th><th>Pairs</th><th>Mean dN</th><th>Mean dS</th>
      <th>Mean dN/dS</th><th>Min dN/dS</th><th>Max dN/dS</th>
    </tr></thead>
    <tbody>{yn00_table_rows}</tbody>
  </table></div>''' if yn00_table_rows else ''}

  <h3 style="margin-top:24px">codeml Branch Model — Population-Specific Selection</h3>
  <p>Foreground ω per population relative to background. Significant (LRT p&lt;0.05) marked ✱.</p>
  {f'<img src="data:image/png;base64,{codeml_fig_b64}" />' if codeml_fig_b64
   else '<p class="muted">codeml results not yet available.</p>'}
  {f'''<div class="table-wrap" style="margin-top:16px"><table>
    <thead><tr>
      <th>Gene</th><th>Population</th><th>lnL (null)</th><th>lnL (branch)</th>
      <th>LRT</th><th>p&lt;0.05</th><th>ω background</th><th>ω foreground</th>
    </tr></thead>
    <tbody>{codeml_table_rows}</tbody>
  </table></div>''' if codeml_table_rows else ''}

  {f'''<h3 style="margin-top:24px">Excluded Populations</h3>
  <div class="table-wrap"><table>
    <thead><tr><th>Population</th><th>n samples</th><th>Reason</th></tr></thead>
    <tbody>{excl_table_rows}</tbody>
  </table></div>''' if excl_table_rows else ''}
</div>

<!-- ── Bayesian Skyline ──────────────────────────────────────────────────── -->
<div class="card" id="bsp">
  <h2>Bayesian Skyline Plot — Demographic History</h2>
  <p>Bayesian Skyline Plots were estimated from the CYTB alignment using BEAST2
     v2.6.3 under an uncorrelated lognormal relaxed clock. The clock rate was given
     a uniform prior of {clock_desc}.
     Substitution was modelled as HKY+Γ<sub>4</sub>. Runs were performed
     separately for each {"included " if IS_FILTERED else ""}population
     {"(n &ge; " + str(MIN_POP_SIZE) + ") " if IS_FILTERED else ""}
     and for all samples combined (50 million generations, sampled every 5,000).
     Convergence was assessed by ESS of the posterior (threshold: 200).</p>

  <h3>BEAST2 Convergence (ESS)</h3>
  {f'<img src="data:image/png;base64,{bsp_ess_fig_b64}" />' if bsp_ess_fig_b64
   else '<p class="muted">BSP ESS results not yet available.</p>'}
  {f'''<div class="table-wrap" style="margin-top:16px"><table>
    <thead><tr>
      <th>Run</th><th>ESS (posterior)</th><th>Status</th>
    </tr></thead>
    <tbody>{bsp_ess_table_rows}</tbody>
  </table></div>
  <p style="margin-top:8px;font-size:0.85em;color:#666">
    ⚠ Runs with ESS &lt; 200 should be re-run with a longer chain (100M generations).
    Use the R script at <code>{OUTDIR}/bsp/plot_bsp.R</code> to generate
    full skyline plots with 95% HPD credible intervals.
  </p>''' if bsp_ess_table_rows else ''}
</div>

<!-- ── Methods ──────────────────────────────────────────────────────────── -->
<div class="card methods" id="methods">
  <h2>Methods</h2>
  <p><b>Read mapping:</b> Paired-end reads were mapped to the reference mitogenome
     (<i>{SPECIES}</i>) using BWA-MEM v0.7.19. Multi-lane samples were concatenated
     prior to mapping. Mapped reads were coordinate-sorted using SAMtools.</p>
  <p><b>Duplicate removal:</b> PCR duplicates were identified and removed using SAMtools
     markdup, following queryname sorting and mate-score annotation with SAMtools fixmate.
     Reads were filtered to retain primary, mapped alignments with mapping quality &ge; 20.</p>
  <p><b>Consensus generation:</b> Variants were called using BCFtools mpileup and BCFtools
     call (multiallelic model). Consensus sequences were generated requiring minimum depth
     of 3x and minimum allele frequency of 0.70.</p>
  <p><b>Alignment:</b> All consensus sequences were aligned using MAFFT v7.525 (--auto).
     A majority-rule species consensus was derived from the alignment.</p>
  <p><b>Phylogenetics:</b> A maximum-likelihood tree was inferred using IQ-TREE v3.0.1
     under the GTR+G model with 1000 ultrafast bootstrap replicates. A population-level
     NJ tree was additionally constructed from per-population consensus sequences.</p>
  <p><b>Mitogenome annotation:</b> The consensus mitogenome was annotated using MitoZ
     v3.5 with the {CLADE} CDS protein database and {gc_desc} genetic code
     (translation table {GENETIC_CODE}). Protein-coding genes were identified using
     GeneWise; ribosomal RNAs were detected using Infernal with MitoZ covariance models.</p>
  <p><b>Selection analysis:</b> Per-gene CDS sequences were extracted from all
     {n_total} individual mitogenomes using coordinates derived from the consensus
     annotation, aligned with MAFFT, and converted to PHYLIP format for PAML v4.
     Pairwise dN/dS ratios were estimated using yn00 (Yang &amp; Nielsen 2000).
     A branch model test (codeml model=2 vs model=0) was performed for each
     gene × population combination for populations with n &ge; {MIN_POP_SIZE if IS_FILTERED else 1};
     significance was assessed by likelihood ratio test (&chi;<sup>2</sup>, 1 d.f.,
     p&lt;0.05 threshold LRT&gt;3.841).</p>
  <p><b>Demographic history:</b> Bayesian Skyline Plots were estimated for each
     {"included " if IS_FILTERED else ""}population
     {"(n &ge; " + str(MIN_POP_SIZE) + ") " if IS_FILTERED else ""}
     and for all samples combined from the CYTB alignment using BEAST2 v2.6.3.
     Analyses assumed an HKY+&Gamma;<sub>4</sub> substitution model and an
     uncorrelated lognormal relaxed clock with a uniform prior on the mean clock
     rate ({clock_desc}). Each run comprised 50 million MCMC generations sampled
     every 5,000 steps; the first 10% were discarded as burnin. Convergence was
     assessed by ESS &ge; 200 for the posterior distribution.</p>
</div>

<!-- ── Output Files ─────────────────────────────────────────────────────── -->
<div class="card" id="outputs">
  <h2>Output Files</h2>
  <table>
    <tr><th>File</th><th>Description</th></tr>
    <tr><td><code>mito_report_all.html</code></td>
        <td>Full report including all populations</td></tr>
    <tr><td><code>mito_report_filtered.html</code></td>
        <td>Report filtered to populations with n &ge; {MIN_POP_SIZE if IS_FILTERED else "N"}</td></tr>
    <tr><td><code>aligned/all_samples_mito_aligned.fa</code></td>
        <td>MAFFT multiple sequence alignment</td></tr>
    <tr><td><code>aligned/consensus_mito.fa</code></td>
        <td>Majority-rule species consensus</td></tr>
    <tr><td><code>tree/mito_tree.treefile</code></td>
        <td>IQ-TREE maximum-likelihood tree (Newick)</td></tr>
    <tr><td><code>tree/mito_tree.contree</code></td>
        <td>Consensus tree with bootstrap values</td></tr>
    <tr><td><code>qc/pipeline_summary.tsv</code></td>
        <td>Per-sample QC metrics</td></tr>
    <tr><td><code>annotation_test/.../summary.txt</code></td>
        <td>MitoZ annotation summary</td></tr>
    <tr><td><code>annotation_test/.../*.gbf</code></td>
        <td>GenBank format annotation</td></tr>
    <tr><td><code>annotation_test/.../circos.png</code></td>
        <td>Circular mitogenome map</td></tr>
    <tr><td><code>dnds/alignments/</code></td>
        <td>Per-gene MAFFT alignments (FASTA + PHYLIP)</td></tr>
    <tr><td><code>dnds/yn00_summary.tsv</code></td>
        <td>Pairwise dN/dS per gene (yn00)</td></tr>
    <tr><td><code>dnds/codeml_branch_summary.tsv</code></td>
        <td>Branch model LRT results (included populations)</td></tr>
    <tr><td><code>dnds/codeml_excluded_pops.tsv</code></td>
        <td>Populations excluded from codeml branch model</td></tr>
    <tr><td><code>dnds/codeml/pop_sizes.json</code></td>
        <td>Population sample counts and filter metadata</td></tr>
    <tr><td><code>bsp/BSP_*.xml</code></td>
        <td>BEAST2 XML input files</td></tr>
    <tr><td><code>bsp/BSP_*/BSP_*.log</code></td>
        <td>BEAST2 MCMC trace logs</td></tr>
    <tr><td><code>bsp/BSP_*/BSP_*_MCC.tree</code></td>
        <td>Maximum clade credibility trees</td></tr>
    <tr><td><code>bsp/ess_summary.tsv</code></td>
        <td>ESS convergence check per BSP run</td></tr>
    <tr><td><code>bsp/bsp_pop_summary.json</code></td>
        <td>BSP population filter metadata</td></tr>
    <tr><td><code>bsp/plot_bsp.R</code></td>
        <td>R script to generate final BSP figures</td></tr>
    <tr><td><code>bam/*_final.bam</code></td>
        <td>Final filtered BAM files</td></tr>
  </table>
</div>

</div>
<div class="footer">
  Generated by mito_pipeline &nbsp;|&nbsp;
  {datetime.now().strftime("%Y-%m-%d %H:%M")} &nbsp;|&nbsp;
  {"Filtered report (n&ge;" + str(MIN_POP_SIZE) + ")" if IS_FILTERED else "Full report (all populations)"}
</div>
</body>
</html>"""

with open(OUTHTML, "w", encoding="utf-8") as f:
    f.write(html)

print(f"Report written: {OUTHTML}")
if excl_display:
    print(f"Excluded populations: {sorted(excl_display.keys())}")
",
      "path": "/mnt/user-data/outputs/mito_report.py"
    },
    "message": "Updated mito_report.py with min-pop-size, dynamic species/clade/clock params, and dual report support",
    "integration_name": null,
    "integration_icon_url": null,
    "icon_name": "file",
    "context": null,
    "display_content": {
      "type": "json_block",
      "json_block": "{"language": "python", "code": "#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime
from io import BytesIO
import base64

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import networkx as nx

from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate mitogenome HTML report")
parser.add_argument("--outdir",       required=True)
parser.add_argument("--species",      default="Unknown species")
parser.add_argument("--output",       default="mito_report.html")
parser.add_argument("--min-pop-size", type=int, default=0,
                    help="Minimum samples per population for per-pop analyses. "
                         "0 = include all (default). Excluded pops shown as greyed "
                         "out in tables and omitted from heatmap/BSP figures.")
parser.add_argument("--genetic-code", type=int, default=2,
                    help="NCBI genetic code used (for methods text)")
parser.add_argument("--clade",        default="Chordata",
                    help="MitoZ clade used (for methods text)")
parser.add_argument("--clock-rate-min", type=float, default=0.005,
                    help="BEAST2 clock rate prior lower bound (for methods text)")
parser.add_argument("--clock-rate-max", type=float, default=0.02,
                    help="BEAST2 clock rate prior upper bound (for methods text)")
args = parser.parse_args()

OUTDIR        = args.outdir
SPECIES       = args.species
OUTHTML       = args.output
MIN_POP_SIZE  = args.min_pop_size
GENETIC_CODE  = args.genetic_code
CLADE         = args.clade
CLOCK_MIN     = args.clock_rate_min
CLOCK_MAX     = args.clock_rate_max

# Determine if this is a filtered report
IS_FILTERED   = MIN_POP_SIZE > 0

# ── File paths ─────────────────────────────────────────────────────────────────
SUMMARY_TSV        = os.path.join(OUTDIR, "qc",     "pipeline_summary.tsv")
CONSENSUS_FA       = os.path.join(OUTDIR, "aligned", "consensus_mito.fa")
ALL_MITO_FA        = os.path.join(OUTDIR, "aligned", "all_samples_mito_aligned.fa")
TREEFILE           = os.path.join(OUTDIR, "tree",    "mito_tree.treefile")
YN00_TSV           = os.path.join(OUTDIR, "dnds",   "yn00_summary.tsv")
CODEML_TSV         = os.path.join(OUTDIR, "dnds",   "codeml_branch_summary.tsv")
CODEML_EXCL_TSV    = os.path.join(OUTDIR, "dnds",   "codeml_excluded_pops.tsv")
BSP_ESS_TSV        = os.path.join(OUTDIR, "bsp",    "ess_summary.tsv")
BSP_POP_JSON       = os.path.join(OUTDIR, "bsp",    "bsp_pop_summary.json")
CODEML_POP_JSON    = os.path.join(OUTDIR, "dnds",   "codeml", "pop_sizes.json")
ANNOTATION_SUMMARY = os.path.join(OUTDIR, "annotation_test",
                     "test.consensus_mito.fa.result", "summary.txt")

# ── Helpers ───────────────────────────────────────────────────────────────────
def median(lst):
    if not lst: return 0
    lst = sorted(lst)
    n = len(lst)
    return lst[n//2] if n % 2 else (lst[n//2-1] + lst[n//2]) / 2

def safe_float(x):
    try: return float(x)
    except: return None

def read_fasta_seq(path):
    seq = ""
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                if not line.startswith(">"): seq += line.strip()
    return seq

def encode_fig(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def get_pop(name):
    m = re.search(r'-([A-Z]+)\d', name)
    return m.group(1) if m else "Unknown"

def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def load_json(path):
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)

# ── Load pop size metadata ─────────────────────────────────────────────────────
bsp_pop_summary    = load_json(BSP_POP_JSON)
codeml_pop_summary = load_json(CODEML_POP_JSON)

# Determine which pops were excluded by each analysis
bsp_excluded    = set(bsp_pop_summary.get('excluded_pops', {}).keys())
codeml_excluded = set(codeml_pop_summary.get('excluded_pops', {}).keys())
all_pop_sizes   = bsp_pop_summary.get('all_pops', codeml_pop_summary.get('all_pops', {}))

# For the filtered report, also apply the --min-pop-size filter at report level
# (in case report is run with a different threshold than the pipeline)
if IS_FILTERED:
    report_excluded = {p for p, n in all_pop_sizes.items() if n < MIN_POP_SIZE}
else:
    report_excluded = set()

# ── Load QC summary ───────────────────────────────────────────────────────────
samples = []
with open(SUMMARY_TSV) as f:
    header = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2 or parts[0] in ("", "0"): continue
        samples.append(dict(zip(header, parts)))

n_total = len(samples)
n_pass  = sum(1 for s in samples if "PASS" in s.get("status",""))
n_warn  = sum(1 for s in samples if "WARN" in s.get("status",""))
n_fail  = sum(1 for s in samples if "FAIL" in s.get("status",""))

depths, breadths, mapped_reads = [], [], []
for s in samples:
    d = safe_float(s.get("mean_depth"))
    b = safe_float(s.get("coverage_breadth_pct"))
    m = safe_float(s.get("mapped_reads"))
    if d is not None: depths.append(d)
    if b is not None: breadths.append(b)
    if m is not None: mapped_reads.append(m)

# ── Consensus ─────────────────────────────────────────────────────────────────
consensus   = read_fasta_seq(CONSENSUS_FA)
cons_len    = len(consensus)
base_counts = Counter(consensus.upper())
gc = (base_counts.get("G",0) + base_counts.get("C",0)) / cons_len * 100 if cons_len else 0

# ── Load alignment + population codes ─────────────────────────────────────────
print("Loading alignment...")
records = list(SeqIO.parse(ALL_MITO_FA, "fasta"))

seqs_by_pop = defaultdict(list)
for rec in records:
    pop = get_pop(rec.id)
    seqs_by_pop[pop].append(str(rec.seq))

# All pops and filtered pops
ALL_POP_ORDER      = sorted(seqs_by_pop.keys())
FILTERED_POP_ORDER = sorted(p for p in ALL_POP_ORDER if p not in report_excluded)
POP_ORDER          = FILTERED_POP_ORDER if IS_FILTERED else ALL_POP_ORDER
N_POPS             = len(POP_ORDER)

print(f"  All populations:      {ALL_POP_ORDER}")
if IS_FILTERED:
    print(f"  Filtered populations: {FILTERED_POP_ORDER}")
    print(f"  Excluded:             {sorted(report_excluded)}")

cmap = matplotlib.colormaps.get_cmap("tab10").resampled(max(len(ALL_POP_ORDER), 1))
pop_colors_map = {pop: cmap(i) for i, pop in enumerate(ALL_POP_ORDER)}

# ── Load annotation summary ───────────────────────────────────────────────────
print("Loading annotation data...")
annot_genes   = []
annot_stats   = {}
annot_related = ""
if os.path.exists(ANNOTATION_SUMMARY):
    with open(ANNOTATION_SUMMARY) as f:
        content = f.read()
    for line in content.split('\n'):
        parts = line.split()
        if len(parts) >= 9 and re.match(r'\S+_mito', parts[0]):
            annot_genes.append({
                'seq_id':    parts[0],
                'start':     parts[1],
                'end':       parts[2],
                'length_bp': parts[3],
                'direction': parts[4],
                'type':      parts[5],
                'gene':      parts[6],
                'product':   ' '.join(parts[7:-1]),
                'freq':      parts[-1],
            })
    for line in content.split('\n'):
        m = re.match(r'(Protein coding|tRNA|rRNA) genes totally found:\s+(\d+)', line)
        if m:
            annot_stats[m.group(1)] = int(m.group(2))
    m2 = re.search(r'Closely_related_species\s*\n\S+\s+\d+\s+\S+\s+(.+)', content)
    if m2:
        annot_related = m2.group(1).strip()

# ── Load dN/dS results ────────────────────────────────────────────────────────
print("Loading dN/dS results...")
yn00_rows   = read_tsv(YN00_TSV)
codeml_rows_all = read_tsv(CODEML_TSV)

# Filter codeml rows for filtered report
if IS_FILTERED:
    codeml_rows = [r for r in codeml_rows_all
                   if r.get('population','') not in report_excluded]
else:
    codeml_rows = codeml_rows_all

# Collect excluded pop info for display
excl_rows = read_tsv(CODEML_EXCL_TSV)
# Merge pipeline exclusions with report-level exclusions
excl_display = {}
for r in excl_rows:
    excl_display[r['population']] = r['n_samples']
if IS_FILTERED:
    for p in report_excluded:
        if p not in excl_display:
            excl_display[p] = all_pop_sizes.get(p, '?')

# ── Load BSP ESS summary ──────────────────────────────────────────────────────
print("Loading BSP ESS data...")
bsp_ess_rows_all = read_tsv(BSP_ESS_TSV)

if IS_FILTERED:
    bsp_ess_rows = [r for r in bsp_ess_rows_all
                    if r.get('run','').replace('BSP_','') not in report_excluded
                    or r.get('run','') == 'BSP_all']
else:
    bsp_ess_rows = bsp_ess_rows_all

# ── Figure 1: Depth histogram ─────────────────────────────────────────────────
print("Generating figures...")
fig1, ax1 = plt.subplots(figsize=(7, 4))
ax1.hist(depths, bins=30, color="#2196F3", edgecolor="white", linewidth=0.5)
ax1.axvline(median(depths), color="#F44336", linestyle="--", linewidth=1.5,
            label=f"Median: {median(depths):.0f}x")
ax1.set_xlabel("Mean depth (x)", fontsize=11)
ax1.set_ylabel("Number of samples", fontsize=11)
ax1.set_title("Coverage Depth Distribution", fontsize=12, fontweight="bold")
ax1.legend()
ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)
depth_hist_b64 = encode_fig(fig1)

# ── Figure 2: Breadth histogram ───────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.hist(breadths, bins=30, color="#4CAF50", edgecolor="white", linewidth=0.5)
ax2.axvline(median(breadths), color="#F44336", linestyle="--", linewidth=1.5,
            label=f"Median: {median(breadths):.1f}%")
ax2.set_xlabel("Coverage breadth (%)", fontsize=11)
ax2.set_ylabel("Number of samples", fontsize=11)
ax2.set_title("Genome Coverage Breadth", fontsize=12, fontweight="bold")
ax2.legend()
ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
breadth_hist_b64 = encode_fig(fig2)

# ── Figure 3: Scatter mapped reads vs depth ───────────────────────────────────
fig3, ax3 = plt.subplots(figsize=(7, 4))
status_colors = []
for s in samples:
    st = s.get("status","")
    if "PASS" in st:   status_colors.append("#4CAF50")
    elif "WARN" in st: status_colors.append("#FF9800")
    else:              status_colors.append("#F44336")
ax3.scatter(mapped_reads, depths, c=status_colors, alpha=0.6, s=30, edgecolors="none")
ax3.set_xlabel("Mapped reads", fontsize=11)
ax3.set_ylabel("Mean depth (x)", fontsize=11)
ax3.set_title("Mapped Reads vs Coverage Depth", fontsize=12, fontweight="bold")
ax3.legend(handles=[
    mpatches.Patch(color="#4CAF50", label="PASS"),
    mpatches.Patch(color="#FF9800", label="WARN"),
    mpatches.Patch(color="#F44336", label="FAIL"),
])
ax3.spines["top"].set_visible(False); ax3.spines["right"].set_visible(False)
scatter_b64 = encode_fig(fig3)

# ── Figure 4: Status pie ──────────────────────────────────────────────────────
fig4, ax4 = plt.subplots(figsize=(5, 4))
pie_vals   = [x for x in [n_pass, n_warn, n_fail] if x > 0]
pie_labels = [f"{l}\n({v})" for l, v in zip(["PASS","WARN","FAIL"],
              [n_pass, n_warn, n_fail]) if v > 0]
pie_colors = [c for c, v in zip(["#4CAF50","#FF9800","#F44336"],
              [n_pass, n_warn, n_fail]) if v > 0]
ax4.pie(pie_vals, labels=pie_labels, colors=pie_colors,
        startangle=90, textprops={"fontsize": 10})
ax4.set_title("QC Status", fontsize=12, fontweight="bold")
pie_b64 = encode_fig(fig4)

# ── Figure 5: Base composition ────────────────────────────────────────────────
fig5, ax5 = plt.subplots(figsize=(5, 4))
bases     = ["A","T","G","C","N"]
base_pcts = [base_counts.get(b,0)/cons_len*100 for b in bases]
bar_cols  = ["#2196F3","#F44336","#4CAF50","#FF9800","#9E9E9E"]
bars = ax5.bar(bases, base_pcts, color=bar_cols, edgecolor="white")
for bar, pct in zip(bars, base_pcts):
    if pct > 0.5:
        ax5.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.3,
                 f"{pct:.1f}%", ha="center", va="bottom", fontsize=9)
ax5.set_xlabel("Base"); ax5.set_ylabel("Percentage (%)")
ax5.set_title("Consensus Base Composition", fontsize=12, fontweight="bold")
ax5.spines["top"].set_visible(False); ax5.spines["right"].set_visible(False)
basecomp_b64 = encode_fig(fig5)

# ── Figure 6: Population NJ tree ─────────────────────────────────────────────
print("Building population NJ tree...")
tree_b64 = None
try:
    pop_cons_records = []
    for pop in POP_ORDER:
        seqs   = seqs_by_pop[pop]
        length = len(seqs[0])
        cons_seq = []
        for i in range(length):
            col = [s[i] for s in seqs if i < len(s) and s[i] not in "-N"]
            cons_seq.append(Counter(col).most_common(1)[0][0] if col else "N")
        pop_cons_records.append(SeqRecord(Seq("".join(cons_seq)), id=pop, description=""))

    tmp_fa = os.path.join(OUTDIR, "tmp_pop_consensus.fa")
    SeqIO.write(pop_cons_records, tmp_fa, "fasta")
    alignment   = AlignIO.read(tmp_fa, "fasta")
    calculator  = DistanceCalculator("identity")
    dm          = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    pop_tree    = constructor.nj(dm)
    os.remove(tmp_fa)

    fig6, ax6 = plt.subplots(figsize=(10, 6))
    Phylo.draw(pop_tree, do_show=False, axes=ax6,
               label_func=lambda x: x.name if x.is_terminal() else "",
               label_colors={pop: pop_colors_map[pop] for pop in POP_ORDER})
    for text_obj in ax6.texts:
        label = text_obj.get_text().strip()
        if label in seqs_by_pop:
            text_obj.set_text(f"{label}  (n={len(seqs_by_pop[label])})")
    ax6.set_title(f"Population-Level Mitogenome NJ Tree\n{SPECIES}",
                  fontsize=13, fontweight="bold")
    ax6.spines["top"].set_visible(False); ax6.spines["right"].set_visible(False)
    ax6.legend(handles=[
        mpatches.Patch(color=pop_colors_map[p],
                       label=f"{p}  (n={len(seqs_by_pop[p])})")
        for p in POP_ORDER
    ], title="Population", loc="lower right", fontsize=9)
    plt.subplots_adjust(left=0.05, right=0.75, top=0.92, bottom=0.05)
    tree_b64 = encode_fig(fig6)
    print("  Population tree done")
except Exception as e:
    print(f"  Tree failed: {e}")

# ── Figure 7: Population network ──────────────────────────────────────────────
print("Building population network...")
haplo_b64 = None
try:
    seq_length = len(records[0].seq)
    dist_between_pops = {}
    for i, p1 in enumerate(POP_ORDER):
        for j, p2 in enumerate(POP_ORDER):
            if i >= j: continue
            dists = []
            for s1 in seqs_by_pop[p1]:
                for s2 in seqs_by_pop[p2]:
                    d = sum(a != b and a not in "-N" and b not in "-N"
                            for a, b in zip(s1, s2)) / seq_length
                    dists.append(d)
            dist_between_pops[(p1, p2)] = np.mean(dists) if dists else 1.0

    G = nx.Graph()
    for pop in POP_ORDER:
        G.add_node(pop, size=len(seqs_by_pop[pop]))
    for (p1, p2), dist in dist_between_pops.items():
        G.add_edge(p1, p2, weight=dist)

    pos      = nx.spring_layout(G, weight="weight", seed=42, k=3)
    max_dist = max(dist_between_pops.values()) if dist_between_pops else 1
    edge_widths = [max(0.5, (1 - G[u][v]["weight"] / max_dist) * 6)
                   for u, v in G.edges()]

    fig7, ax7 = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, ax=ax7, width=edge_widths,
                           alpha=0.5, edge_color="#888888")
    nx.draw_networkx_nodes(G, pos, ax=ax7,
                           node_color=[pop_colors_map[p] for p in G.nodes()],
                           node_size=[len(seqs_by_pop[p])*60+200 for p in G.nodes()],
                           alpha=0.9)
    nx.draw_networkx_labels(G, pos, ax=ax7,
                            labels={p: f"{p}\n(n={len(seqs_by_pop[p])})"
                                    for p in G.nodes()},
                            font_size=9, font_weight="bold")
    ax7.set_title(f"Population Network — {SPECIES}\n"
                  "Node size = sample count, edge thickness = genetic similarity",
                  fontsize=12, fontweight="bold")
    ax7.axis("off")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.02)
    haplo_b64 = encode_fig(fig7)
    print("  Population network done")
except Exception as e:
    print(f"  Network failed: {e}")

# ── Figure 8: dN/dS per gene bar chart ───────────────────────────────────────
print("Generating dN/dS figure...")
dnds_fig_b64 = None
if yn00_rows:
    try:
        genes     = [r['gene'] for r in yn00_rows]
        mean_dnds = [float(r['mean_dNdS']) for r in yn00_rows]
        colors    = ['#F44336' if v > 1 else '#2196F3' for v in mean_dnds]

        fig8, ax8 = plt.subplots(figsize=(12, 5))
        bars = ax8.bar(genes, mean_dnds, color=colors, edgecolor='white')
        ax8.axhline(1.0, color='black', linestyle='--', linewidth=1,
                    label='dN/dS = 1 (neutral)')
        ax8.axhline(0.0, color='grey', linestyle='-', linewidth=0.5)
        for bar, val in zip(bars, mean_dnds):
            ax8.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                     f"{val:.3f}", ha='center', va='bottom', fontsize=8)
        ax8.set_xlabel("Gene", fontsize=11)
        ax8.set_ylabel("Mean pairwise dN/dS (yn00)", fontsize=11)
        ax8.set_title(f"Pairwise dN/dS per Mitochondrial Gene — {SPECIES}",
                      fontsize=12, fontweight="bold")
        ax8.legend()
        ax8.spines["top"].set_visible(False); ax8.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        dnds_fig_b64 = encode_fig(fig8)
        print("  dN/dS figure done")
    except Exception as e:
        print(f"  dN/dS figure failed: {e}")

# ── Figure 9: codeml heatmap ──────────────────────────────────────────────────
print("Generating codeml heatmap...")
codeml_fig_b64 = None
if codeml_rows:
    try:
        genes = sorted(set(r['gene'] for r in codeml_rows))
        pops  = sorted(set(r['population'] for r in codeml_rows))

        omega_mat = np.full((len(pops), len(genes)), np.nan)
        sig_mat   = np.zeros((len(pops), len(genes)), dtype=bool)
        for r in codeml_rows:
            gi = genes.index(r['gene'])
            pi = pops.index(r['population'])
            try:
                omega_mat[pi, gi] = float(r['omega_foreground'])
                sig_mat[pi, gi]   = r['p<0.05'] == 'YES'
            except (ValueError, KeyError):
                pass

        fig9, ax9 = plt.subplots(figsize=(max(10, len(genes)*1.2),
                                          max(5, len(pops)*0.8)))
        im = ax9.imshow(omega_mat, aspect='auto', cmap='RdYlBu_r', vmin=0, vmax=1)
        plt.colorbar(im, ax=ax9, label='Foreground ω (dN/dS)')
        ax9.set_xticks(range(len(genes)))
        ax9.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
        ax9.set_yticks(range(len(pops)))
        ax9.set_yticklabels(pops, fontsize=9)
        ax9.set_title(f"codeml Branch Model — Foreground ω per Population × Gene\n"
                      f"* = significant (LRT p<0.05)", fontsize=11, fontweight="bold")
        for pi in range(len(pops)):
            for gi in range(len(genes)):
                if sig_mat[pi, gi]:
                    ax9.text(gi, pi, '*', ha='center', va='center',
                             fontsize=14, fontweight='bold', color='white')
        plt.tight_layout()
        codeml_fig_b64 = encode_fig(fig9)
        print("  codeml heatmap done")
    except Exception as e:
        print(f"  codeml heatmap failed: {e}")

# ── Figure 9: BSP population size plots ──────────────────────────────────────
print("Generating BSP figures...")
bsp_plots_b64 = {}
if bsp_pop_data:
    POP_COLOURS = {"SI":"#1976D2","NS":"#388E3C","CAI":"#F57C00",
                   "SB":"#7B1FA2","COI":"#D32F2F","SP":"#0097A7",
                   "DH":"#5D4037","all":"#455A64"}
    for pop_label, data in bsp_pop_data.items():
        pop_sizes = data["pop_sizes"]
        if not pop_sizes: continue
        n_groups = len(pop_sizes[0])
        import numpy as np
        ps_arr = np.array(pop_sizes)  # shape: (n_samples, n_groups)
        medians = np.median(ps_arr, axis=0)
        lo95    = np.percentile(ps_arr, 2.5,  axis=0)
        hi95    = np.percentile(ps_arr, 97.5, axis=0)
        groups  = np.arange(1, n_groups + 1)
        colour  = POP_COLOURS.get(pop_label, "#607D8B")
        try:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.fill_between(groups, lo95, hi95, alpha=0.25, color=colour)
            ax.plot(groups, medians, color=colour, linewidth=2)
            ax.set_yscale("log")
            ax.set_xlabel("Coalescent interval (recent → ancient)", fontsize=10)
            ax.set_ylabel("Effective population size (Ne, log scale)", fontsize=10)
            n_samples = bsp_pop_summary.get("included_pops", {}).get(pop_label,
                        bsp_pop_summary.get("all_pops", {}).get(pop_label, ""))
            title_pop = "All samples" if pop_label == "all" else f"Population {pop_label}"
            n_label   = f" (n={n_samples})" if n_samples else ""
            ax.set_title(f"BSP — {title_pop}{n_label}", fontsize=11)
            ax.tick_params(labelsize=9)
            fig.tight_layout()
            bsp_plots_b64[pop_label] = encode_fig(fig)
            plt.close(fig)
            print(f"  BSP plot done for {pop_label}")
        except Exception as e:
            print(f"  BSP plot failed for {pop_label}: {e}")

# ── Figure 10: BSP ESS bar chart ─────────────────────────────────────────────
print("Generating BSP ESS figure...")
bsp_ess_fig_b64 = None
if bsp_ess_rows:
    try:
        runs = [r['run'].replace('BSP_','') for r in bsp_ess_rows]
        ess  = [float(r['ESS_posterior']) for r in bsp_ess_rows]
        cols = ['#4CAF50' if e >= 200 else '#F44336' for e in ess]

        fig10, ax10 = plt.subplots(figsize=(10, 4))
        bars = ax10.bar(runs, ess, color=cols, edgecolor='white')
        ax10.axhline(200, color='black', linestyle='--', linewidth=1,
                     label='ESS = 200 (minimum threshold)')
        for bar, val in zip(bars, ess):
            ax10.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                      f"{val:.0f}", ha='center', va='bottom', fontsize=8)
        ax10.set_xlabel("Run", fontsize=11)
        ax10.set_ylabel("ESS (posterior)", fontsize=11)
        ax10.set_title("BEAST2 ESS — Convergence Check per BSP Run",
                       fontsize=12, fontweight="bold")
        ax10.legend()
        ax10.spines["top"].set_visible(False); ax10.spines["right"].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        bsp_ess_fig_b64 = encode_fig(fig10)
        print("  BSP ESS figure done")
    except Exception as e:
        print(f"  BSP ESS figure failed: {e}")

# ── Build excluded-population banner ─────────────────────────────────────────
excl_banner = ""
if excl_display:
    excl_items = ", ".join(f"<b>{p}</b> (n={n})" for p, n in sorted(excl_display.items()))
    if IS_FILTERED:
        excl_banner = f"""
        <div style="background:#FFF8E1;border-left:4px solid #FF9800;padding:12px 16px;
                    border-radius:4px;margin-bottom:16px;font-size:0.9em;">
          ⚠ <b>Population filter applied (n &ge; {MIN_POP_SIZE}).</b>
          The following populations were excluded from per-population analyses
          (codeml branch model, BSP): {excl_items}.
          They are still included in pairwise yn00 dN/dS and the combined BSP run.
          See <code>mito_report_all.html</code> for the unfiltered report.
        </div>"""
    else:
        excl_banner = f"""
        <div style="background:#E3F2FD;border-left:4px solid #1565C0;padding:12px 16px;
                    border-radius:4px;margin-bottom:16px;font-size:0.9em;">
          ℹ <b>Small populations detected.</b>
          The following populations have fewer samples than the pipeline filter threshold
          and were excluded from per-population analyses during pipeline execution: {excl_items}.
          See <code>mito_report_filtered.html</code> for the report with these excluded.
        </div>"""

# ── Table: excluded populations ───────────────────────────────────────────────
excl_table_rows = ""
for pop, n in sorted(excl_display.items()):
    excl_table_rows += f"""
    <tr>
        <td><b>{pop}</b></td>
        <td>{n}</td>
        <td>Below min_pop_size threshold — excluded from codeml branch model and BSP</td>
    </tr>"""

# ── QC table rows ─────────────────────────────────────────────────────────────
qc_rows = ""
for s in sorted(samples, key=lambda x: x.get("sample","")):
    status = s.get("status","")
    css = "pass" if "PASS" in status else ("warn" if "WARN" in status else "fail")
    qc_rows += f"""
    <tr class="{css}">
        <td>{s.get('sample','')}</td>
        <td>{int(float(s.get('mapped_reads',0))):,}</td>
        <td>{float(s.get('mean_depth',0)):.1f}x</td>
        <td>{float(s.get('coverage_breadth_pct',0)):.1f}%</td>
        <td>{s.get('num_variants_pass','')}</td>
        <td>{int(s.get('consensus_length_bp',0)):,} bp</td>
        <td>{s.get('N_count','0')}</td>
        <td><b>{status}</b></td>
    </tr>"""

# ── Annotation table rows ─────────────────────────────────────────────────────
annot_rows = ""
for g in annot_genes:
    type_col = "#2196F3" if g['type'] == 'CDS' else "#4CAF50"
    annot_rows += f"""
    <tr>
        <td><b>{g['gene']}</b></td>
        <td>{g['product']}</td>
        <td><span style="color:{type_col};font-weight:bold">{g['type']}</span></td>
        <td>{g['start']}–{g['end']}</td>
        <td>{g['length_bp']} bp</td>
        <td>{'&minus;' if g['direction'] == '-' else '+'}</td>
    </tr>"""

# ── yn00 table rows ───────────────────────────────────────────────────────────
yn00_table_rows = ""
for r in yn00_rows:
    try:
        omega = float(r['mean_dNdS'])
        col = "#F44336" if omega > 1 else ("#FF9800" if omega > 0.5 else "#4CAF50")
        yn00_table_rows += f"""
        <tr>
            <td><b>{r['gene']}</b></td>
            <td>{r['n_pairs']}</td>
            <td>{float(r['mean_dN']):.5f}</td>
            <td>{float(r['mean_dS']):.5f}</td>
            <td style="color:{col};font-weight:bold">{omega:.5f}</td>
            <td>{float(r['min_dNdS']):.5f}</td>
            <td>{float(r['max_dNdS']):.5f}</td>
        </tr>"""
    except (ValueError, KeyError):
        pass

# ── codeml table rows ─────────────────────────────────────────────────────────
codeml_table_rows = ""
for r in codeml_rows:
    sig   = r.get('p<0.05','no')
    css   = 'pass' if sig == 'YES' else ''
    badge = "<span style='color:#F44336;font-weight:bold'>YES ✱</span>" \
            if sig == 'YES' else "no"
    n_s   = r.get('n_samples','')
    codeml_table_rows += f"""
    <tr class="{css}">
        <td><b>{r.get('gene','')}</b></td>
        <td>{r.get('population','')} <span style="color:#999;font-size:0.85em">(n={n_s})</span></td>
        <td>{r.get('lnL_null','')}</td>
        <td>{r.get('lnL_branch','')}</td>
        <td>{r.get('LRT','')}</td>
        <td>{badge}</td>
        <td>{r.get('omega_background','')}</td>
        <td>{r.get('omega_foreground','')}</td>
    </tr>"""

# ── BSP ESS table rows ────────────────────────────────────────────────────────
bsp_ess_table_rows = ""
for r in bsp_ess_rows:
    css = 'pass' if r.get('status','').startswith('OK') else 'fail'
    bsp_ess_table_rows += f"""
    <tr class="{css}">
        <td>{r.get('run','').replace('BSP_','')}</td>
        <td>{r.get('ESS_posterior','')}</td>
        <td><b>{r.get('status','')}</b></td>
    </tr>"""

# ── Methods text (dynamic) ────────────────────────────────────────────────────
# Genetic code description
gc_descriptions = {
    1: "standard", 2: "vertebrate mitochondrial", 4: "mold/protozoan mitochondrial",
    5: "invertebrate mitochondrial", 9: "echinoderm mitochondrial",
    13: "ascidian mitochondrial", 14: "alternative flatworm mitochondrial"
}
gc_desc = gc_descriptions.get(GENETIC_CODE, f"code {GENETIC_CODE}")

# Clock rate description
clock_desc = f"U({CLOCK_MIN}, {CLOCK_MAX}) subs/site/Myr"

# Filtered report subtitle
report_subtitle = ""
if IS_FILTERED:
    report_subtitle = f' <span style="background:#FFF8E1;color:#E65100;padding:2px 8px;border-radius:10px;font-size:0.7em;font-weight:bold;">FILTERED (n≥{MIN_POP_SIZE})</span>'
else:
    report_subtitle = f' <span style="background:#E8F5E9;color:#2E7D32;padding:2px 8px;border-radius:10px;font-size:0.7em;font-weight:bold;">ALL POPULATIONS</span>'

# ── Write HTML ────────────────────────────────────────────────────────────────
print("Writing HTML...")
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Mitogenome Report — {SPECIES}{"  [filtered]" if IS_FILTERED else ""}</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 0;
          background: #f4f6f9; color: #222; }}
  .header {{ background: linear-gradient(135deg, #1A237E, #3949AB);
             color: white; padding: 40px; text-align: center; }}
  .header h1 {{ margin: 0; font-size: 2em; }}
  .header p  {{ margin: 8px 0 0; opacity: 0.85; font-size: 1.1em; }}
  .container {{ max-width: 1200px; margin: 30px auto; padding: 0 20px; }}
  .card {{ background: white; border-radius: 8px; padding: 24px;
           margin-bottom: 24px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
  h2 {{ color: #1A237E; border-bottom: 2px solid #E8EAF6;
        padding-bottom: 8px; margin-top: 0; }}
  h3 {{ color: #37474F; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(4, 1fr);
                   gap: 16px; margin-bottom: 8px; }}
  .stat-box {{ background: #F5F7FF; border-radius: 6px; padding: 16px;
               text-align: center; border-left: 4px solid #3949AB; }}
  .stat-box.pass {{ border-left-color: #4CAF50; }}
  .stat-box.warn {{ border-left-color: #FF9800; }}
  .stat-box.fail {{ border-left-color: #F44336; }}
  .stat-box .val {{ font-size: 2em; font-weight: bold; color: #1A237E; }}
  .stat-box.pass .val {{ color: #2E7D32; }}
  .stat-box.warn .val {{ color: #E65100; }}
  .stat-box.fail .val {{ color: #B71C1C; }}
  .stat-box .lbl {{ font-size: 0.85em; color: #666; margin-top: 4px; }}
  .plots-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  img {{ max-width: 100%; border-radius: 4px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; }}
  th {{ background: #1A237E; color: white; padding: 8px 10px;
        text-align: left; position: sticky; top: 0; }}
  td {{ border-bottom: 1px solid #eee; padding: 6px 10px; }}
  .pass td {{ background: #F1F8E9; }}
  .warn td {{ background: #FFF8E1; }}
  .fail td {{ background: #FFEBEE; }}
  .table-wrap {{ max-height: 500px; overflow-y: auto; border-radius: 6px;
                 border: 1px solid #ddd; }}
  pre {{ background: #F5F5F5; padding: 16px; border-radius: 6px;
         font-size: 0.8em; overflow-x: auto; white-space: pre-wrap; }}
  .methods p {{ line-height: 1.7; text-align: justify; }}
  .footer {{ text-align: center; color: #999; font-size: 0.8em; padding: 20px; }}
  .muted {{ color: #999; font-style: italic; }}
  .nav {{ background: white; border-bottom: 1px solid #ddd; padding: 12px 20px;
          position: sticky; top: 0; z-index: 100;
          box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
  .nav a {{ color: #3949AB; text-decoration: none; margin-right: 20px;
            font-size: 0.9em; font-weight: 500; }}
  .nav a:hover {{ text-decoration: underline; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 10px;
            font-size: 0.8em; font-weight: bold; }}
  .badge-blue  {{ background:#E3F2FD; color:#1565C0; }}
  .badge-green {{ background:#E8F5E9; color:#2E7D32; }}
</style>
</head>
<body>

<div class="header">
  <h1>Mitogenome Extraction Pipeline{report_subtitle}</h1>
  <p><i>{SPECIES}</i> &nbsp;|&nbsp; {datetime.now().strftime("%d %B %Y, %H:%M")}
     &nbsp;|&nbsp; {len(POP_ORDER)} populations
     {"&nbsp;|&nbsp; min pop size: " + str(MIN_POP_SIZE) if IS_FILTERED else ""}
  </p>
</div>

<div class="nav">
  <a href="#summary">Summary</a>
  <a href="#coverage">Coverage</a>
  <a href="#qc-table">Per-Sample QC</a>
  <a href="#consensus">Consensus</a>
  <a href="#annotation">Annotation</a>
  <a href="#phylo">Phylogenetics</a>
  <a href="#dnds">dN/dS</a>
  <a href="#bsp">Bayesian Skyline</a>
  <a href="#methods">Methods</a>
  <a href="#outputs">Output Files</a>
</div>

<div class="container">

{excl_banner}

<!-- ── Run Summary ──────────────────────────────────────────────────────── -->
<div class="card" id="summary">
  <h2>Run Summary</h2>
  <div class="summary-grid">
    <div class="stat-box">
      <div class="val">{n_total}</div><div class="lbl">Total Samples</div></div>
    <div class="stat-box pass">
      <div class="val">{n_pass}</div>
      <div class="lbl">PASS ({n_pass/n_total*100:.0f}%)</div></div>
    <div class="stat-box warn">
      <div class="val">{n_warn