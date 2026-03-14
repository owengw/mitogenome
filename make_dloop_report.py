#!/usr/bin/env python3
# =============================================================================
# make_dloop_report.py — D-loop haplotype network + NJ tree HTML report
#
# Generates an interactive HTML report containing:
#   Tab 1: TCS-style haplotype network (all populations, colour-coded)
#   Tab 2: Unrooted NJ tree coloured by population
#
# Usage:
#   python3 make_dloop_report.py \
#       --aligned dloop_aligned.fa \
#       --treefile dloop_nj.treefile \
#       --output dloop_report.html \
#       --pop-defs SI:SI:#3b82f6 NS:NS:#22c55e CAI:CAI:#f97316 SB:SB:#a855f7
# =============================================================================

import argparse
import copy
import math
import re
import sys
from collections import defaultdict, Counter
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--aligned",  required=True)
parser.add_argument("--treefile", required=True)
parser.add_argument("--output",   required=True)
parser.add_argument("--pop-defs", nargs="+", default=[])
args = parser.parse_args()

# ── Population colour map ─────────────────────────────────────────────────────
pop_colours = {}
pop_patterns = []
for pd in args.pop_defs:
    parts = pd.split(":")
    if len(parts) >= 2:
        pattern = parts[0]
        pop     = parts[1]
        colour  = parts[2] if len(parts) > 2 else "#94a3b8"
        pop_colours[pop] = colour
        pop_patterns.append((pattern, pop, colour))

def get_pop(sample_id):
    sid = sample_id.split("|")[1] if "|" in sample_id else sample_id
    for pattern, pop, colour in pop_patterns:
        if pattern.upper() in sample_id.upper():
            return pop, colour
    return "UNK", "#94a3b8"

# ── Load aligned sequences ────────────────────────────────────────────────────
records = list(SeqIO.parse(args.aligned, "fasta"))
print(f"Loaded {len(records)} aligned sequences")

# ── Build haplotypes ──────────────────────────────────────────────────────────
# Group identical sequences; track which populations carry each haplotype
hap_seqs  = {}   # hap_id -> sequence string
hap_pops  = defaultdict(lambda: defaultdict(int))  # hap_id -> {pop: count}
seq_to_hap = {}

hap_counter = 0
for rec in records:
    seq = str(rec.seq).upper()
    pop, colour = get_pop(rec.id)
    if seq not in seq_to_hap:
        hap_counter += 1
        hid = f"H{hap_counter}"
        seq_to_hap[seq] = hid
        hap_seqs[hid]   = seq
    hid = seq_to_hap[seq]
    hap_pops[hid][pop] += 1

print(f"Unique haplotypes: {hap_counter}")

# ── Compute pairwise distances ────────────────────────────────────────────────
def hamming(s1, s2):
    """Count pairwise differences (ignoring gaps and Ns)."""
    diffs = 0
    compared = 0
    for a, b in zip(s1, s2):
        if a in "-N" or b in "-N":
            continue
        compared += 1
        if a != b:
            diffs += 1
    return diffs

hap_ids = sorted(hap_seqs.keys(), key=lambda x: int(x[1:]))
n_haps  = len(hap_ids)

dist = {}
for i, h1 in enumerate(hap_ids):
    for j, h2 in enumerate(hap_ids):
        if i < j:
            d = hamming(hap_seqs[h1], hap_seqs[h2])
            dist[(h1, h2)] = d
            dist[(h2, h1)] = d
        elif i == j:
            dist[(h1, h2)] = 0

# ── TCS network (parsimony connection limit) ──────────────────────────────────
# Connect haplotypes that differ by <= epsilon steps (TCS-style)
# Use minimum spanning network approach: connect nearest neighbours

def mst_network(hap_ids, dist):
    """
    Build minimum spanning tree connecting all haplotypes.
    Returns list of (h1, h2, distance) edges.
    """
    if len(hap_ids) <= 1:
        return []

    in_tree  = {hap_ids[0]}
    edges    = []

    while len(in_tree) < len(hap_ids):
        best = None
        for h1 in in_tree:
            for h2 in hap_ids:
                if h2 in in_tree:
                    continue
                d = dist.get((h1, h2), 9999)
                if best is None or d < best[2]:
                    best = (h1, h2, d)
        if best:
            edges.append(best)
            in_tree.add(best[1])

    # Also add edges within the parsimony limit (max 3 steps) to show loops
    parsimony_limit = 3
    for i, h1 in enumerate(hap_ids):
        for j, h2 in enumerate(hap_ids):
            if i >= j:
                continue
            d = dist.get((h1, h2), 9999)
            if d <= parsimony_limit:
                edge = (h1, h2, d)
                if edge not in edges and (h2, h1, d) not in edges:
                    edges.append(edge)

    return edges

network_edges = mst_network(hap_ids, dist)
print(f"Network edges: {len(network_edges)}")

# ── Layout network using force-directed placement ─────────────────────────────
def layout_network(hap_ids, edges, width=800, height=600, iterations=200):
    """
    Simple force-directed layout (Fruchterman-Reingold style).
    """
    import random
    random.seed(42)

    pos = {h: [random.uniform(100, width-100),
               random.uniform(100, height-100)] for h in hap_ids}

    k   = math.sqrt((width * height) / max(len(hap_ids), 1))
    dt  = 0.1

    for _ in range(iterations):
        disp = {h: [0.0, 0.0] for h in hap_ids}

        # Repulsive forces between all pairs
        for i, h1 in enumerate(hap_ids):
            for j, h2 in enumerate(hap_ids):
                if i >= j:
                    continue
                dx = pos[h1][0] - pos[h2][0]
                dy = pos[h1][1] - pos[h2][1]
                d  = max(math.sqrt(dx*dx + dy*dy), 0.01)
                fr = (k * k) / d
                disp[h1][0] += (dx / d) * fr
                disp[h1][1] += (dy / d) * fr
                disp[h2][0] -= (dx / d) * fr
                disp[h2][1] -= (dy / d) * fr

        # Attractive forces along edges
        for h1, h2, d_val in edges:
            dx = pos[h1][0] - pos[h2][0]
            dy = pos[h1][1] - pos[h2][1]
            d  = max(math.sqrt(dx*dx + dy*dy), 0.01)
            # Stronger attraction for closer haplotypes
            fa = (d * d) / k * (1 + d_val * 0.5)
            disp[h1][0] -= (dx / d) * fa
            disp[h1][1] -= (dy / d) * fa
            disp[h2][0] += (dx / d) * fa
            disp[h2][1] += (dy / d) * fa

        # Apply displacements with cooling
        temp = k * (1 - _ / iterations)
        for h in hap_ids:
            d  = max(math.sqrt(disp[h][0]**2 + disp[h][1]**2), 0.01)
            mv = min(d, temp)
            pos[h][0] += (disp[h][0] / d) * mv * dt
            pos[h][1] += (disp[h][1] / d) * mv * dt
            # Clamp to bounds
            pos[h][0] = max(40, min(width-40,  pos[h][0]))
            pos[h][1] = max(40, min(height-40, pos[h][1]))

    return pos

pos = layout_network(hap_ids, network_edges)

# ── Total count per haplotype (for node sizing) ───────────────────────────────
hap_totals = {h: sum(hap_pops[h].values()) for h in hap_ids}
max_total  = max(hap_totals.values()) if hap_totals else 1

# ── Build network SVG ─────────────────────────────────────────────────────────
def build_network_svg(hap_ids, hap_pops, hap_totals, pos, edges,
                      pop_colours, width=820, height=620):
    """
    Render haplotype network as SVG.
    Node size proportional to haplotype frequency.
    Pie chart sectors show population composition per haplotype.
    Edge labels show number of mutations.
    """
    parts = []

    # Edges
    for h1, h2, d_val in edges:
        x1, y1 = pos[h1]
        x2, y2 = pos[h2]
        # Draw tick marks on edge for each mutation step
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        colour = "#334155" if d_val <= 1 else "#475569"
        width_attr = "1.5" if d_val <= 1 else "1"
        parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{colour}" stroke-width="{width_attr}"/>')
        if d_val > 1:
            parts.append(
                f'<text x="{mid_x:.1f}" y="{mid_y:.1f}" font-size="8" '
                f'fill="#64748b" text-anchor="middle">{d_val}</text>')

    # Nodes — pie chart per haplotype
    for h in hap_ids:
        x, y      = pos[h]
        total     = hap_totals[h]
        r         = max(6, min(22, 5 + math.sqrt(total) * 3))
        pops_here = hap_pops[h]

        if len(pops_here) == 1:
            # Solid circle
            pop_name = list(pops_here.keys())[0]
            col = pop_colours.get(pop_name, "#94a3b8")
            parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" '
                f'fill="{col}" stroke="#0f172a" stroke-width="1.5"/>')
        else:
            # Pie sectors
            sorted_pops = sorted(pops_here.items(), key=lambda kv: -kv[1])
            angle = -math.pi / 2
            for pop_name, count in sorted_pops:
                col     = pop_colours.get(pop_name, "#94a3b8")
                frac    = count / total
                sweep   = frac * 2 * math.pi
                x1s     = x + r * math.cos(angle)
                y1s     = y + r * math.sin(angle)
                angle  += sweep
                x2s     = x + r * math.cos(angle)
                y2s     = y + r * math.sin(angle)
                large   = 1 if sweep > math.pi else 0
                parts.append(
                    f'<path d="M{x:.1f},{y:.1f} L{x1s:.2f},{y1s:.2f} '
                    f'A{r:.1f},{r:.1f} 0 {large},1 {x2s:.2f},{y2s:.2f} Z" '
                    f'fill="{col}" stroke="#0f172a" stroke-width="1"/>')

        # Haplotype label
        label_y = y - r - 3
        fs = 8 if total > 2 else 7
        fw = "bold" if total > 5 else "normal"
        parts.append(
            f'<text x="{x:.1f}" y="{label_y:.1f}" font-size="{fs}" '
            f'fill="#e2e8f0" text-anchor="middle" font-weight="{fw}">{h}</text>')
        if total > 1:
            parts.append(
                f'<text x="{x:.1f}" y="{y+r+10:.1f}" font-size="7" '
                f'fill="#94a3b8" text-anchor="middle">(n={total})</text>')

    return (f'<svg viewBox="0 0 {width} {height}" '
            f'style="width:100%;display:block">\n' +
            '\n'.join(parts) + '\n</svg>')

network_svg = build_network_svg(
    hap_ids, hap_pops, hap_totals, pos, network_edges, pop_colours)

# ── Parse + render NJ tree ────────────────────────────────────────────────────
with open(args.treefile) as f:
    newick = f.read().strip()

def parse_newick(s):
    ancestors, tree = [], {}
    tokens = re.split(r'\s*(;|\(|\)|,|:)\s*', s)
    subtree = tree
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == '(':
            subtree['children'] = [{}]
            ancestors.append(subtree)
            subtree = subtree['children'][0]
        elif tok == ',':
            parent = ancestors[-1]
            parent['children'].append({})
            subtree = parent['children'][-1]
        elif tok == ')':
            subtree = ancestors.pop()
        elif tok == ':':
            if i + 1 < len(tokens):
                try:    subtree['length'] = float(tokens[i + 1])
                except: pass
                i += 1
        elif tok and tok != ';':
            if 'name' not in subtree:
                subtree['name'] = tok.replace('_', ' ')
        i += 1
    return tree

def all_leaves(node):
    if 'children' not in node: return [node]
    out = []
    for c in node['children']: out.extend(all_leaves(c))
    return out

def count_leaves(node):
    return len(all_leaves(node))

def assign_x_bl(node, depth=0):
    """Assign x position using branch lengths."""
    bl = node.get('length', 0.01)
    node['_x'] = depth + bl
    for c in node.get('children', []):
        assign_x_bl(c, node['_x'])

def assign_y(node, counter=None):
    if counter is None: counter = [0]
    if 'children' not in node:
        node['_y'] = counter[0]; counter[0] += 1
    else:
        for c in node['children']: assign_y(c, counter)
        node['_y'] = sum(c['_y'] for c in node['children']) / len(node['children'])
    return counter[0]

def render_nj_tree(tree):
    t = copy.deepcopy(tree)
    assign_x_bl(t, 0)
    total = assign_y(t)

    max_x = max(l['_x'] for l in all_leaves(t))

    PAD_L, PAD_R, PAD_T = 20, 200, 20
    SVG_W  = 900
    ROW_H  = max(6, min(14, 500 // max(total, 1)))
    SVG_H  = max(300, total * ROW_H + PAD_T * 2)
    TREE_W = SVG_W - PAD_L - PAD_R

    def px(node): return PAD_L + (node['_x'] / max(max_x, 0.001)) * TREE_W
    def py(node): return PAD_T + (node['_y'] / max(total - 1, 1)) * (SVG_H - PAD_T * 2)

    lines, labels = [], []

    def walk(node):
        x1, y1 = px(node), py(node)
        if 'children' in node:
            for c in node['children']:
                x2, y2 = px(c), py(c)
                lines.append(f'<line x1="{x1:.1f}" y1="{y2:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#334155" stroke-width="1"/>')
                lines.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x1:.1f}" y2="{y2:.1f}" stroke="#334155" stroke-width="1"/>')
                walk(c)
        else:
            name     = node.get('name', '')
            # Parse pop from name (format: sampleID|POP)
            pop, col = get_pop(name)
            label    = name.split("|")[0] if "|" in name else name
            label    = label[:25]
            labels.append(f'<circle cx="{x1+2:.1f}" cy="{y1:.1f}" r="3" fill="{col}"/>')
            labels.append(f'<text x="{x1+8:.1f}" y="{y1+3:.1f}" font-size="8" fill="{col}" font-style="italic">{label}</text>')

    walk(t)
    content = '\n'.join(lines + labels)
    return (f'<svg viewBox="0 0 {SVG_W} {SVG_H}" style="width:100%;display:block">\n'
            f'{content}\n</svg>', total)

tree_parsed = parse_newick(newick)
tree_svg, n_leaves = render_nj_tree(tree_parsed)
print(f"NJ tree: {n_leaves} leaves")

# ── Legend ────────────────────────────────────────────────────────────────────
legend_html = ""
for pop, col in pop_colours.items():
    legend_html += (
        f'<span style="display:inline-flex;align-items:center;gap:5px;margin:3px 10px 3px 0">'
        f'<span style="width:12px;height:12px;border-radius:50%;background:{col};flex-shrink:0"></span>'
        f'<span style="color:#94a3b8;font-size:12px">{pop}</span></span>'
    )

# Haplotype frequency summary table
hap_table_rows = ""
for h in hap_ids:
    pops_str = "  ".join(f'<span style="color:{pop_colours.get(p,"#94a3b8")}">{p}:{c}</span>'
                          for p, c in sorted(hap_pops[h].items()))
    hap_table_rows += (
        f'<tr><td style="color:#e2e8f0;padding:2px 8px">{h}</td>'
        f'<td style="color:#94a3b8;padding:2px 8px">{hap_totals[h]}</td>'
        f'<td style="padding:2px 8px">{pops_str}</td></tr>'
    )

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>D-loop Analysis — P. longirostris</title>
<style>
  * {{ box-sizing:border-box; margin:0; padding:0 }}
  body {{ background:#0f172a; color:#e2e8f0; font-family:'Courier New',monospace; padding:20px }}
  h1 {{ color:#3b82f6; font-size:20px; letter-spacing:1px; margin-bottom:4px }}
  .subtitle {{ color:#475569; font-size:12px; margin-bottom:16px }}
  .legend {{ margin-bottom:12px }}
  .tab-bar {{ display:flex; gap:8px; margin-bottom:12px }}
  .tab {{ padding:6px 16px; border-radius:4px; cursor:pointer; font-size:12px;
           border:1px solid #334155; background:#1e293b; color:#94a3b8 }}
  .tab.active {{ background:#3b82f6; border-color:#3b82f6; color:#000; font-weight:bold }}
  .tree-wrap {{ background:#0a1628; border:1px solid #1e293b; border-radius:8px;
                overflow:auto; max-height:80vh; padding:10px }}
  .tree-panel {{ display:none }}
  .tree-panel.active {{ display:block }}
  table {{ border-collapse:collapse; font-size:11px; margin-top:10px }}
  th {{ color:#64748b; padding:3px 8px; border-bottom:1px solid #1e293b; text-align:left }}
  .stats {{ background:#0a1628; border:1px solid #1e293b; border-radius:6px;
             padding:12px; margin-bottom:12px; font-size:12px }}
  .stats span {{ color:#3b82f6 }}
</style>
</head>
<body>
<h1>D-loop Analysis — Plestiodon longirostris</h1>
<p class="subtitle">Control region · {len(records)} sequences · {hap_counter} unique haplotypes · unrooted NJ tree</p>

<div class="stats">
  Sequences: <span>{len(records)}</span> &nbsp;|&nbsp;
  Haplotypes: <span>{hap_counter}</span> &nbsp;|&nbsp;
  Network edges: <span>{len(network_edges)}</span>
</div>

<div class="legend">{legend_html}</div>

<div class="tab-bar">
  <div class="tab active" onclick="show('network',this)">Haplotype Network</div>
  <div class="tab" onclick="show('tree',this)">NJ Tree</div>
  <div class="tab" onclick="show('table',this)">Haplotype Table</div>
</div>

<div class="tree-wrap">
  <div id="panel-network" class="tree-panel active">
    <p style="font-size:11px;color:#475569;margin-bottom:8px">
      TCS-style parsimony network. Node size ∝ frequency. 
      Pie charts show population composition. Edge numbers = mutation steps.
    </p>
    {network_svg}
  </div>
  <div id="panel-tree" class="tree-panel">
    <p style="font-size:11px;color:#475569;margin-bottom:8px">
      Unrooted NJ tree. Branch lengths proportional to substitutions.
    </p>
    {tree_svg}
  </div>
  <div id="panel-table" class="tree-panel">
    <table>
      <thead><tr>
        <th>Haplotype</th><th>Total n</th><th>By population</th>
      </tr></thead>
      <tbody>{hap_table_rows}</tbody>
    </table>
  </div>
</div>

<script>
function show(which,el){{
  document.querySelectorAll('.tree-panel').forEach(p=>p.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
  document.getElementById('panel-'+which).classList.add('active');
  el.classList.add('active');
}}
</script>
</body>
</html>"""

with open(args.output, "w") as f:
    f.write(html)

print(f"HTML report written → {args.output}")
print(f"  Network: {hap_counter} haplotypes, {len(network_edges)} edges")
print(f"  NJ tree: {n_leaves} leaves")