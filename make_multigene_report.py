#!/usr/bin/env python3
# =============================================================================
# make_multigene_report.py — HTML report for multigene phylogeny
#
# Displays three trees in tabs:
#   1. MrBayes consensus tree (partitioned 12S + ND4 + CYTB, posterior probs)
#   2. IQ-TREE ML tree (concatenated alignment, bootstrap)
#   3. NJ tree (CYTB only)
#
# Usage:
#   python3 make_multigene_report.py \
#       --outdir /path/to/multigene \
#       --mrbayes-tree mrbayes_output.con.tre \
#       --ml-tree iqtree_concat.treefile \
#       --nj-tree nj_cytb.treefile \
#       --output multigene_report.html
# =============================================================================

import argparse
import copy
import os
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--outdir",       required=True)
parser.add_argument("--mrbayes-tree", default=None)
parser.add_argument("--ml-tree",      default=None)
parser.add_argument("--nj-tree",      default=None)
parser.add_argument("--output",       required=True)
args = parser.parse_args()

# ── Colour / label helpers ────────────────────────────────────────────────────
OUTGROUP_KEYWORDS = [
    "eumeces", "scincella", "chalcides", "trachylepis",
    "mabuya", "lygosoma", "mesoscincus", "eurylepis",
]

def get_color(name):
    if not name: return "#94a3b8"
    n = name.lower()
    if "longirostris" in n:                                            return "#f97316"
    if "reynoldsi"    in n:                                            return "#a78bfa"
    if any(x in n for x in ("fasciatus","laticeps","inexpectatus",
        "anthracinus","multivirgatus","obsoletus","septentrionalis",
        "skiltonianus","tetragrammus","egregius","gilberti")):         return "#34d399"
    if any(x in n for x in ("barbouri","finitimus","japonicus")):      return "#60a5fa"
    if any(x in n for x in ("elegans","takarai","kuchinoshimensis",
        "stimpsonii","marginatus","oshimensis","latiscutatus",
        "leucostictus","kishinouyei","tamdaoensis","chinensis",
        "quadrilineatus","liui","coreensis","tunganus","popei",
        "bilineatus")):                                                return "#38bdf8"
    if any(x in n for x in ("brevirostris","callicephalus","capito",
        "colimensis","copei","dicei","dugesii","indubitus",
        "lagunensis","longiartus","lotus","lynxe","multilineatus",
        "nietoi","ochoteranae","parviauriculatus","parvulus",
        "sumichrasti")):                                               return "#4ade80"
    if any(x in n for x in OUTGROUP_KEYWORDS):                        return "#f43f5e"
    return "#94a3b8"

def get_label(name):
    if not name: return ""
    if "longirostris" in name.lower():
        # e.g. Plestiodon_longirostris_SI -> P. longirostris (SI)
        parts = name.replace("_", " ").split()
        pop = parts[-1] if len(parts) > 2 else ""
        return f"P. longirostris ({pop})" if pop else "P. longirostris"
    name = name.replace("_", " ")
    parts = name.split()
    known = ("Plestiodon","Eumeces","Scincella","Chalcides","Trachylepis",
             "Mabuya","Lygosoma","Mesoscincus","Eurylepis")
    if len(parts) >= 2 and parts[0] in known:
        return f"{parts[0][0]}. {parts[1]}"
    return name[:30]

def is_outgroup(node):
    name = node.get("name","").lower()
    return "children" not in node and any(k in name for k in OUTGROUP_KEYWORDS)

# ── Newick parser ─────────────────────────────────────────────────────────────
def parse_newick(s):
    """
    Parse a newick string into a nested dict.
    Handles:
      - branch lengths after ':'
      - support/posterior values as node names after ')'
      - MrBayes consensus tree format (posterior probs as node labels)
    """
    # Strip comments [&...] which MrBayes includes in its consensus tree
    s = re.sub(r'\[&[^\]]*\]', '', s)
    s = re.sub(r'\[[^\]]*\]', '', s)
    s = s.strip().rstrip(";")

    ancestors, tree = [], {}
    subtree = tree
    i = 0
    tokens = re.split(r'\s*(;|\(|\)|,|:)\s*', s)

    while i < len(tokens):
        tok = tokens[i].strip()
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
                try:
                    subtree['length'] = float(tokens[i+1].strip())
                except (ValueError, IndexError):
                    pass
                i += 1
        elif tok and tok != ';':
            if 'name' not in subtree:
                subtree['name'] = tok
        i += 1
    return tree

# ── Tree utilities ─────────────────────────────────────────────────────────────
def all_leaves(node):
    if 'children' not in node:
        return [node]
    out = []
    for c in node['children']:
        out.extend(all_leaves(c))
    return out

def count_leaves(node):
    return len(all_leaves(node))

def og_fraction(node):
    leaves = all_leaves(node)
    return sum(1 for l in leaves if is_outgroup(l)) / len(leaves) if leaves else 0

def sort_by_outgroup(node):
    """Recursively sort children so outgroup branch renders at bottom."""
    if 'children' not in node:
        return node
    for c in node['children']:
        sort_by_outgroup(c)
    node['children'].sort(key=lambda c: og_fraction(c))
    return node

def collapse_longi(node):
    """Collapse all P. longirostris leaves into a single labelled node."""
    node = copy.deepcopy(node)
    if 'children' not in node:
        return node

    def all_longi(n):
        if 'children' not in n:
            return 'longirostris' in n.get('name','').lower()
        return all(all_longi(c) for c in n['children'])

    new_children = []
    for c in node['children']:
        if all_longi(c):
            n = count_leaves(c)
            new_children.append({
                'name': f'P. longirostris ({n} reps)',
                '_collapsed': True
            })
        else:
            new_children.append(collapse_longi(c))
    node['children'] = new_children
    return node

# ── Support value parser ──────────────────────────────────────────────────────
def parse_support(name):
    """
    Extract support value from node name.
    MrBayes: posterior probability as decimal (e.g. '0.95', '1.00')
    IQ-TREE: bootstrap as integer (e.g. '95', '100')
    Returns (display_string, float_value) or (None, None)
    """
    if not name:
        return None, None
    name = name.strip()
    try:
        v = float(name)
        if 0.0 <= v <= 1.0 and '.' in name:
            # MrBayes posterior probability
            return f"{v:.2f}", v
        elif 1.0 < v <= 100.0:
            # IQ-TREE bootstrap
            return f"{int(v)}", v / 100.0
        elif v == 1.0 and '.' not in name:
            return "100", 1.0
    except (ValueError, TypeError):
        pass
    return None, None

# ── SVG renderer ──────────────────────────────────────────────────────────────
def assign_x(node, depth=0):
    node['_x'] = depth
    for c in node.get('children', []):
        assign_x(c, depth + 1)

def assign_x_bl(node, depth=0):
    """Branch-length aware x assignment."""
    bl = node.get('length', 0.01)
    node['_x'] = depth + bl
    for c in node.get('children', []):
        assign_x_bl(c, node['_x'])

def assign_y(node, counter=None):
    if counter is None:
        counter = [0]
    if 'children' not in node:
        node['_y'] = counter[0]
        counter[0] += 1
    else:
        for c in node['children']:
            assign_y(c, counter)
        node['_y'] = sum(c['_y'] for c in node['children']) / len(node['children'])
    return counter[0]

def get_max_depth(node):
    if 'children' not in node:
        return node.get('_x', 0)
    return max(get_max_depth(c) for c in node['children'])

def render_tree(tree, svg_id, use_bl=False, support_type="bootstrap",
                show_all=True, collapse_longirostris=False):
    """
    Render a phylogenetic tree as SVG.

    support_type: 'bootstrap' (IQ-TREE, threshold 70/95)
                  'posterior' (MrBayes, threshold 0.70/0.95)
    use_bl:       use branch lengths for x positions
    """
    t = copy.deepcopy(tree)
    sort_by_outgroup(t)
    if collapse_longirostris:
        t = collapse_longi(t)

    if use_bl:
        assign_x_bl(t, 0)
    else:
        assign_x(t, 0)
    total = assign_y(t)
    md    = get_max_depth(t)

    PAD_L, PAD_R, PAD_T = 15, 240, 20
    SVG_W  = 1000
    ROW_H  = 7 if show_all else 18
    SVG_H  = max(400, total * ROW_H + PAD_T * 2)
    TREE_W = SVG_W - PAD_L - PAD_R
    TREE_H = SVG_H - PAD_T * 2

    def px(node):
        return PAD_L + (node['_x'] / max(md, 0.001)) * TREE_W

    def py(node):
        return PAD_T + (node['_y'] / max(total - 1, 1)) * TREE_H

    lines, dots, labels = [], [], []

    def walk(node):
        x1, y1 = px(node), py(node)

        if 'children' in node:
            for c in node['children']:
                x2, y2 = px(c), py(c)
                is_l   = 'longirostris' in c.get('name','').lower() or c.get('_collapsed')
                stroke = '#f9731640' if is_l else '#1e3a5f'
                sw     = '0.6' if show_all else '1.2'
                lines.append(
                    f'<line x1="{x1:.1f}" y1="{y2:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                    f'stroke="{stroke}" stroke-width="{sw}"/>')
                lines.append(
                    f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x1:.1f}" y2="{y2:.1f}" '
                    f'stroke="{stroke}" stroke-width="{sw}"/>')
                walk(c)

            # Support value label
            raw_name = node.get('name', '')
            disp, val = parse_support(raw_name)
            if disp is not None and val is not None:
                if support_type == "posterior":
                    hi_thresh, lo_thresh = 0.95, 0.70
                    hi_col, lo_col = "#fbbf24", "#475569"
                    show = val >= lo_thresh
                else:
                    hi_thresh, lo_thresh = 0.95, 0.70
                    hi_col, lo_col = "#fbbf24", "#475569"
                    show = val >= lo_thresh

                if show:
                    col = hi_col if val >= hi_thresh else lo_col
                    fs  = 7 if show_all else 9
                    labels.append(
                        f'<text x="{x1+2:.1f}" y="{y1-2:.1f}" font-size="{fs}" '
                        f'fill="{col}">{disp}</text>')

        else:
            name = node.get('name', '')
            col  = get_color(name)
            lbl  = name if node.get('_collapsed') else get_label(name)
            is_l = 'longirostris' in name.lower() or node.get('_collapsed')
            r    = 4 if node.get('_collapsed') else (1.5 if show_all else 3)
            fw   = 'bold' if is_l else 'normal'
            fs   = 10 if node.get('_collapsed') else (7 if show_all else 9)
            fi   = 'normal' if is_l else 'italic'
            dots.append(
                f'<circle cx="{x1+2:.1f}" cy="{y1:.1f}" r="{r}" fill="{col}"/>')
            labels.append(
                f'<text x="{x1+8:.1f}" y="{y1+3:.1f}" font-size="{fs}" '
                f'fill="{col}" font-weight="{fw}" font-style="{fi}">{lbl}</text>')

    walk(t)
    content = '\n'.join(lines + dots + labels)
    return (
        f'<svg id="{svg_id}" viewBox="0 0 {SVG_W} {SVG_H}" '
        f'style="width:100%;display:block">\n{content}\n</svg>',
        total
    )

# ── Load trees ────────────────────────────────────────────────────────────────
def load_tree(path, label):
    if not path or not os.path.exists(path):
        print(f"  WARNING: {label} tree not found: {path}")
        return None
    with open(path) as f:
        content = f.read().strip()

    # MrBayes .con.tre may contain a NEXUS block — extract the newick
    if content.upper().startswith("#NEXUS") or "BEGIN TREES" in content.upper():
        # Find the tree statement
        match = re.search(r'tree\s+\S+\s*=\s*(\[.*?\]\s*)?([(\[].+;)',
                          content, re.IGNORECASE | re.DOTALL)
        if match:
            newick = match.group(2).strip()
        else:
            # Try to find any line with a newick-like string
            for line in content.split('\n'):
                if '(' in line and ')' in line and ';' in line:
                    newick = line.strip()
                    break
            else:
                print(f"  WARNING: Could not extract newick from {label} tree file")
                return None
    else:
        newick = content

    try:
        tree = parse_newick(newick)
        n = count_leaves(tree)
        print(f"  {label}: {n} leaves")
        return tree
    except Exception as e:
        print(f"  ERROR parsing {label} tree: {e}")
        return None

print("Loading trees...")
mrbayes_tree = load_tree(args.mrbayes_tree, "MrBayes")
ml_tree      = load_tree(args.ml_tree,      "IQ-TREE ML")
nj_tree      = load_tree(args.nj_tree,      "NJ (CYTB)")

# ── Render SVGs ───────────────────────────────────────────────────────────────
tabs     = []
panels   = []
tab_idx  = 0

def add_tab(label, svg_html, n_leaves, notes, active=False):
    global tab_idx
    tid = f"panel-{tab_idx}"
    active_class = " active" if active else ""
    tabs.append(
        f'<div class="tab{active_class}" onclick="show(\'{tid}\',this)">{label}</div>')
    panels.append(
        f'<div id="{tid}" class="tree-panel{active_class}">'
        f'<p class="tree-note">{notes}</p>'
        f'{svg_html}'
        f'</div>'
    )
    tab_idx += 1

# -- MrBayes trees (clean + full)
if mrbayes_tree:
    print("Rendering MrBayes trees...")
    svg_mb_clean, n_mb = render_tree(
        mrbayes_tree, "mb-clean",
        use_bl=True, support_type="posterior",
        show_all=False, collapse_longirostris=True)
    svg_mb_full, _ = render_tree(
        mrbayes_tree, "mb-full",
        use_bl=True, support_type="posterior",
        show_all=True, collapse_longirostris=False)
    add_tab(
        "Bayesian (MrBayes) — collapsed", svg_mb_clean, n_mb,
        "Partitioned MrBayes tree (12S + ND4 + CYTB). Branch lengths proportional to "
        "substitutions. Node labels = posterior probability (yellow ≥0.95, grey ≥0.70). "
        "P. longirostris haplotypes collapsed.",
        active=True)
    add_tab(
        "Bayesian (MrBayes) — full", svg_mb_full, n_mb,
        "Full MrBayes tree showing all P. longirostris population representatives.")
else:
    tabs.append(
        '<div class="tab disabled">Bayesian (MrBayes) — not available</div>')

# -- IQ-TREE ML trees (clean + full)
if ml_tree:
    print("Rendering IQ-TREE ML trees...")
    svg_ml_clean, n_ml = render_tree(
        ml_tree, "ml-clean",
        use_bl=True, support_type="bootstrap",
        show_all=False, collapse_longirostris=True)
    svg_ml_full, _ = render_tree(
        ml_tree, "ml-full",
        use_bl=True, support_type="bootstrap",
        show_all=True, collapse_longirostris=False)
    first = not mrbayes_tree
    add_tab(
        "ML (IQ-TREE) — collapsed", svg_ml_clean, n_ml,
        "Maximum-likelihood tree (concatenated 12S + ND4 + CYTB, GTR+G). "
        "Node labels = ultrafast bootstrap (yellow ≥95, grey 70–94). "
        "P. longirostris collapsed.",
        active=first)
    add_tab(
        "ML (IQ-TREE) — full", svg_ml_full, n_ml,
        "Full IQ-TREE ML tree with all P. longirostris population representatives.")

# -- NJ tree (CYTB only)
if nj_tree:
    print("Rendering NJ tree...")
    svg_nj, n_nj = render_tree(
        nj_tree, "nj",
        use_bl=True, support_type="bootstrap",
        show_all=True, collapse_longirostris=True)
    first = not mrbayes_tree and not ml_tree
    add_tab(
        "NJ (CYTB)", svg_nj, n_nj,
        "Neighbour-joining tree based on CYTB alignment only (GTR+G distances via IQ-TREE). "
        "Included for comparison with Bayesian and ML trees.",
        active=first)

# ── Fetch summary table ───────────────────────────────────────────────────────
fetch_tsv = os.path.join(args.outdir, "fetch_summary.tsv")
fetch_table_html = ""
if os.path.exists(fetch_tsv):
    import csv
    with open(fetch_tsv) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    # Pivot: taxon × gene
    from collections import defaultdict
    by_taxon = defaultdict(dict)
    for r in rows:
        by_taxon[r['taxon']][r['gene']] = int(r.get('n_downloaded', 0))

    header = '<tr><th>Taxon</th><th>Role</th><th>12S</th><th>ND4</th><th>CYTB</th></tr>'
    body   = ""
    for taxon, genes in sorted(by_taxon.items()):
        role_val = next((r['role'] for r in rows if r['taxon'] == taxon), '')
        role_col  = "#f43f5e" if role_val == "outgroup" else "#94a3b8"
        cells = ""
        for g in ("12S","ND4","CYTB"):
            n = genes.get(g, 0)
            col = "#22c55e" if n > 0 else "#475569"
            cells += f'<td style="color:{col};text-align:center">{n if n > 0 else "—"}</td>'
        body += (f'<tr><td style="font-style:italic;color:#e2e8f0">{taxon}</td>'
                 f'<td style="color:{role_col};font-size:0.85em">{role_val}</td>'
                 f'{cells}</tr>')

    fetch_table_html = f"""
<div class="card" style="margin-top:20px">
  <h3>GenBank Coverage — Sequences Downloaded</h3>
  <div class="table-wrap">
    <table><thead>{header}</thead><tbody>{body}</tbody></table>
  </div>
</div>"""

# ── Model summary ─────────────────────────────────────────────────────────────
model_html = ""
for gene in ("12S","ND4","CYTB"):
    mt_file = os.path.join(args.outdir, f"modeltest_{gene}.out")
    if os.path.exists(mt_file):
        with open(mt_file) as f:
            mt_content = f.read()
        # Try to extract the BIC best model line
        match = re.search(r'BIC\s+\|\s+(\S+)', mt_content)
        if match:
            model_html += f"<li><b>{gene}:</b> {match.group(1)} (BIC)</li>"

if model_html:
    model_html = f'<h3>Best-Fit Substitution Models (modeltest-ng, BIC)</h3><ul>{model_html}</ul>'

# ── Legend ────────────────────────────────────────────────────────────────────
legend_items = [
    ('#f97316', 'P. longirostris (your data)'),
    ('#a78bfa', 'P. reynoldsi'),
    ('#34d399', 'North American Plestiodon'),
    ('#60a5fa', 'East Asian (barbouri group)'),
    ('#38bdf8', 'East Asian (elegans group)'),
    ('#4ade80', 'Mexican / Central American'),
    ('#f43f5e', 'Outgroups'),
    ('#94a3b8', 'Unclassified'),
]
legend_html = "".join(
    f'<span style="display:inline-flex;align-items:center;gap:5px;margin:3px 10px 3px 0">'
    f'<span style="width:10px;height:10px;border-radius:50%;background:{c};flex-shrink:0"></span>'
    f'<span style="color:#94a3b8;font-size:11px">{l}</span></span>'
    for c, l in legend_items
)

support_note = (
    "<b>MrBayes:</b> posterior probability shown at nodes "
    "(<span style='color:#fbbf24'>yellow ≥0.95</span>, "
    "<span style='color:#475569'>grey ≥0.70</span>, hidden &lt;0.70). &nbsp;|&nbsp; "
    "<b>IQ-TREE / NJ:</b> ultrafast bootstrap "
    "(<span style='color:#fbbf24'>yellow ≥95</span>, "
    "<span style='color:#475569'>grey 70–94</span>, hidden &lt;70)."
)

n_trees = sum(1 for t in [mrbayes_tree, ml_tree, nj_tree] if t is not None)
tab_bar  = "\n  ".join(tabs)
panel_html = "\n".join(panels)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Multigene Phylogeny — P. longirostris</title>
<style>
  * {{ box-sizing:border-box; margin:0; padding:0 }}
  body {{ background:#0f172a; color:#e2e8f0;
          font-family:'Courier New',monospace; padding:20px }}
  h1   {{ color:#f97316; font-size:20px; letter-spacing:1px; margin-bottom:4px }}
  h2   {{ color:#e2e8f0; font-size:15px; margin:16px 0 8px }}
  h3   {{ color:#94a3b8; font-size:13px; margin:12px 0 6px }}
  .subtitle {{ color:#475569; font-size:12px; margin-bottom:14px }}
  .legend   {{ margin-bottom:10px; line-height:2 }}
  .bs-note  {{ font-size:11px; color:#94a3b8; margin-bottom:12px }}
  .tab-bar  {{ display:flex; gap:6px; margin-bottom:10px; flex-wrap:wrap }}
  .tab {{ padding:5px 14px; border-radius:4px; cursor:pointer; font-size:11px;
           border:1px solid #334155; background:#1e293b; color:#94a3b8 }}
  .tab.active  {{ background:#f97316; border-color:#f97316;
                  color:#000; font-weight:bold }}
  .tab.disabled{{ opacity:0.4; cursor:default }}
  .tree-wrap {{ background:#0a1628; border:1px solid #1e293b;
                border-radius:8px; overflow:auto; max-height:85vh; padding:8px }}
  .tree-panel {{ display:none }}
  .tree-panel.active {{ display:block }}
  .tree-note  {{ font-size:11px; color:#475569; margin-bottom:6px }}
  .card {{ background:#0a1628; border:1px solid #1e293b; border-radius:8px;
            padding:16px; margin-top:16px }}
  .table-wrap {{ overflow-x:auto; margin-top:8px }}
  table {{ border-collapse:collapse; font-size:11px; width:100% }}
  th    {{ color:#64748b; padding:4px 8px; border-bottom:1px solid #1e293b;
           text-align:left }}
  td    {{ padding:3px 8px; border-bottom:1px solid #0f172a }}
  ul    {{ padding-left:20px; color:#94a3b8; font-size:12px; line-height:1.8 }}
  li b  {{ color:#e2e8f0 }}
  .note {{ font-size:11px; color:#475569; margin-top:8px }}
</style>
</head>
<body>
<h1>Multigene Phylogeny — Plestiodon longirostris</h1>
<p class="subtitle">Partitioned 12S + ND4 + CYTB ·
   MrBayes (Bayesian) · IQ-TREE (ML) · NJ (CYTB) ·
   {n_trees} tree(s) loaded</p>

<div class="legend">{legend_html}</div>
<p class="bs-note">{support_note}</p>

<div class="tab-bar">
  {tab_bar}
</div>
<div class="tree-wrap">
  {panel_html}
</div>
<p class="note">Scroll within the panel to navigate the tree.</p>

{model_html and f'<div class="card">{model_html}</div>' or ''}

{fetch_table_html}

<script>
function show(which, el) {{
  document.querySelectorAll('.tree-panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  document.getElementById(which).classList.add('active');
  el.classList.add('active');
}}
</script>
</body>
</html>"""

with open(args.output, "w") as f:
    f.write(html)

print(f"\nHTML report written → {args.output}")
print(f"  Trees rendered: {n_trees}")
print(f"  Tabs: {len(tabs)}")