#!/usr/bin/env python3
# =============================================================================
# fetch_multigene.py — Download 12S, ND4, and CYTB sequences from GenBank
#                      for Plestiodon multigene phylogeny
# =============================================================================
# Usage:
#   python3 fetch_multigene.py --email your@email.ac.uk --outdir multigene/
#
# Replaces fetch_cytb.py. Fetches three mitochondrial genes independently,
# writes per-gene combined FASTAs and a summary TSV.
#
# Requires: Biopython (pip install biopython)
# =============================================================================

import argparse
import os
import sys
import time
from Bio import Entrez, SeqIO

# ── Target taxa ───────────────────────────────────────────────────────────────
INGROUP_TAXA = [
    # North American
    "Plestiodon anthracinus", "Plestiodon egregius", "Plestiodon fasciatus",
    "Plestiodon gilberti", "Plestiodon inexpectatus", "Plestiodon laticeps",
    "Plestiodon longirostris", "Plestiodon multivirgatus", "Plestiodon obsoletus",
    "Plestiodon reynoldsi", "Plestiodon septentrionalis", "Plestiodon skiltonianus",
    "Plestiodon tetragrammus",
    # Mexican / Central American
    "Plestiodon brevirostris", "Plestiodon callicephalus", "Plestiodon capito",
    "Plestiodon colimensis", "Plestiodon copei", "Plestiodon dicei",
    "Plestiodon dugesii", "Plestiodon indubitus", "Plestiodon lagunensis",
    "Plestiodon longiartus", "Plestiodon lotus", "Plestiodon lynxe",
    "Plestiodon multilineatus", "Plestiodon nietoi", "Plestiodon ochoteranae",
    "Plestiodon parviauriculatus", "Plestiodon parvulus", "Plestiodon sumichrasti",
    # East Asian
    "Plestiodon barbouri", "Plestiodon chinensis", "Plestiodon coreensis",
    "Plestiodon elegans", "Plestiodon finitimus", "Plestiodon japonicus",
    "Plestiodon kishinouyei", "Plestiodon kuchinoshimensis", "Plestiodon latiscutatus",
    "Plestiodon leucostictus", "Plestiodon liui", "Plestiodon marginatus",
    "Plestiodon oshimensis", "Plestiodon popei", "Plestiodon quadrilineatus",
    "Plestiodon stimpsonii", "Plestiodon takarai", "Plestiodon tamdaoensis",
    "Plestiodon tunganus", "Plestiodon bilineatus",
]

OUTGROUP_TAXA = [
    "Eumeces schneiderii",
    "Scincella lateralis",
]

# Old Eumeces names for fallback searching
EUMECES_FALLBACK = {
    "Plestiodon anthracinus":      "Eumeces anthracinus",
    "Plestiodon callicephalus":    "Eumeces callicephalus",
    "Plestiodon brevirostris":     "Eumeces brevirostris",
    "Plestiodon capito":           "Eumeces capito",
    "Plestiodon colimensis":       "Eumeces colimensis",
    "Plestiodon copei":            "Eumeces copei",
    "Plestiodon dicei":            "Eumeces dicei",
    "Plestiodon dugesii":          "Eumeces dugesii",
    "Plestiodon fasciatus":        "Eumeces fasciatus",
    "Plestiodon gilberti":         "Eumeces gilberti",
    "Plestiodon indubitus":        "Eumeces indubitus",
    "Plestiodon inexpectatus":     "Eumeces inexpectatus",
    "Plestiodon lagunensis":       "Eumeces lagunensis",
    "Plestiodon laticeps":         "Eumeces laticeps",
    "Plestiodon longirostris":     "Eumeces longirostris",
    "Plestiodon lynxe":            "Eumeces lynxe",
    "Plestiodon multilineatus":    "Eumeces multilineatus",
    "Plestiodon multivirgatus":    "Eumeces multivirgatus",
    "Plestiodon nietoi":           "Eumeces nietoi",
    "Plestiodon obsoletus":        "Eumeces obsoletus",
    "Plestiodon ochoteranae":      "Eumeces ochoteranae",
    "Plestiodon parviauriculatus": "Eumeces parviauriculatus",
    "Plestiodon parvulus":         "Eumeces parvulus",
    "Plestiodon reynoldsi":        "Eumeces reynoldsi",
    "Plestiodon septentrionalis":  "Eumeces septentrionalis",
    "Plestiodon skiltonianus":     "Eumeces skiltonianus",
    "Plestiodon sumichrasti":      "Eumeces sumichrasti",
    "Plestiodon tetragrammus":     "Eumeces tetragrammus",
}

# ── Gene definitions ──────────────────────────────────────────────────────────
# Each gene has:
#   name        : display name
#   gene_names  : list of synonyms to match in GenBank gene/product qualifiers
#   title_terms : terms to search in sequence titles
#   min_len     : minimum acceptable extracted length (bp)
#   max_len     : maximum acceptable extracted length (bp)
#   is_rrna     : True for ribosomal RNA genes (use rRNA feature type)

GENES = {
    "12S": {
        "name":       "12S rRNA",
        "gene_names": ["12s", "12s rrna", "12s ribosomal rna", "s-rrna", "rrns"],
        "title_terms":["12s", "12s rrna", "small subunit ribosomal"],
        "min_len":    700,
        "max_len":    1100,
        "is_rrna":    True,
    },
    "ND4": {
        "name":       "ND4",
        "gene_names": ["nd4", "nadh4", "nadh dehydrogenase subunit 4",
                       "nadh-ubiquinone oxidoreductase chain 4"],
        "title_terms":["nd4", "nadh4", "nadh dehydrogenase subunit 4"],
        "min_len":    1200,
        "max_len":    1400,
        "is_rrna":    False,
    },
    "CYTB": {
        "name":       "cytochrome b",
        "gene_names": ["cytb", "cyt b", "cytochrome b", "cytb gene"],
        "title_terms":["cytochrome b", "cytb", "cyt b"],
        "min_len":    900,
        "max_len":    1200,
        "is_rrna":    False,
    },
}

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Fetch 12S, ND4, and CYTB sequences from GenBank for Plestiodon phylogeny")
parser.add_argument("--email",          required=True,
                    help="Email address for NCBI Entrez")
parser.add_argument("--outdir",         default="multigene",
                    help="Output directory (default: multigene/)")
parser.add_argument("--max-per-species",type=int, default=3,
                    help="Max sequences per species per gene (default: 3)")
parser.add_argument("--outgroups-only", action="store_true",
                    help="Only download outgroup sequences")
parser.add_argument("--genes",          default="12S,ND4,CYTB",
                    help="Comma-separated list of genes to fetch (default: 12S,ND4,CYTB)")
parser.add_argument("--delay",          type=float, default=0.4,
                    help="Delay between NCBI requests in seconds (default: 0.4)")
args = parser.parse_args()

Entrez.email = args.email
os.makedirs(args.outdir, exist_ok=True)

GENES_TO_FETCH = [g.strip().upper() for g in args.genes.split(",")]
for g in GENES_TO_FETCH:
    if g not in GENES:
        print(f"ERROR: unknown gene '{g}'. Valid options: {', '.join(GENES.keys())}")
        sys.exit(1)

# ── Helpers ───────────────────────────────────────────────────────────────────
def search_taxon_gene(taxon, gene_key, max_per):
    """
    Three-strategy search for a specific gene in a specific taxon.
    Returns list of GenBank IDs (deduplicated).
    """
    gdef = GENES[gene_key]
    min_len = gdef["min_len"]
    max_len = gdef["max_len"]

    # Build title search terms
    title_clause = " OR ".join(f'{t}[Title]' for t in gdef["title_terms"])

    queries = [
        # Strategy 1: standalone gene record by title
        (f'"{taxon}"[Organism] AND ({title_clause}) AND {min_len}:{max_len}[SLEN]'),
        # Strategy 2: gene feature annotation search
        (f'"{taxon}"[Organism] AND {gdef["gene_names"][0]}[Gene] AND {min_len}:{max_len}[SLEN]'),
        # Strategy 3: complete mitogenome (extract feature later)
        (f'"{taxon}"[Organism] AND mitochondri*[Title] AND complete[Title] AND 10000:20000[SLEN]'),
    ]

    all_ids = []
    seen = set()
    for query in queries:
        try:
            handle = Entrez.esearch(db="nuccore", term=query, retmax=200)
            record = Entrez.read(handle)
            handle.close()
            for i in record["IdList"]:
                if i not in seen:
                    seen.add(i)
                    all_ids.append(i)
            time.sleep(0.2)
        except Exception as e:
            print(f"  Search error ({taxon}, {gene_key}): {e}")
    return all_ids


def extract_gene_from_record(rec, gene_key):
    """
    Extract a specific gene feature from a GenBank record.
    For standalone short records: return as-is if length within range.
    For mitogenomes: find and extract the CDS/rRNA feature.
    Returns a SeqRecord or None.
    """
    gdef   = GENES[gene_key]
    min_len = gdef["min_len"]
    max_len = gdef["max_len"]

    # Already a standalone gene-length record
    if min_len <= len(rec.seq) <= max_len:
        return rec

    # Too short even for a gene — skip
    if len(rec.seq) < min_len:
        return None

    # Mitogenome — search features
    feature_types = ["rRNA"] if gdef["is_rrna"] else ["CDS", "gene"]
    gene_synonyms = gdef["gene_names"]

    for feat in rec.features:
        if feat.type not in feature_types:
            continue
        qualifiers = feat.qualifiers
        gene_val    = " ".join(qualifiers.get("gene",    [""])).lower()
        product_val = " ".join(qualifiers.get("product", [""])).lower()
        combined    = gene_val + " " + product_val

        if any(syn in combined for syn in gene_synonyms):
            try:
                sub = feat.extract(rec)
                sub.id          = rec.id
                sub.description = rec.description
                if min_len <= len(sub.seq) <= max_len:
                    return sub
            except Exception:
                pass
    return None


def fetch_records(ids):
    """Fetch GenBank records for a list of IDs in one call."""
    try:
        handle  = Entrez.efetch(db="nuccore", id=",".join(ids),
                                rettype="gb", retmode="text")
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        return records
    except Exception as e:
        print(f"  Fetch error: {e}")
        return []


def clean_name(organism):
    return organism.replace(" ", "_").replace("/", "_")


# ── Main loop ─────────────────────────────────────────────────────────────────
taxa = OUTGROUP_TAXA if args.outgroups_only else INGROUP_TAXA + OUTGROUP_TAXA

# gene_key -> list of SeqRecord
all_sequences = {g: [] for g in GENES_TO_FETCH}
# gene_key -> list of summary dicts
summary_rows  = {g: [] for g in GENES_TO_FETCH}

print(f"{'='*70}")
print(f"  Multigene CYTB/ND4/12S fetch — {len(taxa)} taxa × {len(GENES_TO_FETCH)} genes")
print(f"  Max per species per gene: {args.max_per_species}")
print(f"  Output: {args.outdir}/")
print(f"{'='*70}\n")

for taxon in taxa:
    role = "outgroup" if taxon in OUTGROUP_TAXA else "ingroup"
    print(f"\n{'─'*60}")
    print(f"  [{role}] {taxon}")

    for gene_key in GENES_TO_FETCH:
        gname = GENES[gene_key]["name"]
        print(f"    {gene_key:6} ...", end=" ", flush=True)

        ids       = search_taxon_gene(taxon, gene_key, args.max_per_species)
        used_name = taxon
        time.sleep(args.delay)

        # Fallback to old Eumeces name
        if not ids and taxon in EUMECES_FALLBACK:
            fallback = EUMECES_FALLBACK[taxon]
            ids      = search_taxon_gene(fallback, gene_key, args.max_per_species)
            time.sleep(args.delay)
            if ids:
                used_name = fallback
                print(f"(retried as {fallback}) ", end="", flush=True)

        if not ids:
            print("not found")
            summary_rows[gene_key].append({
                "taxon": taxon, "role": role, "searched_as": used_name,
                "n_found": 0, "n_downloaded": 0, "accessions": ""
            })
            continue

        # Fetch more records than needed to account for mitogenome extraction failures
        fetch_ids = ids[:args.max_per_species * 5]
        records   = fetch_records(fetch_ids)
        time.sleep(args.delay)

        if not records:
            print("fetch failed")
            summary_rows[gene_key].append({
                "taxon": taxon, "role": role, "searched_as": used_name,
                "n_found": len(ids), "n_downloaded": 0, "accessions": ""
            })
            continue

        accessions = []
        kept       = 0
        for rec in records:
            if kept >= args.max_per_species:
                break
            processed = extract_gene_from_record(rec, gene_key)
            if processed is None:
                continue

            safe_name         = clean_name(taxon)   # always modern Plestiodon name
            orig_id           = rec.id
            processed.id      = f"{safe_name}_{orig_id}"
            processed.name    = processed.id
            processed.description = ""

            all_sequences[gene_key].append(processed)
            accessions.append(processed.id)

            # Per-species per-gene FASTA
            per_sp_fa = os.path.join(args.outdir, f"{clean_name(taxon)}_{gene_key}.fa")
            with open(per_sp_fa, "a") as fh:
                SeqIO.write(processed, fh, "fasta")
            kept += 1

        print(f"{kept} seq  ({', '.join(a.split('_')[-1] for a in accessions)})")
        summary_rows[gene_key].append({
            "taxon": taxon, "role": role, "searched_as": used_name,
            "n_found": len(ids), "n_downloaded": kept,
            "accessions": "; ".join(accessions)
        })

# ── Write per-gene combined FASTAs ────────────────────────────────────────────
print(f"\n{'='*70}")
print("  Writing combined FASTAs")
print(f"{'='*70}")

for gene_key in GENES_TO_FETCH:
    out_fa = os.path.join(args.outdir, f"all_{gene_key}_genbank.fa")
    with open(out_fa, "w") as fh:
        SeqIO.write(all_sequences[gene_key], fh, "fasta")
    n = len(all_sequences[gene_key])
    print(f"  {gene_key:6}: {n:3} sequences → {out_fa}")

# ── Write summary TSV ─────────────────────────────────────────────────────────
summary_tsv = os.path.join(args.outdir, "fetch_summary.tsv")
with open(summary_tsv, "w") as fh:
    fh.write("gene\ttaxon\trole\tsearched_as\tn_found\tn_downloaded\taccessions\n")
    for gene_key in GENES_TO_FETCH:
        for row in summary_rows[gene_key]:
            fh.write(
                f"{gene_key}\t{row['taxon']}\t{row['role']}\t"
                f"{row.get('searched_as', row['taxon'])}\t"
                f"{row['n_found']}\t{row['n_downloaded']}\t{row['accessions']}\n"
            )
print(f"\n  Summary TSV → {summary_tsv}")

# ── Coverage report ───────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print("  Coverage summary")
print(f"{'='*70}")
for gene_key in GENES_TO_FETCH:
    rows    = summary_rows[gene_key]
    found   = sum(1 for r in rows if r["n_downloaded"] > 0)
    missing = [r["taxon"] for r in rows if r["n_downloaded"] == 0]
    print(f"\n  {gene_key} ({GENES[gene_key]['name']}): {found}/{len(rows)} taxa")
    if missing:
        for m in missing:
            print(f"    ✗  {m}")