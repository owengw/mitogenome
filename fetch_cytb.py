#!/usr/bin/env python3
# =============================================================================
# fetch_cytb.py — Download CYTB sequences from GenBank for Plestiodon phylogeny
# =============================================================================
# Usage:
#   python3 fetch_cytb.py --email your@email.ac.uk --outdir cytb_sequences/
#
# Requires: Biopython (pip install biopython)
# =============================================================================

import argparse
import os
import sys
import time
from Bio import Entrez, SeqIO

# ── Target taxa ───────────────────────────────────────────────────────────────
# Ingroup: all Plestiodon species from Brandley et al. + any additional
# Outgroups: closest outgroups used in Brandley et al.
INGROUP_TAXA = [
    # North American species
    "Plestiodon anthracinus",
    "Plestiodon egregius",
    "Plestiodon fasciatus",
    "Plestiodon gilberti",
    "Plestiodon inexpectatus",
    "Plestiodon laticeps",
    "Plestiodon longirostris",
    "Plestiodon multivirgatus",
    "Plestiodon obsoletus",
    "Plestiodon reynoldsi",
    "Plestiodon septentrionalis",
    "Plestiodon skiltonianus",
    "Plestiodon tetragrammus",
    # Mexican / Central American species
    "Plestiodon brevirostris",
    "Plestiodon callicephalus",
    "Plestiodon capito",
    "Plestiodon colimensis",
    "Plestiodon copei",
    "Plestiodon dicei",
    "Plestiodon dugesii",
    "Plestiodon indubitus",
    "Plestiodon lagunensis",
    "Plestiodon longiartus",
    "Plestiodon lotus",
    "Plestiodon lynxe",
    "Plestiodon multilineatus",
    "Plestiodon nietoi",
    "Plestiodon ochoteranae",
    "Plestiodon parviauriculatus",
    "Plestiodon parvulus",
    "Plestiodon sumichrasti",
    # East Asian species
    "Plestiodon barbouri",
    "Plestiodon chinensis",
    "Plestiodon coreensis",
    "Plestiodon elegans",
    "Plestiodon finitimus",
    "Plestiodon japonicus",
    "Plestiodon kishinouyei",
    "Plestiodon kuchinoshimensis",
    "Plestiodon latiscutatus",
    "Plestiodon leucostictus",
    "Plestiodon liui",
    "Plestiodon marginatus",
    "Plestiodon oshimensis",
    "Plestiodon popei",
    "Plestiodon quadrilineatus",
    "Plestiodon stimpsonii",
    "Plestiodon takarai",
    "Plestiodon tamdaoensis",
    "Plestiodon tunganus",
    # Recently described species (post-2010, not in Brandley et al.)
    "Plestiodon bilineatus",
]

OUTGROUP_TAXA = [
    # Two closest outgroups only — more distant ones cause long-branch attraction
    "Eumeces schneiderii",
    "Scincella lateralis",
    "Trachylepis quinquetaeniata",
]

# ── Eumeces fallback mapping ──────────────────────────────────────────────────
# Many North/Central American Plestiodon were formerly Eumeces.
# If a Plestiodon search returns nothing, retry with the old name.
EUMECES_FALLBACK = {
    "Plestiodon anthracinus":    "Eumeces anthracinus",
    "Plestiodon callicephalus":  "Eumeces callicephalus",
    "Plestiodon brevirostris":   "Eumeces brevirostris",
    "Plestiodon capito":         "Eumeces capito",
    "Plestiodon colimensis":     "Eumeces colimensis",
    "Plestiodon copei":          "Eumeces copei",
    "Plestiodon dicei":          "Eumeces dicei",
    "Plestiodon dugesii":        "Eumeces dugesii",
    "Plestiodon fasciatus":      "Eumeces fasciatus",
    "Plestiodon gilberti":       "Eumeces gilberti",
    "Plestiodon indubitus":      "Eumeces indubitus",
    "Plestiodon inexpectatus":   "Eumeces inexpectatus",
    "Plestiodon lagunensis":     "Eumeces lagunensis",
    "Plestiodon laticeps":       "Eumeces laticeps",
    "Plestiodon longirostris":   "Eumeces longirostris",
    "Plestiodon lynxe":          "Eumeces lynxe",
    "Plestiodon multilineatus":  "Eumeces multilineatus",
    "Plestiodon multivirgatus":  "Eumeces multivirgatus",
    "Plestiodon nietoi":         "Eumeces nietoi",
    "Plestiodon obsoletus":      "Eumeces obsoletus",
    "Plestiodon ochoteranae":    "Eumeces ochoteranae",
    "Plestiodon parviauriculatus": "Eumeces parviauriculatus",
    "Plestiodon parvulus":       "Eumeces parvulus",
    "Plestiodon reynoldsi":      "Eumeces reynoldsi",
    "Plestiodon septentrionalis":"Eumeces septentrionalis",
    "Plestiodon skiltonianus":   "Eumeces skiltonianus",
    "Plestiodon sumichrasti":    "Eumeces sumichrasti",
    "Plestiodon tetragrammus":   "Eumeces tetragrammus",
}

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Fetch CYTB sequences from GenBank for Plestiodon phylogeny")
parser.add_argument("--email",  required=True,
                    help="Email address for NCBI Entrez (required by NCBI)")
parser.add_argument("--outdir", default="cytb_sequences",
                    help="Output directory (default: cytb_sequences/)")
parser.add_argument("--min-len", type=int, default=500,
                    help="Minimum sequence length in bp (default: 500)")
parser.add_argument("--max-len", type=int, default=1200,
                    help="Maximum sequence length in bp (default: 1200)")
parser.add_argument("--max-per-species", type=int, default=3,
                    help="Max sequences to download per species (default: 3)")
parser.add_argument("--outgroups-only", action="store_true",
                    help="Only download outgroup sequences")
parser.add_argument("--delay", type=float, default=0.4,
                    help="Delay between NCBI requests in seconds (default: 0.4)")
args = parser.parse_args()

Entrez.email = args.email
os.makedirs(args.outdir, exist_ok=True)

# ── Helpers ───────────────────────────────────────────────────────────────────
def search_taxon(taxon, min_len, max_len):
    """
    Search NCBI for CYTB sequences. Uses two strategies:
    1. Standalone CYTB sequences (title-based, fast)
    2. Mitochondrial genome records (gene feature-based, catches older deposits)
    Returns combined deduplicated ID list.
    """
    queries = [
        # Strategy 1: standalone CYTB gene record
        (f'"{taxon}"[Organism] AND '
         f'(cytochrome b[Title] OR cytb[Title] OR "cyt b"[Title]) AND '
         f'{min_len}:{max_len}[SLEN]'),
        # Strategy 2: broader — any mitochondrial record with cytb annotation
        (f'"{taxon}"[Organism] AND '
         f'cytochrome b[Gene] AND '
         f'{min_len}:{max_len}[SLEN]'),
        # Strategy 3: complete/partial mitochondrial genome (extract CYTB later)
        (f'"{taxon}"[Organism] AND '
         f'mitochondri*[Title] AND complete[Title] AND '
         f'10000:20000[SLEN]'),
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
            print(f"  Search error for {taxon}: {e}")
    return all_ids

def extract_cytb_from_record(rec, min_len, max_len):
    """
    If record is a mitogenome, extract just the CYTB feature as a new SeqRecord.
    If record is already a short standalone sequence, return as-is.
    """
    if len(rec.seq) <= max_len:
        return rec  # already a standalone CYTB-length record

    # Search for CYTB feature in mitogenome
    for feat in rec.features:
        if feat.type == "CDS":
            qualifiers = feat.qualifiers
            gene = qualifiers.get("gene", [""])[0].lower()
            product = qualifiers.get("product", [""])[0].lower()
            if "cytb" in gene or "cytochrome b" in product or "cyt b" in gene:
                try:
                    sub = feat.extract(rec)
                    sub.id = rec.id
                    sub.description = rec.description
                    if min_len <= len(sub.seq) <= max_len:
                        return sub
                except Exception:
                    pass
    return None  # no CYTB feature found or wrong length

def fetch_records(ids):
    """Fetch GenBank records for a list of IDs."""
    try:
        handle = Entrez.efetch(db="nuccore", id=",".join(ids),
                               rettype="gb", retmode="text")
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        return records
    except Exception as e:
        print(f"  Fetch error: {e}")
        return []

def clean_name(organism):
    """Convert organism name to safe filename component."""
    return organism.replace(" ", "_").replace("/", "_")

# ── Main loop ─────────────────────────────────────────────────────────────────
taxa = OUTGROUP_TAXA if args.outgroups_only else INGROUP_TAXA + OUTGROUP_TAXA

all_sequences = []   # collected SeqRecord objects for combined FASTA
summary_rows  = []   # for the TSV summary

print(f"{'='*70}")
print(f"  CYTB sequence fetch — {len(taxa)} taxa")
print(f"  Length filter: {args.min_len}–{args.max_len} bp")
print(f"  Max per species: {args.max_per_species}")
print(f"  Output: {args.outdir}/")
print(f"{'='*70}\n")

for taxon in taxa:
    role = "outgroup" if taxon in OUTGROUP_TAXA else "ingroup"
    print(f"  [{role}] {taxon} ...", end=" ", flush=True)

    ids = search_taxon(taxon, args.min_len, args.max_len)
    time.sleep(args.delay)
    used_name = taxon

    # Retry with old Eumeces name if no results
    if not ids and taxon in EUMECES_FALLBACK:
        fallback = EUMECES_FALLBACK[taxon]
        ids = search_taxon(fallback, args.min_len, args.max_len)
        time.sleep(args.delay)
        if ids:
            used_name = fallback
            print(f"(retried as {fallback}) ", end="", flush=True)

    if not ids:
        print("no sequences found")
        summary_rows.append({
            "taxon": taxon, "role": role, "searched_as": used_name,
            "n_found": 0, "n_downloaded": 0, "accessions": ""
        })
        continue

    # Fetch in batches — limit IDs but fetch more to account for mitogenomes
    # where we may need to extract the CYTB feature (not all will have it)
    ids = ids[:args.max_per_species * 5]   # fetch more, filter after extraction
    records = fetch_records(ids)
    time.sleep(args.delay)

    if not records:
        print("fetch failed")
        summary_rows.append({
            "taxon": taxon, "role": role, "searched_as": used_name,
            "n_found": len(ids), "n_downloaded": 0, "accessions": ""
        })
        continue

    accessions = []
    kept = 0
    for rec in records:
        if kept >= args.max_per_species:
            break
        # Extract CYTB feature if this is a mitogenome
        processed = extract_cytb_from_record(rec, args.min_len, args.max_len)
        if processed is None:
            continue   # mitogenome with no CYTB feature, or wrong length

        # Rename to canonical Plestiodon name + accession for clarity in tree
        safe_name = clean_name(taxon)   # always use modern Plestiodon name
        orig_id = rec.id
        processed.id          = f"{safe_name}_{orig_id}"
        processed.name        = processed.id
        processed.description = ""
        all_sequences.append(processed)
        accessions.append(processed.id)

        # Also write per-species FASTA
        per_sp_fa = os.path.join(args.outdir, f"{clean_name(taxon)}.fa")
        with open(per_sp_fa, "a") as fh:
            SeqIO.write(processed, fh, "fasta")
        kept += 1

    print(f"{kept} sequences  ({', '.join(accessions)})")
    summary_rows.append({
        "taxon": taxon, "role": role, "searched_as": used_name,
        "n_found": len(ids), "n_downloaded": kept,
        "accessions": "; ".join(accessions)
    })

# ── Write combined FASTA ──────────────────────────────────────────────────────
combined_fa = os.path.join(args.outdir, "all_cytb_outgroups.fa")
with open(combined_fa, "w") as fh:
    SeqIO.write(all_sequences, fh, "fasta")
print(f"\nCombined FASTA → {combined_fa}  ({len(all_sequences)} sequences)")

# ── Write summary TSV ─────────────────────────────────────────────────────────
summary_tsv = os.path.join(args.outdir, "fetch_summary.tsv")
with open(summary_tsv, "w") as fh:
    fh.write("taxon\trole\tsearched_as\tn_found\tn_downloaded\taccessions\n")
    for row in summary_rows:
        fh.write(f"{row['taxon']}\t{row['role']}\t{row.get('searched_as', row['taxon'])}\t"
                 f"{row['n_found']}\t{row['n_downloaded']}\t{row['accessions']}\n")
print(f"Summary TSV    → {summary_tsv}")

# ── Print coverage report ─────────────────────────────────────────────────────
found    = [r for r in summary_rows if r["n_downloaded"] > 0]
missing  = [r for r in summary_rows if r["n_downloaded"] == 0]

print(f"\n{'='*70}")
print(f"  Coverage: {len(found)}/{len(taxa)} taxa with sequences")
if missing:
    print(f"\n  Missing taxa (no CYTB found):")
    for r in missing:
        print(f"    {r['role']:10}  {r['taxon']}")
print(f"{'='*70}")
print("\nNext steps:")
print(f"  1. Add your P. longirostris consensus CYTB sequences to {combined_fa}")
print(f"  2. Align with MAFFT: mafft --auto all_cytb_combined.fa > all_cytb_aligned.fa")
print(f"  3. Run IQ-TREE:  iqtree -s all_cytb_aligned.fa -m GTR+G -B 1000 -T AUTO")