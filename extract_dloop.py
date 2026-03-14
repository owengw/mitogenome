#!/usr/bin/env python3
# =============================================================================
# extract_dloop.py — Extract D-loop / control region from MitoZ GFF annotation
#
# For each sample:
#   1. Find the MitoZ GFF file in the annotation directory
#   2. Locate the D-loop / control region feature
#   3. Extract the sequence from the assembled mitogenome FASTA
#   4. Label with sample ID and population code
#   5. Write all sequences to a combined FASTA
#
# Falls back to positional extraction (last ~900bp of mitogenome) if no
# control region feature is annotated.
# =============================================================================

import argparse
import os
import re
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# D-loop feature names as annotated by MitoZ
DLOOP_FEATURE_NAMES = [
    "control region", "d-loop", "dloop", "d loop",
    "CR", "control_region", "non-coding region",
    "rep_origin",
]

parser = argparse.ArgumentParser(description="Extract D-loop from MitoZ GFF files")
parser.add_argument("--annot-dir",    required=True,
                    help="MitoZ annotation output directory (contains per-sample subdirs)")
parser.add_argument("--assembly-dir", required=True,
                    help="Directory with assembled mitogenome FASTAs (one per sample)")
parser.add_argument("--outfile",      required=True,
                    help="Output FASTA file for all D-loop sequences")
parser.add_argument("--pop-defs",     nargs="+", default=[],
                    help="Population definitions: PATTERN:POP_CODE:HEX_COLOUR ...")
parser.add_argument("--min-len",      type=int, default=400,
                    help="Minimum D-loop length to accept (default: 400)")
parser.add_argument("--max-len",      type=int, default=1500,
                    help="Maximum D-loop length to accept (default: 1500)")
args = parser.parse_args()

# Parse population definitions
pop_defs = []
for pd in args.pop_defs:
    parts = pd.split(":")
    if len(parts) == 3:
        pop_defs.append({"pattern": parts[0], "pop": parts[1], "colour": parts[2]})
    elif len(parts) == 2:
        pop_defs.append({"pattern": parts[0], "pop": parts[1], "colour": "#94a3b8"})

def get_population(sample_id):
    """Match sample ID to population from pop_defs."""
    for pd in pop_defs:
        if pd["pattern"].upper() in sample_id.upper():
            return pd["pop"], pd["colour"]
    return "UNK", "#94a3b8"

def find_gff_files(annot_dir):
    """Find all GFF files in MitoZ annotation subdirectories."""
    gff_files = []
    annot_path = Path(annot_dir)
    # MitoZ writes: annot_dir/<sample>/<prefix>.result/<prefix>.genes.gff or similar
    for gff in annot_path.rglob("*.gff"):
        gff_files.append(gff)
    for gff in annot_path.rglob("*.gtf"):
        gff_files.append(gff)
    return gff_files

def find_fasta_for_sample(assembly_dir, sample_id):
    """Find assembled mitogenome FASTA for a given sample ID."""
    asm_path = Path(assembly_dir)
    # Try exact name match first
    candidates = list(asm_path.glob(f"*{sample_id}*.fa")) + \
                 list(asm_path.glob(f"*{sample_id}*.fasta")) + \
                 list(asm_path.glob(f"*{sample_id}*.fas"))
    if candidates:
        return candidates[0]
    return None

def parse_gff_dloop(gff_path):
    """
    Parse GFF/GTF file and return (start, end, strand) of D-loop feature.
    Returns None if not found.
    Coordinates are 1-based inclusive (GFF3 standard).
    """
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2].lower()
            attributes   = parts[8].lower()
            start        = int(parts[3])
            end          = int(parts[4])
            strand       = parts[6]

            # Check feature type and attributes for D-loop keywords
            is_dloop = False
            for kw in DLOOP_FEATURE_NAMES:
                if kw.lower() in feature_type or kw.lower() in attributes:
                    is_dloop = True
                    break

            if is_dloop:
                return (start - 1, end, strand)  # convert to 0-based for Python slicing
    return None

def extract_dloop_from_fasta(fasta_path, start, end, strand):
    """Extract subsequence from FASTA given 0-based start, 1-based end."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None
    # Take first record (assembled mitogenome)
    rec = records[0]
    sub = rec.seq[start:end]
    if strand == "-":
        sub = sub.reverse_complement()
    return sub

def positional_fallback(fasta_path, approx_len=900):
    """
    Fallback: extract last ~900bp of mitogenome as approximate D-loop region.
    The D-loop is typically at the 3' end of the conventional mitogenome map.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None
    seq = records[0].seq
    return seq[-approx_len:]

# ── Main extraction loop ──────────────────────────────────────────────────────
gff_files = find_gff_files(args.annot_dir)
print(f"Found {len(gff_files)} GFF files in {args.annot_dir}")

if not gff_files:
    print("WARNING: No GFF files found — checking assembly directory for FASTAs")

# Also scan assembly directory directly for FASTA files
asm_fastas = list(Path(args.assembly_dir).glob("*.fa")) + \
             list(Path(args.assembly_dir).glob("*.fasta"))
print(f"Found {len(asm_fastas)} FASTA files in {args.assembly_dir}")

extracted = []
used_fallback = []
failed = []

# Build a mapping from sample ID to GFF path
# MitoZ subdirectory names typically match sample IDs
sample_gff_map = {}
for gff in gff_files:
    # Extract sample name from path: annot_dir/<sample>/.../<prefix>.gff
    parts = gff.parts
    annot_parts = Path(args.annot_dir).parts
    # Sample dir is immediately below annot_dir
    if len(parts) > len(annot_parts):
        sample_dir = parts[len(annot_parts)]
        if sample_dir not in sample_gff_map:
            sample_gff_map[sample_dir] = gff

print(f"\nProcessing {len(asm_fastas)} samples...")
print(f"{'─'*60}")

for fasta_path in sorted(asm_fastas):
    sample_id = fasta_path.stem
    # Remove common suffixes
    for suffix in ["_mito", "_consensus", "_assembly", "_mitogenome"]:
        sample_id = sample_id.replace(suffix, "")

    pop, colour = get_population(sample_id)

    # Find matching GFF
    gff_path = None
    for sid, gff in sample_gff_map.items():
        if sample_id.upper() in sid.upper() or sid.upper() in sample_id.upper():
            gff_path = gff
            break

    dloop_seq = None
    method    = "unknown"

    if gff_path and gff_path.exists():
        coords = parse_gff_dloop(gff_path)
        if coords:
            start, end, strand = coords
            dloop_seq = extract_dloop_from_fasta(fasta_path, start, end, strand)
            method = f"GFF ({start+1}-{end}, {strand})"

    if dloop_seq is None or len(dloop_seq) < args.min_len:
        # Fallback: positional extraction
        dloop_seq = positional_fallback(fasta_path)
        method    = "positional_fallback"
        used_fallback.append(sample_id)

    if dloop_seq is None or len(dloop_seq) < args.min_len:
        print(f"  ✗  {sample_id:<20} — extraction failed")
        failed.append(sample_id)
        continue

    # Truncate if too long
    if len(dloop_seq) > args.max_len:
        dloop_seq = dloop_seq[:args.max_len]

    # Build record with population info in header
    rec_id = f"{sample_id}|{pop}"
    rec = SeqRecord(
        Seq(str(dloop_seq)),
        id=rec_id,
        name=rec_id,
        description=f"pop={pop} colour={colour} method={method}"
    )
    extracted.append(rec)
    print(f"  ✔  {sample_id:<20} pop={pop:<5} len={len(dloop_seq)}  [{method}]")

# Write output FASTA
with open(args.outfile, "w") as fh:
    SeqIO.write(extracted, fh, "fasta")

print(f"\n{'='*60}")
print(f"  Extracted:      {len(extracted)} sequences")
print(f"  GFF-based:      {len(extracted) - len(used_fallback)}")
print(f"  Positional fb:  {len(used_fallback)}")
print(f"  Failed:         {len(failed)}")
if used_fallback:
    print(f"\n  Fallback samples (no GFF annotation):")
    for s in used_fallback:
        print(f"    {s}")
if failed:
    print(f"\n  Failed samples:")
    for s in failed:
        print(f"    {s}")
print(f"\n  Output → {args.outfile}")