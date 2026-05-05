#!/usr/bin/env python3
"""
M. pneumoniae MLVA typer — Dégrange et al. 2009 amplicon-size method.

Per locus:
  1. Locate the F primer (forward strand) and R primer (as RC on forward
     strand) in the assembly, on either strand of any contig. Allow up to
     `primer_max_mismatches` Hamming-distance mismatches.
  2. Amplicon = sequence from F primer 5' end to R primer 3' end.
  3. copies = round((amplicon_length − constant) / motif_length)

The `constant` per locus is calibrated from M129 in the build step.

Usage:
  python mlva_type.py --assembly contigs.fasta --loci loci.yaml [--out profile.json]
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq


@dataclass
class LocusCall:
    name: str
    copies: int | None
    motif_length: int
    amplicon_length: int | None
    array_length: int | None
    contig: str | None
    strand: str | None
    confidence: str   # "high" | "low" | "missing"
    note: str = ""


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def all_hits(seq: str, query: str, max_mm: int) -> list[tuple[int, int]]:
    """Return all (position, hamming-distance) hits of query in seq within max_mm."""
    n, k = len(seq), len(query)
    return [(i, hamming(seq[i:i + k], query)) for i in range(n - k + 1)
            if hamming(seq[i:i + k], query) <= max_mm]


def find_amplicon(
    contigs: dict[str, str], fwd: str, rev: str, max_mm: int,
    max_amplicon_bp: int = 1000,
) -> tuple[str, str, int, int] | None:
    """Find the F + RC(R) primer pair on either strand of any contig.
    Enumerates all (F, RC(R)) pairs within max_amplicon_bp of each other and
    picks the one with the fewest combined mismatches — guards against
    false-positive primer hits elsewhere in the genome (e.g. the Mpn14 primer
    has 4 sub-threshold hits before the real one)."""
    rev_rc = str(Seq(rev).reverse_complement())
    pairs: list[tuple[int, str, str, int, int]] = []  # (mm_total, contig, strand, start, end)
    for contig_id, seq in contigs.items():
        for strand, strand_seq in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
            f_hits = all_hits(strand_seq, fwd, max_mm)
            r_hits = all_hits(strand_seq, rev_rc, max_mm)
            for f_pos, f_mm in f_hits:
                for r_pos, r_mm in r_hits:
                    if r_pos < f_pos + len(fwd):
                        continue
                    if r_pos + len(rev_rc) - f_pos > max_amplicon_bp:
                        continue
                    pairs.append((f_mm + r_mm, contig_id, strand, f_pos,
                                  r_pos + len(rev_rc)))
    if not pairs:
        return None
    _, contig_id, strand, amp_start, amp_end = min(pairs, key=lambda p: p[0])
    return contig_id, strand, amp_start, amp_end


def type_sample(assembly_path: Path, config: dict) -> dict:
    contigs = {rec.id: str(rec.seq).upper()
               for rec in SeqIO.parse(assembly_path, "fasta")}
    max_mm = config.get("primer_max_mismatches", 2)

    def call(locus: dict) -> LocusCall:
        name = locus["name"]
        mlen = locus["motif_length"]
        constant = locus["constant"]

        hit = find_amplicon(contigs, locus["f_primer"], locus["r_primer"], max_mm)
        if hit is None:
            return LocusCall(name, None, mlen, None, None, None, None,
                             "missing", "primers not found in assembly")
        contig_id, strand, amp_start, amp_end = hit
        amp_len = amp_end - amp_start
        array_len = amp_len - constant

        if array_len < 0:
            return LocusCall(name, None, mlen, amp_len, array_len,
                             contig_id, strand, "low",
                             "negative array length — primer amplicon shorter than reference constant")

        copies = round(array_len / mlen)
        # Confidence: amplicon length must be (close to) an exact multiple of
        # motif length above the constant. ±2 bp tolerance for sequencing slop.
        residual = abs(array_len - copies * mlen)
        confidence = "high" if residual <= 2 else "low"
        note = ("" if confidence == "high"
                else f"array length {array_len} not a clean multiple of {mlen} (residual {residual} bp)")
        return LocusCall(name, copies, mlen, amp_len, array_len,
                         contig_id, strand, confidence, note)

    calls = [call(locus) for locus in config["loci"]]
    profile = "-".join(str(c.copies) if c.copies is not None else "X" for c in calls)
    return {"profile": profile, "loci": [asdict(c) for c in calls]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--loci", type=Path,
                        default=Path(__file__).parent / "loci.yaml")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.loci.read_text())
    result = type_sample(args.assembly, config)
    output = json.dumps(result, indent=2)
    if args.out:
        args.out.write_text(output)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
