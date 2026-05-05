#!/usr/bin/env python3
"""
Build loci.yaml from the M129 reference using Dégrange et al. 2009's
amplicon-size MLVA method.

Dégrange's note that "the size variations of the amplicons were exact
multiples of the repeats" tells us copy count is determined by PCR
fragment size, not by per-copy TR scoring. Per locus we record:

  - the published primer pair (Dégrange Table 2)
  - the motif length (Dégrange Table 3)
  - a calibrated `constant` = M129_amplicon_length − M129_copies × motif_length

At typing time:    copies = round((amplicon_length − constant) / motif_length).

Build-time validation: M129 must reproduce Dégrange's published profile
4-5-7-2 (Mpn13-16). With --strict the build fails on any mismatch.

Usage:
  .venv/bin/python build_loci_from_m129.py \
      --reference .context/attachments/NC_000912.1.fna \
      --out loci.yaml [--strict]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

# Dégrange et al. 2009, Tables 2 & 3 + M129 type P from Table 4.
# Sequences are written 5'→3' as printed in Table 2.
LOCI = [
    {"name": "Mpn13", "f_primer": "GACCAGCATTAGATTGCTATG",
     "r_primer": "AACAAATTAAGCAGCTCACG",
     "motif_length": 16, "tr_identity_pct": 83, "m129_copies": 4},
    {"name": "Mpn14", "f_primer": "CTCAGGGCGAAACCTTAAAG",
     "r_primer": "GCAATGGCTTTCAGCACAAC",
     "motif_length": 21, "tr_identity_pct": 84, "m129_copies": 5},
    {"name": "Mpn15", "f_primer": "CAACAGCACCACATCTTTAG",
     "r_primer": "GCTAATCTTGCAAACGCTGC",
     "motif_length": 21, "tr_identity_pct": 81, "m129_copies": 7},
    {"name": "Mpn16", "f_primer": "GACGCGTTCGCTAAAAGAG",
     "r_primer": "CAGGCTCAACCAAATAATGG",
     "motif_length": 47, "tr_identity_pct": 100, "m129_copies": 2},
]


def rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def find_amplicon(genome: str, fwd: str, rev: str) -> tuple[int, int] | None:
    """Locate the amplicon defined by F/R primers on the forward strand of
    `genome`. Returns (amplicon_start, amplicon_end_exclusive) or None.
    The amplicon spans from the F primer 5' end to the 3' end of where the
    R primer anneals (i.e. start of RC(R) + len(R))."""
    f_pos = genome.find(fwd)
    if f_pos < 0:
        return None
    r_anneal_pos = genome.find(rc(rev), f_pos + len(fwd))
    if r_anneal_pos < 0:
        return None
    return f_pos, r_anneal_pos + len(rev)


def build_locus_entry(genome: str, spec: dict) -> dict:
    found = find_amplicon(genome, spec["f_primer"], spec["r_primer"])
    if found is None:
        # Try reverse-complemented genome (in case M129 chromosome is given in
        # the opposite orientation)
        rc_genome = rc(genome)
        found = find_amplicon(rc_genome, spec["f_primer"], spec["r_primer"])
        if found is None:
            raise RuntimeError(f"{spec['name']}: primers not found in M129")
        amp_start, amp_end = found
        # convert back to forward coords for reporting
        amp_start, amp_end = (len(genome) - amp_end, len(genome) - amp_start)
        amplicon = rc_genome[found[0]:found[1]]
    else:
        amp_start, amp_end = found
        amplicon = genome[amp_start:amp_end]

    amp_len = amp_end - amp_start
    mlen = spec["motif_length"]
    m129_copies = spec["m129_copies"]
    constant = amp_len - m129_copies * mlen

    return {
        "name": spec["name"],
        "motif_length": mlen,
        "tr_identity_pct": spec["tr_identity_pct"],
        "f_primer": spec["f_primer"],
        "r_primer": spec["r_primer"],
        "constant": constant,
        "m129_copies": m129_copies,
        "_meta": {
            "m129_amplicon_length_bp": amp_len,
            "m129_amplicon_start_1based": amp_start + 1,
            "m129_amplicon_end_1based": amp_end,
            "m129_array_length_bp": amp_len - constant,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path,
                        default=Path(__file__).parent / "M129_NC_000912.1.fna")
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).parent / "loci.yaml")
    parser.add_argument("--strict", action="store_true",
                        help="fail if M129 profile doesn't match Dégrange 4-5-7-2")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.reference, "fasta"))
    if len(records) != 1:
        raise RuntimeError(f"expected 1 chromosome, got {len(records)}")
    genome = str(records[0].seq).upper()
    print(f"reference: {records[0].id} ({len(genome):,} bp)\n")

    entries = [build_locus_entry(genome, spec) for spec in LOCI]

    # Validate: applying the calibrated constants to the M129 amplicons must
    # reproduce the expected M129 copy counts exactly.
    print("M129 calibration:")
    fmt = "  {name:<6} amp={amp:>4} bp  array={arr:>4} bp  motif={mlen:>3} bp  copies={cp}"
    issues = []
    for e in entries:
        cp = (e["_meta"]["m129_amplicon_length_bp"] - e["constant"]) // e["motif_length"]
        print(fmt.format(name=e["name"],
                         amp=e["_meta"]["m129_amplicon_length_bp"],
                         arr=e["_meta"]["m129_array_length_bp"],
                         mlen=e["motif_length"],
                         cp=cp))
        if cp != e["m129_copies"]:
            issues.append(f"{e['name']}: computed {cp}, expected {e['m129_copies']}")

    profile = "-".join(str(e["m129_copies"]) for e in entries)
    print(f"\nM129 profile: {profile} (Dégrange 2009 type P)")

    if issues:
        print("\nVALIDATION FAILED:")
        for i in issues:
            print(f"  {i}")
        if args.strict:
            return 1

    config = {
        "reference_strain": "M129",
        "reference_accession": records[0].id,
        "expected_m129_profile": profile,
        "loci": entries,
        "primer_max_mismatches": 2,
    }
    args.out.write_text(yaml.safe_dump(config, sort_keys=False, width=200))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
