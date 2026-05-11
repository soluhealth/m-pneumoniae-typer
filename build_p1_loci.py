#!/usr/bin/env python3
"""
Build p1_loci.yaml from M129 (P1-1) and FH (P1-2) using Kenri et al. 2008's
multiplex amplicon-size genotyping scheme.

Two loci, each with an outer gene-specific anchor primer plus a multiplex
forward primer and a set of subtype-specific reverse primers:

  RepMP4    anchor = ADH2R (unique, 1 genome-wide hit → P1 RepMP4 region)
            multiplex F = ADH4F
            R = N1   (P1-1)        → 343 bp
                2N2C (P1-2)        → 560 bp
  RepMP2/3  anchor = ADH3  (unique, 1 genome-wide hit → P1 RepMP2/3 region)
            multiplex F = MP2/3-F3
            R = R3-1  (P1-1)       → 394 bp
                R3-2  (P1-2)       → 617 bp

(Kenri 2008's R3-2V variant-2a primer is not used — variant calling is
out of scope for this tool.)

Why the anchor: M. pneumoniae carries multiple RepMP4/2-3 PARALOGS (3 ADH4F
sites, 12 MP2/3-F3 sites in M129) — most have mixed-subtype sequences and
amplify with both subtype reverses. In vitro, the nested PCR's first round
isolates the actual P1 gene before the typing multiplex; in silico the
same job is done by anchoring on a uniquely-binding outer primer (ADH2R or
ADH3) and only searching the typing primers within ±ANCHOR_RADIUS bp of
that anchor.

Build-time validation: M129 anchors → P1-1 call on both loci;
FH anchors → P1-2 call on both loci.

Usage:
  .venv/bin/python build_p1_loci.py --m129 M129_NC_000912.1.fna \\
      --fh FH_NC_017504.1.fna --out p1_loci.yaml [--strict]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq


# Kenri et al. 2008, J Med Microbiol 57:469-475, Table 1.
LOCI = [
    {
        "name": "RepMP4",
        "anchor": {"primer_name": "ADH2R", "primer": "GAACTTAGCGCCAGCAACTGCCAT"},
        "f_primer": {"primer_name": "ADH4F",
                     "primer": "GACCGCATCAACCACCTTTGCGTTACG"},
        "subtype_reverses": {
            "P1-1": {"primer_name": "N1",   "primer": "CCCGGTGGTGGAAGTATTTT", "expected_bp": 343},
            "P1-2": {"primer_name": "2N2C", "primer": "TGCCTTGGTCACCGGAGTTG", "expected_bp": 560},
        },
    },
    {
        "name": "RepMP2_3",
        "anchor": {"primer_name": "ADH3", "primer": "CGAGTTTGCTGCTAACGAGT"},
        "f_primer": {"primer_name": "MP2/3-F3", "primer": "TCGACCAAGCCAACCTCCAG"},
        "subtype_reverses": {
            "P1-1": {"primer_name": "R3-1", "primer": "TTGGAATCGGACCCACTTCG", "expected_bp": 394},
            "P1-2": {"primer_name": "R3-2", "primer": "CGACGTTGTGTTTGTGCCAC", "expected_bp": 617},
        },
    },
]

ANCHOR_RADIUS = 5000   # bp window around anchor where P1 gene is searched
MAX_AMPLICON_BP = 1500
PRIMER_MAX_MM = 2
SIZE_TOLERANCE_BP = 5  # ± bp from paper-published amplicon size


def rc(s: str) -> str:
    return str(Seq(s).reverse_complement())


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def all_hits(seq: str, query: str, max_mm: int) -> list[tuple[int, int]]:
    n, k = len(seq), len(query)
    return [(i, hamming(seq[i:i + k], query)) for i in range(n - k + 1)
            if hamming(seq[i:i + k], query) <= max_mm]


def find_anchor_hits(contigs: dict[str, str], anchor: str, max_mm: int):
    """Return list of (contig_id, strand, position_on_strand_seq, mm)."""
    hits = []
    for cid, seq in contigs.items():
        for strand, ss in (("+", seq), ("-", rc(seq))):
            for pos, mm in all_hits(ss, anchor, max_mm):
                hits.append((cid, strand, pos, mm))
    return hits


def find_amplicons_in_window(window: str, fwd: str, rev: str, max_mm: int,
                             max_amplicon_bp: int):
    """Return list of (mm_total, amp_start_in_window, amp_end_in_window)."""
    rev_rc = rc(rev)
    out = []
    for fp, fmm in all_hits(window, fwd, max_mm):
        for rp, rmm in all_hits(window[fp + len(fwd):fp + max_amplicon_bp], rev_rc, max_mm):
            actual_rp = fp + len(fwd) + rp
            out.append((fmm + rmm, fp, actual_rp + len(rev_rc)))
    return out


def type_locus_from_assembly(
    contigs: dict[str, str], locus: dict,
    max_mm: int = PRIMER_MAX_MM,
    anchor_radius: int = ANCHOR_RADIUS,
    max_amplicon_bp: int = MAX_AMPLICON_BP,
    size_tolerance: int = SIZE_TOLERANCE_BP,
):
    """Type a single P1 locus from an assembly.

    Flow:
      1. Find every anchor-primer hit on both strands of every contig.
      2. Require a uniquely-best anchor (lowest mm, no ties). Tied best →
         no-call: silently picking one risks typing off a paralog.
      3. Build a ±anchor_radius window on the + strand around the anchor's
         binding site, then scan it in both orientations so the F+R pair
         is found regardless of gene direction.
      4. For each subtype reverse, keep the best-mm amplicon whose size is
         within size_tolerance of the paper-published expected_bp.
      5. Resolve: exactly one subtype matches → call; zero or both → no-call
         with the conflict surfaced in 'note'.

    Returns dict: call, confidence (high/low/missing), note, evidence
    (chosen amplicon + anchor), and anchor_hits (all considered)."""
    # find_amplicons_in_window scans up to max_amplicon_bp past the F primer;
    # if the window is shorter we'd silently truncate genuine amplicons.
    assert max_amplicon_bp <= anchor_radius * 2, \
        f"max_amplicon_bp ({max_amplicon_bp}) must fit within window 2*anchor_radius ({anchor_radius * 2})"

    # --- Step 1: locate the anchor primer ---
    anchor_hits = find_anchor_hits(contigs, locus["anchor"]["primer"], max_mm)
    if not anchor_hits:
        return {"call": None, "confidence": "missing",
                "note": f"anchor {locus['anchor']['primer_name']} not found",
                "evidence": None, "anchor_hits": []}
    anchor_hits.sort(key=lambda h: h[3])
    best_mm = anchor_hits[0][3]
    n_best = sum(1 for h in anchor_hits if h[3] == best_mm)
    # Step 2: a merely non-unique anchor (one uniquely best plus weaker
    # also-rans) is fine and just gets a note. Tied best → refuse.
    if n_best > 1:
        return {"call": None, "confidence": "low",
                "note": (f"anchor {locus['anchor']['primer_name']} ambiguous: "
                         f"{n_best} hits tied at mm={best_mm}"),
                "evidence": None, "anchor_hits": anchor_hits}
    note = ""
    if len(anchor_hits) > 1:
        note = (f"anchor {locus['anchor']['primer_name']} non-unique "
                f"({len(anchor_hits)} hits); used uniquely-best hit at mm={best_mm}")

    # --- Step 3: build the search window ---
    cid, strand, anchor_pos, anchor_mm = anchor_hits[0]
    contig_seq = contigs[cid]
    # Map anchor_pos back to + strand coords: a hit at position p on rc(seq)
    # of length L corresponds to the anchor's binding site [L-p-k, L-p] on +.
    if strand == "+":
        center_fwd = anchor_pos
    else:
        center_fwd = len(contig_seq) - anchor_pos - len(locus["anchor"]["primer"])
    win_start = max(0, center_fwd - anchor_radius)
    win_end = min(len(contig_seq), center_fwd + len(locus["anchor"]["primer"]) + anchor_radius)
    window_fwd = contig_seq[win_start:win_end]
    # PCR doesn't care which strand the primers bind, so the genuine amplicon
    # may sit on either orientation of the window. Scan both. (A real amplicon
    # appears in only one — these primers aren't palindromic.)
    windows = (("+", window_fwd), ("-", rc(window_fwd)))

    # --- Step 4: find best amplicon per subtype reverse ---
    fwd = locus["f_primer"]["primer"]
    matches = {}
    for name, sp in locus["subtype_reverses"].items():
        close: list[tuple[int, int, int, str]] = []
        for w_strand, w_seq in windows:
            amps = find_amplicons_in_window(w_seq, fwd, sp["primer"], max_mm, max_amplicon_bp)
            for a in amps:
                size_dev = abs((a[2] - a[1]) - sp["expected_bp"])
                if size_dev <= size_tolerance:
                    close.append(a + (w_strand,))
        if close:
            mm, fp, ep, w_strand = min(close, key=lambda a: a[0])
            size_dev = abs((ep - fp) - sp["expected_bp"])
            matches[name] = {"amp_bp": ep - fp, "mm": mm, "size_dev_bp": size_dev,
                             "win_strand": w_strand, "f_pos_in_window": fp,
                             "primer_name": sp["primer_name"]}

    if len(matches) == 0:
        return {"call": None, "confidence": "missing",
                "note": (note + "; " if note else "") + "no subtype amplicon at expected size",
                "evidence": None, "anchor_hits": anchor_hits}
    if len(matches) > 1:
        return {"call": None, "confidence": "low",
                "note": (note + "; " if note else "") + f"both subtypes amplify: {list(matches)}",
                "evidence": {"matches": matches}, "anchor_hits": anchor_hits}

    [(call, ev)] = matches.items()
    # Tolerance mirrors mlva_type.py's ±2 bp size-deviation rule, plus allows
    # up to 1 primer mismatch (clinical isolates routinely carry the odd SNP
    # in primer-binding regions; 2 mm is at the search threshold and worth
    # flagging).
    confidence = "high" if (ev["mm"] <= 1 and ev["size_dev_bp"] <= 2) else "low"

    return {
        "call": call,
        "confidence": confidence,
        "note": note,
        "evidence": {"chosen_subtype": call, "subtype_amplicon": ev,
                     "anchor": {"contig": cid, "strand": strand,
                                "anchor_pos_in_strand": anchor_pos,
                                "anchor_mm": anchor_mm}},
        "anchor_hits": anchor_hits,
    }


def load_assembly(path: Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m129", type=Path, default=Path(__file__).parent / "M129_NC_000912.1.fna")
    parser.add_argument("--fh", type=Path, default=Path(__file__).parent / "FH_NC_017504.1.fna")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "p1_loci.yaml")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    refs = {"M129": (load_assembly(args.m129), "P1-1"),
            "FH":   (load_assembly(args.fh),   "P1-2")}

    issues: list[str] = []
    out_loci = []
    for locus in LOCI:
        print(f"== {locus['name']} (anchor={locus['anchor']['primer_name']}, "
              f"F={locus['f_primer']['primer_name']}) ==")
        per_ref = {}
        for ref_name, (contigs, expect) in refs.items():
            res = type_locus_from_assembly(contigs, locus)
            per_ref[ref_name] = res
            print(f"  {ref_name} (expect {expect}): call={res['call']}, "
                  f"anchor_hits={len(res['anchor_hits'])}")
            if res["evidence"] and res["evidence"].get("subtype_amplicon"):
                a = res["evidence"]["subtype_amplicon"]
                print(f"    {res['call']} amp={a['amp_bp']} bp (mm={a['mm']})")
            if res["note"]:
                print(f"    note: {res['note']}")
            if res["call"] != expect:
                issues.append(f"{locus['name']}/{ref_name}: called {res['call']}, expected {expect}")
        out_loci.append({
            "name": locus["name"],
            "anchor": locus["anchor"],
            "f_primer": locus["f_primer"],
            "subtype_reverses": locus["subtype_reverses"],
            "_meta": {ref: {"call": per_ref[ref]["call"],
                            "evidence": per_ref[ref]["evidence"]}
                      for ref in refs},
        })
        print()

    if issues:
        print("VALIDATION ISSUES:")
        for i in issues:
            print(f"  - {i}")
        if args.strict:
            return 1
    else:
        print("validation: M129 → P1-1 and FH → P1-2 on both loci.")

    config = {
        "scheme": "Kenri et al. 2008 (J Med Microbiol 57:469-475) multiplex amplicon-size, anchored on outer primer",
        "references": {"P1-1": "M129 (NC_000912.1)", "P1-2": "FH (NC_017504.1)"},
        "primer_max_mismatches": PRIMER_MAX_MM,
        "max_amplicon_bp": MAX_AMPLICON_BP,
        "anchor_radius_bp": ANCHOR_RADIUS,
        "size_tolerance_bp": SIZE_TOLERANCE_BP,
        "loci": out_loci,
    }
    args.out.write_text(yaml.safe_dump(config, sort_keys=False, width=200))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
