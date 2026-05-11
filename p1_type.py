#!/usr/bin/env python3
"""
M. pneumoniae P1 typer — Kenri et al. 2008 multiplex amplicon-size method,
anchored on outer gene-specific primers to avoid the multiple RepMP4/2-3
paralogs scattered across the genome.

Per locus:
  1. Locate the outer anchor primer (ADH2R for RepMP4, ADH3 for RepMP2/3).
     The anchor is unique enough to pin down the actual P1 gene region.
  2. Search the multiplex F primer + each subtype-specific R primer within
     a ±ANCHOR_RADIUS window around the anchor, on both strand orientations.
  3. The R primer that produces an amplicon of the paper-expected size
     determines the subtype call.

The overall P1 subtype is the consensus of the two loci's calls
(disagreement is reported as 'mixed').

Usage:
  .venv/bin/python p1_type.py --assembly contigs.fasta \\
      --loci p1_loci.yaml --out p1_call.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from build_p1_loci import type_locus_from_assembly, load_assembly


def combine_locus_calls(per_locus: dict) -> tuple[str, str, str]:
    """Reduce per-locus calls to an overall P1 subtype, overall confidence
    ('high' / 'low' / 'missing'), and a note."""
    n_total = len(per_locus)
    called = {n: r for n, r in per_locus.items() if r["call"] is not None}
    if not called:
        return "unknown", "missing", "no locus produced a call"
    unique = set(r["call"] for r in called.values())
    if len(unique) > 1:
        return "mixed", "low", \
            f"loci disagreed: {dict((k, v['call']) for k, v in per_locus.items())}"
    call = next(iter(unique))
    if len(called) < n_total:
        return call, "low", f"partial call: only {len(called)}/{n_total} loci produced a call"
    if any(r["confidence"] != "high" for r in called.values()):
        return call, "low", "low-confidence amplicon (mismatches or off-by-bp size)"
    return call, "high", ""


def type_sample(assembly_path: Path, config: dict) -> dict:
    contigs = load_assembly(assembly_path)
    locus_kwargs = {
        "max_mm": config.get("primer_max_mismatches", 2),
        "anchor_radius": config.get("anchor_radius_bp", 5000),
        "max_amplicon_bp": config.get("max_amplicon_bp", 1500),
        "size_tolerance": config.get("size_tolerance_bp", 5),
    }
    per_locus = {}
    for locus in config["loci"]:
        per_locus[locus["name"]] = type_locus_from_assembly(contigs, locus, **locus_kwargs)

    subtype, confidence, note = combine_locus_calls(per_locus)
    return {
        "subtype": subtype,
        "confidence": confidence,
        "note": note,
        "loci": {name: {
            "call": r["call"],
            "confidence": r["confidence"],
            "note": r["note"],
            "amplicon_bp": (r["evidence"]["subtype_amplicon"]["amp_bp"]
                            if r.get("evidence") and r["evidence"].get("subtype_amplicon") else None),
            "mismatches": (r["evidence"]["subtype_amplicon"]["mm"]
                           if r.get("evidence") and r["evidence"].get("subtype_amplicon") else None),
            "size_dev_bp": (r["evidence"]["subtype_amplicon"]["size_dev_bp"]
                            if r.get("evidence") and r["evidence"].get("subtype_amplicon") else None),
            "anchor_hits": len(r["anchor_hits"]),
        } for name, r in per_locus.items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True, type=Path)
    parser.add_argument("--loci", type=Path,
                        default=Path(__file__).parent / "p1_loci.yaml")
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
