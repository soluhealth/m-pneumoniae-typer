# M. pneumoniae MLVA typer (scratch)

In silico MLVA typer for *Mycoplasma pneumoniae* implementing the
amplicon-size method of Dégrange et al. 2009 (4-locus Mpn13/14/15/16
scheme).

## Layout

```
.
├── README.md
├── requirements.txt
├── M129_NC_000912.1.fna       # vendored M. pneumoniae M129 reference
├── loci.yaml                  # built from M129
├── build_loci_from_m129.py    # one-time: calibrate constants from M129
└── mlva_type.py               # the typer (assembly → profile)
```

## Setup

```sh
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

(Tested on Python 3.9.)

## Running

```sh
.venv/bin/python mlva_type.py --assembly path/to/contigs.fasta --out profile.json
```

## How it works

Dégrange et al. 2009 (the original MLVA paper) types isolates by
multiplex PCR and capillary fragment sizing. The crucial sentence in
their methods:

> "The size variations of the amplicons were exact multiples of the repeats."

So MLVA copy count is determined by the **PCR fragment size**, not by
running TRF on each amplicon and counting copies. Per-copy detection (as
TRF does) is irrelevant to the type call — what matters is the amplicon
length.

Per locus we record (from Dégrange Tables 2 & 3):

- F primer (5'→3'), R primer (5'→3', as printed)
- motif length
- a calibrated `constant = M129_amplicon_length − M129_copies × motif_length`

At typing time:

```
amplicon_length = position(R 3' end) − position(F 5' start)
copies          = round((amplicon_length − constant) / motif_length)
```

Confidence is `high` if the array length is within 2 bp of an exact
multiple of the motif (allowing for sequencing slop), `low` otherwise.

### Primer search

Both primers are searched on both strands of every contig, allowing up to
`primer_max_mismatches` (default 2) Hamming-distance substitutions per
primer. **All** F-RC(R) pairs within 1 kb of each other are enumerated;
the pair with the fewest combined mismatches wins. This guards against
sub-threshold false-positive primer hits elsewhere in the genome — Mpn14-F
in particular has 4 such hits in M129 before the real one, and a naive
"first hit wins" implementation pairs a false F with a false R-RC,
returning the wrong amplicon.

## M129 self-validation

```
profile: 4-5-7-2 (Dégrange 2009 type P)
  Mpn13:  4 copies (amplicon 415 bp, high)
  Mpn14:  5 copies (amplicon 399 bp, high)
  Mpn15:  7 copies (amplicon 241 bp, high)
  Mpn16:  2 copies (amplicon 353 bp, high)
```

Pass `--strict` to `build_loci_from_m129.py` to make the build fail if
this calibration ever regresses.

## Why this replaces the per-copy approach

An earlier version of this typer ran pytrf (`GTRFinder`/`ATRFinder`) on
each amplicon, plus a custom canonical-motif sliding finder, and counted
copies directly. That gave **4-4-6-2** for M129 — which is what's
literally detectable as tandem repeats in NC_000912.1, but is *not* the
profile reported in the literature.

The discrepancy isn't a TRF-settings issue: in NC_000912.1, the M129
Mpn14 array has 3 perfect copies + 1 imperfect (~80%) copy + flanking
non-repeat sequence. By per-copy criteria there are 4 copies. But on
agarose / capillary, the PCR product is 399 bp, which equals
`constant + 5 × 21 bp`, so the *typing call* is 5. Dégrange's MLVA
classifies by the second number, not the first.

The amplicon-size method also handles the heavily degenerate Mpn15 array
correctly: regardless of how many "copies" pass an identity threshold,
the 241 bp amplicon corresponds to 7 motif units of 21 bp.

## Known limitations

- **Hamming-only primer matching** (no indels). Adequate for most
  *M. pneumoniae* strains; would need replacing with a proper aligner for
  highly divergent isolates.
- **Assembly-based** — short-read assemblies can collapse longer VNTRs,
  giving an undersized amplicon and an undercount. Long-read or hybrid
  assemblies are more reliable.
- **5-locus scheme not implemented.** Dégrange's original includes Mpn1
  (12 bp motif, hsdS gene), which Chalker et al. 2011 dropped due to
  instability across passages. Easy to add if needed — same primer-pair
  pattern.
