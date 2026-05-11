# M. pneumoniae typer (scratch)

In silico typers for *Mycoplasma pneumoniae*:

- **MLVA** — Dégrange et al. 2009, 4-locus Mpn13/14/15/16 amplicon-size scheme.
- **P1** — Kenri et al. 2008 multiplex amplicon-size genotyping, P1-1 vs P1-2.

## Layout

```
.
├── README.md
├── requirements.txt
├── M129_NC_000912.1.fna       # vendored M. pneumoniae M129 reference (P1-1)
├── FH_NC_017504.1.fna         # vendored M. pneumoniae FH reference   (P1-2)
├── loci.yaml                  # MLVA loci, built from M129
├── build_loci_from_m129.py    # MLVA: one-time calibration from M129
├── mlva_type.py               # MLVA typer (assembly → profile)
├── p1_loci.yaml               # P1 loci, built from M129 + FH
├── build_p1_loci.py           # P1: one-time calibration from M129 + FH
└── p1_type.py                 # P1 typer (assembly → subtype)
```

The two `*_loci.yaml` files are **build artifacts** — primers and locus
specs are authored in the corresponding `build_*.py` script (the source
of truth) and the YAML is regenerated each time the build runs. Don't
hand-edit the YAML; edit the script and rebuild.

## Setup

```sh
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

(Tested on Python 3.9.)

## Running

MLVA:

```sh
.venv/bin/python mlva_type.py --assembly path/to/contigs.fasta --out profile.json
```

P1:

```sh
.venv/bin/python p1_type.py --assembly path/to/contigs.fasta --out p1_call.json
```

## MLVA typing

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

### M129 self-validation

```
profile: 4-5-7-2 (Dégrange 2009 type P)
  Mpn13:  4 copies (amplicon 415 bp, high)
  Mpn14:  5 copies (amplicon 399 bp, high)
  Mpn15:  7 copies (amplicon 241 bp, high)
  Mpn16:  2 copies (amplicon 353 bp, high)
```

Pass `--strict` to `build_loci_from_m129.py` to make the build fail if
this calibration ever regresses.

### Why this replaces the per-copy MLVA approach

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

## P1 typing

Implements Kenri et al. 2008 (*J Med Microbiol* 57:469-475), which types
isolates by multiplex amplicon-size PCR at two regions of the P1 adhesin
gene (MPN141): **RepMP4** (typing reverses N1 / 2N2C) and **RepMP2/3**
(typing reverses R3-1 / R3-2). Expected amplicon sizes from the paper:

| locus    | P1-1 | P1-2 |
|----------|------|------|
| RepMP4   | 343  | 560  |
| RepMP2/3 | 394  | 617  |

The overall call combines both loci's calls. Kenri's R3-2V "variant 2a"
sub-classification is intentionally not implemented — we report only
the P1-1 / P1-2 dichotomy.

### Why we anchor on outer primers

The *M. pneumoniae* genome carries multiple **RepMP4 / RepMP2-3 paralogs**
(3 ADH4F binding sites and 12 MP2/3-F3 binding sites in M129) — most of
those paralogs contain mixed-subtype sequences and amplify with both the
P1-1 and P1-2 reverse primers, which would defeat naïve primer search. In
vitro, Kenri et al. handle this with a nested PCR: a first round of
gene-specific primers isolates the actual P1 gene before the typing
multiplex runs.

In silico we get the same specificity by **anchoring on a uniquely-binding
outer primer**:

- **ADH2R** for RepMP4 (1 genome-wide hit in both M129 and FH)
- **ADH3** for RepMP2/3 (1 genome-wide hit in both M129 and FH)

Each locus's typing primers are only searched within ±5 kb of the anchor.
Amplicons are accepted only when their size matches the paper-published
size within ±5 bp.

### Confidence

`high` requires `mismatches ≤ 1` and `size_dev_bp ≤ 2` (matches MLVA's
±2 bp size tolerance, plus allows a single primer-region SNP). Anything
looser falls to `low`. A locus that finds no amplicon at all is `missing`.
The overall call is `high` only if every locus called and every locus is
high.

### M129 + FH self-validation

```
M129 → P1-1   RepMP4: 343 bp (mm=0, high)   RepMP2/3: 394 bp (mm=0, high)
FH   → P1-2   RepMP4: 560 bp (mm=0, high)   RepMP2/3: 617 bp (mm=0, high)
```

Pass `--strict` to `build_p1_loci.py` to fail the build if either
calibration regresses.

## Known limitations

- **Hamming-only primer matching** (no indels). Adequate for most
  *M. pneumoniae* strains; would need replacing with a proper aligner for
  highly divergent isolates.
- **Assembly-based** — short-read assemblies can collapse longer VNTRs,
  giving an undersized amplicon and an undercount. Long-read or hybrid
  assemblies are more reliable.
- **Cross-contig amplicons not supported** (both typers). If the assembly
  breaks within the P1 gene or a VNTR locus such that the F primer and R
  primer end up on different contigs, the amplicon can't be reconstructed
  and the call returns `unknown` / `missing`. Long-read or hybrid
  assemblies avoid this.
- **MLVA 5-locus scheme not implemented.** Dégrange's original includes
  Mpn1 (12 bp motif, hsdS gene), which Chalker et al. 2011 dropped due to
  instability across passages. Easy to add if needed — same primer-pair
  pattern.
