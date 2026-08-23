# Training-data overlap of the benchmarking list

[`benchmarks.md`](benchmarks.md) reports that ssign recovers 38 of 85 validated secreted
proteins. This page measures how much of that the detection tools already knew.

**56 of the 85 (66%) are in at least one predictor's training set, and so are 29 of the 38
ssign predicts.** Read the recall figures as an upper bound.

Membership overstates it, though. What matters is whether the tool that *made* each call
had seen the protein: **27 of the 38 have a caller that never saw it; 11 rest entirely on a
tool trained on them.**

Per-protein: [`training_data_overlap.csv`](training_data_overlap.csv) (membership),
[`per_tool_call_vs_training.csv`](per_tool_call_vs_training.csv) (per call).

## The four tools

| Tool | Version | Reference | Training data |
|---|---|---|---|
| MacSyFinder + TXSScan | v2, TXSScan 1.1.4 | [Néron 2023](https://doi.org/10.24072/pcjournal.250), [Abby 2016](https://doi.org/10.1038/srep23080) | None. Rule-based: 280 HMM profiles (9 PFAM, 6 TIGRFAM, 265 in-house) plus co-localisation rules |
| SignalP | 6.0 | [Teufel 2022](https://doi.org/10.1038/s41587-021-01156-3) | 20,290 sequences, first 70 aa only, from UniProtKB 2018_04 extended to 2020-11-07, PROSITE and TOPDB |
| DeepLocPro | 1.0 | [Moreno 2024](https://doi.org/10.1093/bioinformatics/btae677) | 11,970 proteins: PSORTdb 4.0 experimental subset (10,241) and UniProtKB 2023_03 (2,982) |
| DeepSecE | 1.0 | [Zhang 2023](https://doi.org/10.34133/research.0258) | 2,918 proteins: 1,341 validated effectors from SecReT4, SecReT6 and BastionHub, plus 1,577 non-effectors |

MacSyFinder has no training set, so its exposure is measured differently
([below](#macsyfinder-and-txsscan)).

## Which proteins are in which training set

In training means the accession appears in the training set, or the sequence is identical
to a training sequence (for SignalP, to its 70-residue truncation).

| Tool | T1SS | T2SS | T3SS | T4SS | T5SS | T6SS | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| DeepLocPro | 14/18 | 3/8 | 6/19 | 1/8 | 10/19 | 4/13 | **38/85** |
| DeepSecE (training) | 3/18 | 4/8 | 3/19 | 4/8 | 1/19 | 11/13 | **26/85** |
| SignalP 6.0 | 0/18 | 1/8 | 0/19 | 0/8 | 6/19 | 0/13 | **7/85** |
| DeepSecE (hold-out test) | 6/18 | 0/8 | 0/19 | 1/8 | 0/19 | 1/13 | **8/85** |

Across the three predictors: 56 in training, 23 with only a close homolog there, 6 with
neither. Fourteen are in two or more sets and one, paAP, is in all three. DeepSecE's
hold-out row is separate because those 8 were its test set, not trained on, which is weaker
exposure but not independence.

Sequence matching is what makes the DeepSecE and SignalP rows right: DeepSecE ships an
anonymised training FASTA in a mixed accession namespace, so 15 of its 26 hits have no
matching accession, and hpmA reaches SignalP's set under a different accession again.

## Which tool made each call

A protein in DeepLocPro's training set only matters if DeepLocPro carried its call, so each
of the 38 predictions was traced to the tools that fired.

| SS type | predicted | caller had not seen it | which tool |
|---|---:|---:|---|
| T1SS | 16 | 13 | DeepSecE 13, DeepLocPro 3 |
| T2SS | 1 | 1 | DeepSecE, SignalP |
| T3SS | 5 | 4 | DeepLocPro 4, DeepSecE 3 |
| T4SS | 0 | 0 | |
| T5SS | 10 | 6 | SignalP 6, DeepLocPro 1 |
| T6SS | 6 | 3 | DeepLocPro 3, DeepSecE 1 |
| **Total** | **38** | **27** | DeepSecE 18, DeepLocPro 11, SignalP 7 |

Nine of the 27 have two independent callers.

The 11 without form two groups: every T5aSS autotransporter call (espP, flu, pic, sat)
comes from SignalP, trained on all four; every T6SS one (TseL, TplE, Tle1) from DeepSecE,
trained on all three. The rest are rsaA, hasA, ltxA and VirA. DeepSecE cuts both ways: the
largest source of independent evidence, and the reason most T6SS calls are not independent.

Proteins are bridged to the run by genome coordinates using the join that produced
`found_by_ssign`, which it reproduces for all 85. The run does not record which tool fired,
so the flags re-apply ssign's thresholds to the emitted evidence columns.

## Does DeepLocPro depend on having seen them

DeepLocPro publishes held-out predictions from its cross-validation, so every training
protein also has predictions from models that never saw it. Scored that way against ssign's
rule (`extracellular_prob >= 0.8`, widened to extracellular-or-outer-membrane for the
T5aSS/T5cSS autotransporters TXSScan models as components), **26 of 38 clear the threshold
and 12 do not** (7 T5SS, 4 T6SS, paAP). yadA comes back at 0.081 and nadA at 0.001.

This says nothing about what the shipped model outputs, which was not measured, only that
the localisation signal for those 12 is not reproducible without the training exposure.

## One label conflict

FhaB (P12255), a validated T5bSS exoprotein and row `T5SS_09`, sits in DeepSecE's
non-effector class as `Non-effector-1520`, byte-identical over 3,590 residues. DeepSecE has
no T5 class, so this is within its scope, but it was trained to call a real secreted protein
non-secreted. It is the only one of the 85 in the negative set.

## MacSyFinder and TXSScan

**Eight benchmark proteins are themselves modelled as machinery components.** All 85 were
scanned with the 280 TXSScan v1.1.4 profiles under MacSyFinder's criteria as macsylib
implements them: gathering-threshold cutoffs, then i-evalue <= 1e-3 and profile coverage
>= 50%.

| Protein | Type | Profile |
|---|---|---|
| espP, flu, iga, pic, sat | T5aSS | `T5aSS_PF03797` |
| nadA, yadA | T5cSS | `T5cSS_PF03895` |
| VgrG-3 | T6SS | `T6SSi_tssI` |

This is biology, not a flaw. Autotransporters carry the translocator domain on the same
chain as the passenger, and VgrG is both a structural spike and a secreted effector, so for
these eight TXSScan finds the protein by construction.

The scan also bounds the effect. No T5bSS exoprotein matches, because `T5bSS_translocator`
models the TpsB pore rather than the secreted TpsA. Neither plpD (T5dSS) nor eae (T5eSS)
match, as v1.1.4 has no model for those subtypes. No T1SS substrate matches; only ABC, MFP
and OMF are modelled.

**Thirty-five pairs have their system in TXSScan's reference or validation set.** Abby 2016
Table S1 lists the systems whose components seeded the in-house profiles and Table S2 the
strains the models were validated on; matching on species and system type, 22 appear in S1
and 21 in S2. S2 is explicit in places: its T1SS entry for *P. aeruginosa* PAO1 names "the
aprDEF system (PA1246-1248, next to the secreted protein aprA, PA1249)", which is row
`T1SS_02`. S1 has no strain rows for T3SS, Flagellum or T4SS, so those 27 proteins score
`no` for want of an entry rather than a checked absence.

## What this means for the recall figure

Two things are unaffected. **Reachability** is a fact about genome layout: that proximity
reaches only 47 of 85, and essentially no T2SS or T4SS substrate, does not depend on any
tool. **Pairing** is a task no tool performs. DeepLocPro predicts localisation, SignalP
signal peptides, DeepSecE effector class without locating the system, TXSScan machinery
without naming substrates. Knowing AprA is extracellular does not tell a tool it belongs to
the aprDEF system one gene away.

What is affected is per-type confidence. T1SS and T3SS hold up best, 17 of their 21 calls
having an independent caller. T5SS is weakest, and specifically the autotransporters, where
4 of 10 calls rest on a SignalP model trained on them.

## Method and limits

Sequences for 84 of 85 came from UniProt; YezP (`T6SS_13`) has no entry and is translated
from the RefSeq CDS at the listed coordinates (`WP_002210468.1`, 117 aa). Matching used
identifier (UniProt primary and secondary, RefSeq, EMBL protein IDs), exact sequence or
truncated prefix, and Smith-Waterman against every training record with BLOSUM62, gap open
11, extend 1. The `in training` counts depend only on the first two routes.

- Homolog tiers are advisory. DeepSecE de-duplicates at 60% identity, SignalP and
  DeepLocPro partition at 30%, so a close homolog may already have been separated by the
  tools' own procedures.
- DeepSecE's 1,727 non-effectors carry no accessions, so were matched by sequence only. A
  near-identical variant there would be missed.
- TXSScan's seed alignments are not distributed, so system overlap is matched on species
  and system type, not sequence.
- DeepLocPro's FASTA does not record PSORTdb vs UniProt provenance.

## Reproduce

Scripts are in [`scripts/`](scripts/), which also carry the training-data download URLs.
Point `SSIGN_OVERLAP_DIR` at a working directory; it defaults to the current one.

```bash
pip install -r scripts/requirements.txt
export SSIGN_OVERLAP_DIR=/path/to/workdir    # holds raw/ and work/
python scripts/fetch_benchmark_uniprot.py
python scripts/match_training_set.py --db <training.fasta> --tool <name> --out <out.tsv>
python scripts/scan_txsscan_profiles.py
python scripts/match_txsscan_organisms.py
python scripts/deeplocpro_heldout_check.py
python scripts/aggregate_overlap.py            # writes training_data_overlap.csv
python scripts/per_tool_call_vs_training.py    # writes per_tool_call_vs_training.csv
```
