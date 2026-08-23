# Training-data overlap of the benchmarking list

[`benchmarks.md`](benchmarks.md) reports that ssign recovers 38 of 85 experimentally
validated secreted proteins. This page asks how much of that the detection tools already
knew, by cross-referencing every benchmark protein against the published training data of
all four.

**56 of the 85 (66%) are in at least one predictor's training set, and so are 29 of the 38
ssign predicts.** Read the recall figures in `benchmarks.md` as an upper bound.

Membership alone overstates the problem, though. What matters is whether the tool that
*made* each call had seen the protein: **for 27 of the 38, at least one tool that called it
had never seen it. 11 rest entirely on a tool trained on them.**

Per-protein results: [`training_data_overlap.csv`](training_data_overlap.csv) (membership)
and [`per_tool_call_vs_training.csv`](per_tool_call_vs_training.csv) (per call).

## The four tools

| Tool | Version | Reference | Training data |
|---|---|---|---|
| MacSyFinder + TXSScan | v2 engine, TXSScan 1.1.4 | [Néron 2023](https://doi.org/10.24072/pcjournal.250), [Abby 2016](https://doi.org/10.1038/srep23080) | None. Rule-based: 280 HMM profiles (9 PFAM, 6 TIGRFAM, 265 in-house) plus co-localisation rules |
| SignalP | 6.0 | [Teufel 2022](https://doi.org/10.1038/s41587-021-01156-3) | 20,290 sequences, first 70 aa only, from UniProtKB 2018_04 extended to 2020-11-07, PROSITE and TOPDB |
| DeepLocPro | 1.0 | [Moreno 2024](https://doi.org/10.1093/bioinformatics/btae677) | 11,970 proteins: PSORTdb 4.0 experimental subset (10,241) and UniProtKB 2023_03 (2,982) |
| DeepSecE | 1.0 | [Zhang 2023](https://doi.org/10.34133/research.0258) | 2,918 proteins: 1,341 validated effectors from SecReT4, SecReT6 and BastionHub, plus 1,577 non-effectors |

MacSyFinder has no training set, no fitted weights and no train/test split, so its
exposure is measured differently. That is [below](#macsyfinder-and-txsscan).

## Which proteins are in which training set

A protein counts as in training if its accession appears in the training set, or if its
sequence is identical to a training sequence (for SignalP, to that sequence's 70-residue
truncation).

| Tool | T1SS | T2SS | T3SS | T4SS | T5SS | T6SS | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| DeepLocPro | 14/18 | 3/8 | 6/19 | 1/8 | 10/19 | 4/13 | **38/85** |
| DeepSecE (training) | 3/18 | 4/8 | 3/19 | 4/8 | 1/19 | 11/13 | **26/85** |
| SignalP 6.0 | 0/18 | 1/8 | 0/19 | 0/8 | 6/19 | 0/13 | **7/85** |
| DeepSecE (hold-out test) | 6/18 | 0/8 | 0/19 | 1/8 | 0/19 | 1/13 | **8/85** |

Across the three predictors: 56 in training, 23 with only a close homolog there, 6 with
neither. Fourteen are in two or more training sets and one, paAP, is in all three.

DeepSecE's hold-out row is separate because those 8 were not trained on. They were its
published test set, so they are not independent evidence either.

Matching on sequence as well as accession is what makes the DeepSecE and SignalP rows
right. DeepSecE ships an anonymised training FASTA and identifies proteins in a mixed
namespace of UniProt, RefSeq, GenBank, UniParc and bare gene names, so 15 of its 26 hits
have no matching accession; hpmA reaches SignalP's set under a different accession again.

## Which tool made each call

ssign calls a protein secreted on evidence from DeepLocPro, DeepSecE or SignalP. A protein
in DeepLocPro's training set only matters if DeepLocPro is what carried its call, so each
of the 38 predictions was traced back to the tools that actually fired.

| SS type | predicted | at least one caller had not seen it | which tool supplied that call |
|---|---:|---:|---|
| T1SS | 16 | 13 | DeepSecE 13, DeepLocPro 3 |
| T2SS | 1 | 1 | DeepSecE, SignalP |
| T3SS | 5 | 4 | DeepLocPro 4, DeepSecE 3 |
| T4SS | 0 | 0 | |
| T5SS | 10 | 6 | SignalP 6, DeepLocPro 1 |
| T6SS | 6 | 3 | DeepLocPro 3, DeepSecE 1 |
| **Total** | **38** | **27** | DeepSecE 18, DeepLocPro 11, SignalP 7 |

Nine of the 27 have two independent callers rather than one.

The remaining 11 fall into two clean groups. Every T5aSS autotransporter call (espP, flu,
pic, sat) comes from SignalP, which was trained on all four. Every one of the T6SS calls
(TseL, TplE, Tle1) comes from DeepSecE, which was trained on all three. The other four are
rsaA, hasA, ltxA and VirA.

DeepSecE cuts both ways: it is the largest source of independent evidence (18 of 27) and
the reason most T6SS calls are not independent.

Proteins are bridged to the run by genome coordinates using the same join that produced
`found_by_ssign`, which it reproduces for all 85 rows. The run does not record which tool
fired, so the per-tool flags re-apply ssign's own thresholds to the emitted evidence
columns.

## Does DeepLocPro depend on having seen them

DeepLocPro publishes held-out predictions from its cross-validation, so every training
protein also has predictions from models that never saw it. Scoring the 38 that way,
against ssign's own rule (`extracellular_prob >= 0.8`, widened to
extracellular-or-outer-membrane for the T5aSS/T5cSS autotransporters TXSScan models as
components):

- **26 of 38 clear the threshold** on models that had not seen them
- **12 do not**: 7 T5SS, 4 T6SS and paAP

yadA comes back at 0.081 and nadA at 0.001, so an unexposed DeepLocPro does not call them
extracellular at all. This says nothing about what the shipped model outputs for them,
which was not measured; it says only that the localisation signal for those 12 is not
reproducible without the training exposure.

## One label conflict

FhaB (P12255), a validated T5bSS secreted exoprotein and row `T5SS_09`, sits in DeepSecE's
non-effector class as `Non-effector-1520`, byte-identical over all 3,590 residues. DeepSecE
has no T5 class, so this is consistent with its own scope, but it does mean DeepSecE was
trained to call a real secreted protein non-secreted. It is the only one of the 85 in the
negative set.

## MacSyFinder and TXSScan

**Eight benchmark proteins are themselves modelled as machinery components.** All 85 were
scanned with the 280 TXSScan v1.1.4 profiles under MacSyFinder's own criteria as macsylib
implements them: gathering-threshold cutoffs, then i-evalue <= 1e-3 and profile coverage
>= 50%.

| Protein | Type | Profile |
|---|---|---|
| espP, flu, iga, pic, sat | T5aSS | `T5aSS_PF03797` |
| nadA, yadA | T5cSS | `T5cSS_PF03895` |
| VgrG-3 | T6SS | `T6SSi_tssI` |

This is biology, not a flaw. Autotransporters carry the translocator domain on the same
chain as the passenger, so for T5aSS and T5cSS the substrate and the component are one
protein, and VgrG is both a structural spike and a secreted effector. For these eight,
TXSScan finds the protein by construction.

The scan also bounds the effect. None of the 10 T5bSS exoproteins match, because
`T5bSS_translocator` models the TpsB pore rather than the secreted TpsA partner. Neither
plpD (T5dSS) nor eae (T5eSS) match, because TXSScan v1.1.4 has no model for those
subtypes. No T1SS substrate matches; only the ABC, MFP and OMF components are modelled.

**Thirty-five benchmark pairs have their system in TXSScan's reference or validation set.**
Abby 2016 Table S1 lists the validated systems whose components seeded the in-house
profiles and Table S2 the strains the models were validated on. Matching on species and
system type, 22 appear in S1 and 21 in S2. Table S2 is explicit in places: its T1SS entry
for *P. aeruginosa* PAO1 names "the aprDEF system (PA1246-1248, next to the secreted
protein aprA, PA1249)", which is row `T1SS_02`.

Table S1 has no strain rows for T3SS, Flagellum or T4SS, citing prior publications instead,
so all 19 T3SS and 8 T4SS proteins score `no` for want of an entry rather than because
their systems were checked and found absent.

## What this means for the recall figure

Two things are unaffected by any of the above.

**Reachability.** Whether a secreted protein sits within ±3 genes of its own machinery is a
fact about genome layout. That proximity can reach only 47 of 85, and essentially no T2SS
or T4SS substrate, does not depend on any tool.

**Pairing.** No tool in the stack predicts system-to-substrate pairs. DeepLocPro predicts
localisation, SignalP predicts signal peptides, DeepSecE predicts effector class without
locating the system, and TXSScan finds machinery without naming substrates. Knowing AprA is
extracellular does not tell a tool it belongs to the aprDEF system one gene away.

What is affected is per-type confidence in the recall table. T1SS and T3SS hold up best:
17 of their 21 calls have an independent caller. T5SS is the weakest, and specifically the
autotransporters, where 4 of 10 calls rest on a SignalP model trained on them.

## Method and limits

Sequences for 84 of 85 came from UniProt. YezP (`T6SS_13`) has no UniProt entry and is
translated from the RefSeq CDS at the coordinates in the benchmarking list
(`WP_002210468.1`, 117 aa). Matching used three routes: identifier (UniProt primary and
secondary accessions, RefSeq, EMBL protein IDs), exact sequence or truncated prefix, and
Smith-Waterman alignment against every training record.

Alignment uses BLOSUM62 with gap open 11 and extend 1, and homology tiers combine identity
over aligned columns with the fraction of the shorter sequence that is an identical
residue, because span coverage alone overstates gappy alignments between long repeat-rich
proteins such as the T1SS RTX toxins. The `in training` counts do not depend on any of
this; they come from accession and exact-sequence matching.

- The homolog tiers are advisory. DeepSecE de-duplicates at 60% identity and SignalP and
  DeepLocPro partition at 30%, so a close homolog may already have been separated from the
  benchmark protein by the tools' own procedures.
- DeepSecE's 1,727 non-effector sequences carry no accessions, so they were matched by
  sequence only. A near-identical variant in that negative set would be missed.
- TXSScan's seed alignments are not distributed, only the built HMMs, so system-level
  overlap is matched on species and system type rather than sequence.
- DeepLocPro's released FASTA does not record whether an entry came from PSORTdb or
  UniProt, so the two sources cannot be separated.

## Reproduce

Scripts are in [`scripts/`](scripts/). Downloaded training data stays out of git, so point
`SSIGN_OVERLAP_DIR` at a working directory; it defaults to the current one.

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

Training data sources, all public:

| Tool | Source |
|---|---|
| SignalP 6.0 | `https://services.healthtech.dtu.dk/services/SignalP-6.0/public_data/train_set.fasta` |
| DeepLocPro | `https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/data/full_dataset.fasta` |
| DeepSecE | `https://github.com/zhangyumeng1sjtu/DeepSecE`, with accessions from Data File S1 |
| TXSScan | `https://github.com/macsy-models/TXSScan` tag 1.1.4, and the Abby 2016 supplement |
