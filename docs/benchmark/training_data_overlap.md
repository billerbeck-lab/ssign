# Training-data overlap of the benchmarking list

[`benchmarks.md`](benchmarks.md) closes with a caveat: some of the 85 secretion-system
and secreted-protein pairs may have been used to train the tools ssign calls. This page
resolves that caveat. Every benchmark protein was cross-referenced against the actual
published training data of all four detection tools.

**The short answer: 56 of 85 (66%) are in at least one predictor's training set, and
among the 38 proteins ssign successfully predicts, 29 (76%) are.** The recall figure in
`benchmarks.md` should be read as an upper bound, not an estimate of performance on
novel proteins.

Per-protein results: [`training_data_overlap.csv`](training_data_overlap.csv).

## The four tools and what "training data" means for each

| Tool | Version | Reference | What it was built from |
|---|---|---|---|
| MacSyFinder + TXSScan | v2 engine, TXSScan 1.1.4 | [Néron 2023](https://doi.org/10.24072/pcjournal.250), [Abby 2016](https://doi.org/10.1038/srep23080) | Not machine learning. 280 HMM profiles (9 PFAM, 6 TIGRFAM, 265 built in-house) plus co-localisation rules |
| SignalP | 6.0 | [Teufel 2022](https://doi.org/10.1038/s41587-021-01156-3) | 20,290 sequences, first 70 aa only. UniProtKB 2018_04 extended to 2020-11-07, plus PROSITE and TOPDB |
| DeepLocPro | 1.0 | [Moreno 2024](https://doi.org/10.1093/bioinformatics/btae677) | 11,970 proteins. PSORTdb 4.0 experimental subset (10,241) plus UniProtKB 2023_03 with experimental localisation evidence (2,982) |
| DeepSecE | 1.0 | [Zhang 2023](https://doi.org/10.34133/research.0258) | 2,918 proteins: 1,341 experimentally validated effectors from SecReT4, SecReT6 and BastionHub, plus 1,577 non-effectors |

MacSyFinder is the odd one out. It has no training set, no fitted weights and no
train/test split, so the equivalent question is which reference systems seeded its HMM
profiles and which strains the finished models were validated on. That is answered
separately below.

## Overlap with the three machine-learning predictors

A protein counts as **in training** if its accession appears in the training set, or if
its amino-acid sequence is identical to a training sequence (or, for SignalP, to that
sequence's 70-residue truncation). Matching on sequence as well as accession matters
twice over: DeepSecE ships an anonymised training FASTA and identifies proteins in a
mixed namespace of UniProt, RefSeq, GenBank, UniParc and bare gene names, so 15 of its 26
hits are invisible to accession matching; and hpmA reaches SignalP's set under a
different accession entirely.

| Tool | T1SS | T2SS | T3SS | T4SS | T5SS | T6SS | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| DeepLocPro | 14/18 | 3/8 | 6/19 | 1/8 | 10/19 | 4/13 | **38/85** |
| DeepSecE (training) | 3/18 | 4/8 | 3/19 | 4/8 | 1/19 | 11/13 | **26/85** |
| SignalP 6.0 | 0/18 | 1/8 | 0/19 | 0/8 | 6/19 | 0/13 | **7/85** |
| DeepSecE (hold-out test) | 6/18 | 0/8 | 0/19 | 1/8 | 0/19 | 1/13 | **8/85** |

Combined across the three predictors: 56 direct, 23 with only a close homolog in
training, 6 with neither. Fourteen are in two or more predictors' training sets, and one,
paAP (`T2SS_04`), is in all three.

DeepSecE's hold-out row is listed separately because those 8 proteins were not trained
on. They were its published test set, so they are still not independent evidence, but
the leakage is weaker.

SignalP's 7 are worth singling out. Six are T5SS, which is where SignalP is load-bearing:
ssign calls a T5 component a substrate only if DeepLocPro localises it **or** SignalP
finds a Sec signal peptide.

### Does DeepLocPro actually depend on having seen them?

DeepLocPro publishes held-out predictions from its nested cross-validation, so for each
training protein there are four predictions from models that never saw it. Applying
ssign's own rule to those held-out predictions (`extracellular_prob >= 0.8`, widened to
extracellular-or-outer-membrane for the T5aSS/T5cSS autotransporters that TXSScan models
as components):

- all 38 have held-out predictions
- **26 of 38 (68%) still pass.** DeepLocPro would have called them secreted anyway
- **12 of 38 flip to negative.** For these, ssign's DeepLocPro signal depends on training exposure

The 12 that flip are 7 T5SS, 4 T6SS and paAP. Those are the two types where ssign's
substrate call leans hardest on localisation.

### One label conflict worth knowing

FhaB (P12255), a validated T5bSS secreted exoprotein and row `T5SS_09` of the benchmark
list, sits in DeepSecE's **non-effector** class as `Non-effector-1520`. The sequences are
byte-identical over all 3,590 residues. DeepSecE has no T5 class at all, so this is
consistent with its own scope, but it does mean DeepSecE was trained to call a real
secreted protein non-secreted. It is the only one of the 85 in the negative set.

## MacSyFinder and TXSScan

Two separate exposures, both established from the shipped models and the Abby 2016
supplement rather than argued from the papers' prose.

**1. Eight benchmark proteins are themselves modelled as machinery components.** All 85
sequences were scanned with the 280 TXSScan v1.1.4 profiles, reproducing MacSyFinder's
hit criteria as macsylib actually implements them: gathering-threshold cutoffs
(`--cut_ga` is on by default, and 278 of the 280 profiles carry a GA), then
i-evalue <= 1e-3 and profile coverage >= 50%.

| Protein | Type | Profile |
|---|---|---|
| espP, flu, iga, pic, sat | T5aSS | `T5aSS_PF03797` |
| nadA, yadA | T5cSS | `T5cSS_PF03895` |
| VgrG-3 | T6SS | `T6SSi_tssI` |

This is expected biology, not a flaw. Autotransporters carry the translocator domain on
the same polypeptide as the passenger, so for T5aSS and T5cSS the "substrate" and the
"component" are one protein. VgrG is both a structural spike and a secreted effector.
For these eight, detection is not an independent prediction: TXSScan finds the protein
by design.

The scan also confirms the limits. None of the 10 T5bSS exoproteins match, because
`T5bSS_translocator` models TpsB (the pore), not TpsA (the secreted partner). Neither
plpD (T5dSS) nor eae (T5eSS) match, because TXSScan v1.1.4 has no model for those
subtypes. No T1SS substrate matches; TXSScan models only the ABC, MFP and OMF components.

**2. Thirty-five benchmark pairs have their system in TXSScan's reference or validation
set.** Abby 2016 Table S1 lists the experimentally validated systems whose components
seeded the in-house profiles, and Table S2 lists the strains the finished models were
validated on. Matching on species and system type, 22 pairs appear in Table S1 and 21 in
Table S2 (35 in either). Table S2 is explicit in places: its T1SS entry for
*P. aeruginosa* PAO1 names "the aprDEF system (PA1246-1248, next to the secreted protein
aprA, PA1249)", which is benchmark row `T1SS_02`.

So for over a third of the list, the machinery half of the pair was a development case
for the tool that detects it. One caveat on reading a `no` here: Table S1 has no strain
rows at all for T3SS, Flagellum or T4SS, citing prior publications instead, so all 19
T3SS and 8 T4SS proteins score `no` for want of a table entry rather than because their
systems were checked and found absent.

## What this does to the recall figure

Splitting the 38 proteins ssign predicts by whether a predictor was trained on them:

| SS type | predicted | in training data | not in training data |
|---|---:|---:|---:|
| T1SS | 16 | 13 | 3 |
| T2SS | 1 | 1 | 0 |
| T3SS | 5 | 2 | 3 |
| T4SS | 0 | 0 | 0 |
| T5SS | 10 | 8 | 2 |
| T6SS | 6 | 5 | 1 |
| **Total** | **38** | **29** | **9** |

The nine in the right-hand column are not exposure-free: eight still have a close
homolog in a training set, and the ninth, Tle1 of the *E. coli* 042 Sci-1 system, has its
system in both TXSScan's profile-seed set and its validation set. Counting all four tools, **not one of the 38 proteins
ssign predicts is free of every documented exposure.** Only 3 of the 85 are (NleG, BepA
and Ats-1), and ssign predicts none of them.

This does not invalidate the benchmark, and two things about it still hold:

- **Reachability is unaffected.** Whether a secreted protein sits within ±3 genes of its
  own machinery is a fact about genome layout, not about any tool. The finding that
  proximity can reach only 47 of 85, and essentially no T2SS or T4SS substrate, is
  untouched by training-data overlap.
- **Pairing is still a novel task.** No tool in the stack predicts system-to-substrate
  pairs. DeepLocPro predicts localisation, SignalP predicts signal peptides, DeepSecE
  predicts effector class without locating the system, and TXSScan finds machinery
  without naming substrates. Memorising that AprA is extracellular does not tell a tool
  that AprA belongs to the aprDEF system one gene away.

What it does mean is that **the 45% recall figure is an upper bound.** On genuinely novel
proteins, the DeepLocPro held-out result suggests roughly a third of the localisation
calls would not survive, and T5SS and T6SS would degrade most.

## Method and limits

Sequences for 84 of 85 proteins came from UniProt. YezP (`T6SS_13`) has no UniProt entry
and is translated from the RefSeq CDS at the coordinates in the benchmarking list
(`WP_002210468.1`, 117 aa). Matching used three routes: identifier (UniProt primary and
secondary accessions, RefSeq, EMBL protein IDs), exact sequence or truncated prefix, and
Smith-Waterman alignment against every training record.

Three details that changed the numbers, recorded because they are easy to get wrong:

- SignalP distributes its datasets in a **3-line** format: header, sequence, then an
  equal-length per-residue annotation string. Parsed as ordinary FASTA, every training
  sequence silently gains 70 junk characters, which kills the prefix route that exists
  precisely because SignalP truncates to 70 aa. Handling it correctly moves SignalP from
  6 hits to 7.
- Alignment is Smith-Waterman (pyopal 0.7.3) under BLAST's protein defaults, **BLOSUM62
  with gap open 11, extend 1**. pyopal's own defaults are BLOSUM50 with gap open 3, cheap
  enough that local alignments stitch through unrelated segments and inflate the homolog
  tiers.
- Candidates for full alignment are the union of two rankings, by raw score and by score
  per residue of the shorter sequence. Raw score alone favours long targets and pushes
  the true best hit out of the candidate set. On the DeepSecE training set this
  reproduces an exhaustive all-against-all search verdict for verdict.

Homology tiers use identity over aligned columns together with the fraction of the
shorter sequence that is an identical residue, because span coverage alone overstates
gappy alignments between long repeat-rich proteins, which is exactly what the T1SS RTX
toxins and T5bSS exoproteins are. The `in training` counts do not depend on any of the
alignment settings: they come from accession and exact-sequence matching.

Known limits:

- DeepSecE's 1,727 non-effector sequences carry no accessions anywhere, so they were
  matched by sequence only. A benchmark protein present as a near-identical variant in
  that negative set would be missed.
- TXSScan's seed alignments are not distributed, only the built HMMs. Table S1 gives
  strain-level granularity, not the accession list behind each profile, so the
  system-level overlap is matched on species and system type rather than on sequence.
- DeepLocPro's released FASTA does not record whether an entry came from PSORTdb or
  UniProt, so the two sources cannot be separated.
- The homolog tiers are advisory. DeepSecE de-duplicates at 60% identity and SignalP and
  DeepLocPro partition at 30%, so a "close homolog" may or may not have been separated
  from the benchmark protein by the tools' own procedures.

## Reproduce

Scripts are in [`scripts/`](scripts/). Downloaded training data is large and stays out of
git, so point `SSIGN_OVERLAP_DIR` at a working directory; it defaults to the current one.

```bash
pip install -r scripts/requirements.txt      # pyopal is pinned; results depend on it
export SSIGN_OVERLAP_DIR=/path/to/workdir    # holds raw/ (downloads) and work/ (intermediates)
python scripts/fetch_benchmark_uniprot.py      # sequences + accessions for the 85
python scripts/match_training_set.py --db <training.fasta> --tool <name> --out <out.tsv>
python scripts/scan_txsscan_profiles.py        # 280 TXSScan HMMs vs the 85 proteins
python scripts/match_txsscan_organisms.py      # benchmark organisms vs Abby 2016 S1/S2
python scripts/deeplocpro_heldout_check.py     # held-out CV predictions vs ssign's rule
python scripts/aggregate_overlap.py            # writes training_data_overlap.csv
```

Training data sources, all public:

- SignalP 6.0: `https://services.healthtech.dtu.dk/services/SignalP-6.0/public_data/train_set.fasta`
- DeepLocPro: `https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/data/full_dataset.fasta`
- DeepSecE: `https://github.com/zhangyumeng1sjtu/DeepSecE` (`data/Train-2918.fasta`), with
  accessions from Data File S1 of the paper
- TXSScan: `https://github.com/macsy-models/TXSScan` tag 1.1.4, and the Abby 2016
  supplement for Tables S1 and S2
