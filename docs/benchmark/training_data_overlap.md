# Training-data overlap of the benchmarking list

[`benchmarks.md`](benchmarks.md) reports that ssign recovers 38 of 83 validated secreted
proteins. This page measures how much of that the detection tools were exposed to in training.

**55 of the 83 are in at least one predictor's training set, and so are 29 of the 38
ssign predicts.** 

However, what matters is whether the tool that *made* each call
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

| Tool | T1SS | T2SS | T3SS | T4SS | T5SS | T6SS | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| DeepLocPro | 14/18 | 3/8 | 6/19 | 1/8 | 9/17 | 4/13 | **37/83** |
| DeepSecE (training) | 3/18 | 4/8 | 3/19 | 4/8 | 1/17 | 11/13 | **26/83** |
| SignalP 6.0 | 0/18 | 1/8 | 0/19 | 0/8 | 6/17 | 0/13 | **7/83** |
| DeepSecE (hold-out test) | 6/18 | 0/8 | 0/19 | 1/8 | 0/17 | 1/13 | **8/83** |

## Which tool made each call

| SS type | predicted | caller had not seen it | which tool |
|---|---:|---:|---|
| T1SS | 16 | 13 | DeepSecE 13, DeepLocPro 3 |
| T2SS | 1 | 1 | DeepSecE, SignalP |
| T3SS | 5 | 4 | DeepLocPro 4, DeepSecE 3 |
| T4SS | 0 | 0 | |
| T5SS | 10 | 6 | SignalP 6, DeepLocPro 1 |
| T6SS | 6 | 3 | DeepLocPro 3, DeepSecE 1 |
| **Total** | **38** | **27** | DeepSecE 18, DeepLocPro 11, SignalP 7 |

## Does DeepLocPro depend on having seen them

DeepLocPro publishes held-out predictions from its cross-validation, so every training
protein also has predictions from models that never saw it. Scored that way against ssign's rule
(`extracellular_prob >= 0.8`, or for the T5aSS/T5cSS autotransporters either extracellular or
outer-membrane above 0.8, since the barrel stays in the membrane its passenger threads through),
**26 of 37 clear the threshold and 11 do not** (6 T5SS, 4 T6SS, paAP). yadA comes back at 0.081
and nadA at 0.048.

This says nothing about what the shipped model outputs, which was not measured, only that
the localisation signal for those 11 is not reproducible without the training exposure.

## MacSyFinder and TXSScan

**Eight benchmark proteins are themselves modelled as machinery components.** All 83 were
scanned with the 280 TXSScan v1.1.4 profiles under MacSyFinder's criteria as macsylib
implements them: gathering-threshold cutoffs, then i-evalue <= 1e-3 and profile coverage
>= 50%.

The following present:

| Protein | Type | Profile |
|---|---|---|
| espP, flu, iga, pic, sat | T5aSS | `T5aSS_PF03797` |
| nadA, yadA | T5cSS | `T5cSS_PF03895` |
| VgrG-3 | T6SS | `T6SSi_tssI` |

## What this means for the benchmark test

Though 29/38 of the found pairs has at least one tool trained on the protein found, 27/38 had at least one tool called that wasn't trained on it, meaning the tool that was trained on those 27 proteins could be removed from ssign and the pair would still be found. Additionally, the **pairing** a secretion system to secreted protein is a task no tool performs on its own and the test of that is unaffected by the fact that some of these tested proteins were in training datasets.

T1SS and T3SS hold up best, 17 of their 21 calls
having an independent caller. T5SS is weakest, and specifically the autotransporters, where
4 of 10 calls rest on a SignalP model trained on them.

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

`match_txsscan_organisms.py` and `deeplocpro_heldout_check.py` read the cached UniProt
records, so re-run `fetch_benchmark_uniprot.py` first whenever the list changes.
