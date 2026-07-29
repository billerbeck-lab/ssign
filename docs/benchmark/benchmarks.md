# Benchmarking ssign: effector recovery

Does ssign find real secretion systems and secreted protein pairs? This page reports the
**effector-recovery benchmark**. A test run on real genomes to find out how many
experimentally-validated secretion-system and secreted proteins ssign actually finds.

By design, ssign only calls a substrate if a
secreted-looking protein sits within a few genes of detected secretion-system
machinery ([`pipeline_overview.md`](../explanation/pipeline_overview.md#phase-4-substrate-identification)).
This means a secreted protein encoded far from its own secretion system machinery is always unreachable.

## The benchmarking list

[`ssign_benchmarking_list.csv`](ssign_benchmarking_list.csv)
holds **85 experimentally-validated effectors** one per secretion-system
instance, spanning
T1–T6SS. Every row carries a UniProt accession, the source genome + contig + CDS
coordinates, and the primary reference (DOI). This was created via mining the literature.

The ceiling is the fraction of each type's listed effectors whose own machinery
sits within ±N genes.

| SS type | in list | ±1 | ±3 | ±5 | ±7 | ±9 |
|---|---:|---:|---:|---:|---:|---:|
| T1SS | 18 | 56% | 89% | 89% | 89% | 89% |
| T2SS | 8 | 13% | 13% | 13% | 13% | 13% |
| T3SS | 19 | 11% | 21% | 32% | 32% | 37% |
| T4SS | 8 | 0% | 0% | 0% | 0% | 0% |
| T5SS | 19 | 84% | 89% | 95% | 95% | 95% |
| T6SS | 13 | 54% | 69% | 77% | 77% | 77% |
| **All** | **85** | **42%** | **55%** | **60%** | **60%** | **61%** |

## Actual recall

Each type splits into how many secreted proteins proximity could
reach at all (**reachable**: within ±3 genes of their own machinery) and how many
ssign actually finds (**recovered**).

Run end-to-end on the 52 source genomes, ssign predicts **38 of the total 85 secreted proteins (44.7%) and 38 of the 48 reachable secreted proteins (79.2%)**.

| SS type | in list | reachable ±3 | recovered |
|---|---:|---:|---:|
| T1SS | 18 | 16 (89%) | 16 (89%) |
| T2SS | 8 | 1 (13%) | 1 (13%) |
| T3SS | 19 | 4 (21%) | 5 (26%)* |
| T4SS | 8 | 0 (0%) | 0 (0%) |
| T5SS | 19 | 17 (89%) | 10 (53%) |
| T6SS | 13 | 9 (69%) | 6 (46%) |
| **Total** | **85** | **47 (55%)** | **38 (45%)** |

*T3SS recovered (5) exceeds reachable (4) by one: VirA (Shigella) is recovered
through a neighbouring T5aSS autotransporter's window (icsA, one gene away), not
its own T3SS machinery (~24 genes away). It is the only cross-system case in the
list. This is a strange case of serendipity.

![Effector recovery by secretion-system type. Green is recovered; amber is reachable but not recovered; grey is unreachable (the effector's own machinery is more than three genes away).](effector_recovery_by_type.png)

An effector counts as recovered if its CDS span overlaps an emitted secreted
protein. Because Bakta renames every locus on re-annotation, effectors are bridged
to the run by genome coordinates, not locus tags.

## Scope

This test shows that ssign is indeed effective at identifying secretion systems and the proteins they secrete for the majority of those that are in the proximity window. It also shows that it cannot identify all types of systems well, especially T2SS and T4SS are out of scope for this tool, due to how far away the proteins they secrete sit from the system machinery.

It is entirely possible that some of the systems and secreted proteins being tested here were used for the training of the tools ssign employs. This is a known limitation of this test and not accounted for. However the test still shows that the identification of secretion systems _and_ the proteins they secrete is a novel use ssign succeeds at.

## Reproduce

The benchmarking list is [`docs/benchmark/ssign_benchmarking_list.csv`](ssign_benchmarking_list.csv).
The recall matcher tests each listed effector's coordinates against the run's emitted
secreted proteins (`emitted_overlap`), with three verified RefSeq↔INSDC contig
aliases reconciled first so all 85 rows resolve. The recall table and figure are
`reachable_within_3` and `found_by_ssign` tallied by `ss_type` over that same list.

---

*Robustness sweeps (proximity-window sensitivity, assembly fragmentation, gene
copy number) and compute-scaling numbers are being re-run on this benchmark's
52-genome panel for consistency and will be added here once complete.*
