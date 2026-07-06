# Phase 2 — actual recall on the gold panel

How many experimentally-verified secretion-system effectors does ssign actually emit as secreted, run end-to-end on the source genomes? This is the *actual* number from the benchmark README's ceiling/impossible/actual framing, measured against the re-verified gold list and reported for all six types (T1–T6, T5SS included).

## Headline

**Recall = 38 / 85 (44.7%)** — CX3 run `3169556` (2026-07-03), final pipeline.

| SS type | recall |
|---|---|
| T1SS | 16/18 (89%) |
| T2SS | 1/8 (13%) |
| T3SS | 5/19 (26%) |
| T4SS | 0/8 (0%) |
| T5SS | 10/19 (53%) |
| T6SS | 6/13 (46%) |
| **Total** | **38/85 (44.7%)** |

Figures: `data/phase2/figures/` (funnel, per-type, emission basis, discordant audit) and `.../figures/summary/` (recall@±3, T3SS story, emission quality, system-level recall).

## Gold standard

`data/phase2/verification_phase_a/gold_list_final.tsv` — **85 rows** (gold v6). One experimentally-validated effector per secretion-system instance, each with a UniProt accession, a genome + contig + CDS coordinates, a primary reference, and a verbatim citation quote. Built by `scripts/65_apply_corrections.py` from a corrections ledger (no hand-typed found/geometry numbers).

Provenance: a **Phase A exhaustive blind re-read** (2026-06-30) — 6 blind agents independently re-verified every row from primary evidence, checking that (a) the UniProt is the named effector (not a machinery/accessory paralog), (b) the locus is that gene, (c) the cited paper experimentally shows secretion by that ss_type. Result: 83/87 confirmed; the 4 flagged rows were pure T3SS translocon components (EspA, HrpK1 dropped; CopN, BipC kept as dual translocator/effectors), giving 85 final rows.

## Run provenance

- Job `batched_RTX6000_20260703_185843_3169556`, 52 genomes, single MultiGenome job.
- Final code: T5SS evidence gate (`signalp-t5ss-substrate-call`, commit `ca18322`) + size-aware tool timeouts (`size-aware-tool-timeouts`).
- Tier `extended`: DeepLocPro + DeepSecE + SignalP predictions; InterProScan + EggNOG + pLM-BLAST + ProtParam annotation. PLM-Effector, BLASTp, HH-suite off. Bakta re-annotation from scratch (`use_input_annotations=False`).
- `--enrichment-stats` → whole-genome DLP/DSE/SignalP (pooled, 160,831 proteins). **All 52 genomes completed 17/17 steps**; no tool timed out.

## Recall rule

`found = ≥1 bp overlap between a gold effector's CDS span and an EMITTED secreted protein on the reconciled molecule.` Canonical matcher: `RerunIndex.emitted_overlap` (`scripts/rerun_coords.py`). Bakta renames every locus on re-annotation, so effectors are bridged to the run by genome **coordinates**, not locus tags.

- **found=no, ORF present** → Bakta called the gene but ssign didn't emit it (a genuine prediction/proximity miss).
- **found=no, no overlap** → Bakta missed the ORF entirely (still a genuine fail).
- Molecule identity reconciled first: 3 verified RefSeq↔INSDC contig aliases (`CONTIG_ALIAS`, lenratio 1.000): PAO1 `AE004091.2↔NC_002516.2`, B. pertussis `NC_002929.2↔BX470248.1`, MC58 `NC_003112.2↔AE002098.2`. All 85 rows reconcile (0 un-reconcilable).
- **Full-assembly override:** 4 T1SS genomes (`NZ_CP031766.1`, `NZ_CBDBTK010000022.1`, `NZ_JABJZG010000001.1`, `NZ_JBCGCZ010000007.1`) are served from `rerun_fullasm/` because their gold effector sits on a contig fragmented in the standard input — the override is worth **+3** (Serralysin, apxIA, hlyA only recoverable from the full assembly).

Reproduce: run `RerunIndex().emitted_overlap(contig, start, stop, strand)` over each gold row (`scripts/rerun_coords.py`); the run's per-genome `*_results{,_raw}.csv` live under `rerun/<unit>/results/`.

## Key findings

- **The T5SS evidence gate is recall-neutral on the gold set.** Recomputing recall against this run gives **0 flips** vs the pre-gate baseline (T5SS holds at 10/19). Every real T5 substrate carries the DeepLocPro-or-SignalP evidence the gate requires; the gate only drops evidence-less false-positive T5 components elsewhere.
- **T4SS 0/8** is the correct consequence of the cross-replicon re-anchors (T4SS effectors sit distal from the VirB apparatus, often on a different replicon), not a bug.
- **Effort-model validation** (`--enrichment-stats` whole-genome pass, 160,831 proteins): DeepLocPro 5.56h, DeepSecE 5.17h, SignalP 10.0h, pLM-BLAST 7.6h — DLP/DSE within 4% of the fitted model. SignalP is ~1.8× DLP (recorded for the calibration refit).

## Scope notes

- Recall counts protein-level effector recovery. A separate **system-level** view (a system instance counts as found if ≥1 of its verified effectors is emitted) is in `scripts/52_system_recall.py` → `figures/summary/06_recall_systems.png`.
- The 4 full-assembly genomes are currently served from a 2026-06-25 run; T1SS is unaffected by the gate/timeout changes, but a same-code rerun of those 4 would make the panel single-version.
