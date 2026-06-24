## Context

T3SS was excluded by default because DeepSecE is an unreliable T3SS predictor: across a 74-genome fleet MacSyFinder validated 0 T3SS while DeepSecE predicted 1,808, the bulk being flagellar proteins misclassified as injectisome T3SS (critical-bug-fix #4). The blunt mitigation, putting `T3SS` in `excluded_systems`, suppresses the unreliable predictor but also discards genuine MacSyFinder/TXSScan T3SS detection, which is itself reliable.

`cross_validate_predictions.py` already contains a partial, conditional guard: `_dse_flag` sets `dse_T3SS_flagged = (dse_type == "T3SS") and not has_t3ss`, i.e. it only distrusts DeepSecE's T3SS calls in genomes that lack a validated T3SS. With T3SS excluded today, `_genome_has_t3ss` is effectively always False (the excluded T3SS never reaches `valid_systems`), so every DSE T3SS call is flagged, which is why the flagellar FPs have stayed suppressed.

The enrichment test (post enrichment-circular-shift-per-run) is already T3SS-aware: `ENRICH_WINDOW_TYPES` includes `T3SS` and `ENRICH_DSE_NO_WINDOW = {T3SS}`, so it would emit a DLP-window T3SS row and skip the DSE T3SS row, it just never sees a T3SS because the system is excluded upstream.

## Goals / Non-Goals

**Goals:**
- T3SS detected and substrate-called by default (default `excluded_systems = [Flagellum, Tad]`).
- DeepSecE never counts toward a T3SS call, in cross-validation and in enrichment, independent of genome content.
- No reintroduction of the flagellar T3SS false positives.

**Non-Goals:**
- Changing how DeepSecE is used for non-T3SS types.
- Changing SignalP or PLM-Effector behaviour (SignalP-negative is correct for T3SS effectors; PLM-E stays off by default).
- Improving T3SS effector recall beyond what DLP + proximity already provide. T3SS substrate calls are expected to be conservative (see Risks).

## Decisions

**1. Make the DSE-T3SS flag unconditional.** Change `_dse_flag` so `t3ss_flagged = (dse_type == "T3SS")` with no dependence on `has_t3ss`. Rationale: DeepSecE cannot distinguish injectisome from flagellar T3SS, so its T3SS calls are untrustworthy *whether or not* the genome has a validated T3SS, the presence of a real T3SS does not make DSE's per-protein T3SS labels reliable. Keeping the `has_t3ss` condition would, the moment T3SS detection is on, start trusting DSE T3SS calls genome-wide in T3SS-bearing genomes, reintroducing the exact flagellar FPs the exclusion was protecting against. *Alternative considered:* trust DSE T3SS only for proteins inside a validated T3SS neighbourhood, rejected as more complex and still exposed to flagellar genes that happen to fall near an injectisome locus.

**2. Retire `_genome_has_t3ss`/`has_t3ss` from the flagging path.** Once the flag is unconditional, `has_t3ss` is no longer needed to decide flagging. Remove it from `_dse_flag`'s signature and its call sites unless it is consumed elsewhere (a grep confirms its only consumer; remove the helper if it becomes dead).

**3. Flip the default in all five locations at once.** `excluded_systems` is duplicated as a literal default in `runner.py`, `cli.py`, two spots in `Home.py`, and `system_filtering.py`. All must change together or the effective default differs by entry point (CLI vs GUI vs internal). This is the cross-file-drift class of bug; grep guarantees completeness.

**4. Enrichment needs no code change.** It is already DLP-only for T3SS; the change is purely that T3SS components now exist to test. Verify by test rather than edit.

## Risks / Trade-offs

- **[Flagellar FPs return] →** Mitigated by Decision 1 (unconditional DSE-T3SS flag) plus a regression test asserting DSE T3SS is flagged in a genome that *has* a validated T3SS, the case the old conditional would have let through.
- **[Few/no T3SS substrates called] →** With DSE off, SignalP-negative (no Sec signal on T3SS effectors), and DLP-extracellular often weak for injected (not freely secreted) effectors, T3SS substrate calls will lean on DLP + proximity and may be sparse. This is acceptable and safe (conservative > flagellar blowup); the enrichment DLP row quantifies whether DLP carries any T3SS signal. Surfaced as an expected outcome, not a bug.
- **[Behaviour change for existing users] →** BREAKING default. Mitigated by the documented `--excluded-systems Flagellum,Tad,T3SS` escape hatch and a CHANGELOG/CLAUDE.md note. No external users yet.
- **[Stale docs] →** CLAUDE.md critical-bug-fix #4 currently says "T3SS excluded by default"; left unedited it becomes actively wrong. Reframed as part of this change.

## Migration Plan

1. Land code + docs + unit tests; full unit suite green.
2. Validate on CX3: re-run Salmonella LT2 (SPI-1/SPI-2) with `--enrichment-stats`; confirm (a) a non-absurd T3SS substrate count, (b) a `T3SS / DLP / window` enrichment row and no `T3SS / DSE` row, (c) no flagellar blowup.
3. Rollback: restore `T3SS` to the default `excluded_systems` (single-commit revert); no data migration involved.

## Open Questions

- None blocking. The "is DLP even meaningful for injected T3SS effectors?" question is empirical and answered by the Salmonella validation row, not a design decision.
