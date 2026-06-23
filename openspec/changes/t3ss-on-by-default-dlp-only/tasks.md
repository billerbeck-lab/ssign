## 1. Flip the default (all entry points)

- [x] 1.1 Dropped `T3SS` from the default `excluded_systems` in `runner.py` (`PipelineConfig`), `cli.py` (flag default + help), `Home.py` (multiselect default + `excluded_systems=` fallback), `system_filtering.py`, and ALSO `validate_macsyfinder_systems.py` + the `DEFAULT_EXCLUDED_SYSTEMS` constant (both missed by the original plan; found via grep). New default: `Flagellum, Tad`.
- [x] 1.2 Grep confirms no remaining `Flagellum,Tad,T3SS` default in shipping code; only the cli.py help text ("pass 'Flagellum Tad T3SS' to restore") references the old triple intentionally.

## 2. DeepSecE never trusted for T3SS

- [x] 2.1 `cross_validate_predictions._dse_flag`: `t3ss_flagged = (dse_type == "T3SS")`, unconditional.
- [x] 2.2 Removed `has_t3ss` from `_dse_flag` + `cross_validate` signature + caller; deleted `_genome_has_t3ss` and the now-unused `--valid-systems` arg (updated the runner call site + 3 shell smoke tests that passed it). Module docstring rewritten to "always flagged".

## 3. Tests

- [x] 3.1 `test_dse_t3ss_flagged_even_alongside_real_secretion`: DSE T3SS flagged + uncounted even when another tool carries the protein (the case the old conditional let through).
- [x] 3.2 `test_dse_t3ss_always_flagged_and_excluded` + `test_non_t3ss_dse_calls_count` (T1SS DSE still counts). `_run` helper de-threaded.
- [x] 3.3 `test_excluded_systems_default` now asserts `["Flagellum","Tad"]`; added `test_excluded_systems_default_matches_constant` (PipelineConfig == `DEFAULT_EXCLUDED_SYSTEMS`, no T3SS). CLI parser is built inline in `main()` (not importable), so covered via the constant.
- [x] 3.4 Already covered by existing enrichment tests: `test_dse_t3ss_excluded` + per-type aggregation asserting `("T3SS","DSE")` is skipped while T3SS DLP/window is produced (logic pre-wired via `ENRICH_DSE_NO_WINDOW`).
- [x] 3.5 Full unit suite green: 1377 passed.

## 4. Docs

- [x] 4.1 CLAUDE.md: `excluded_systems` default + Key Parameters + reframed critical-bug-fix #4 (T3SS detected by default; DeepSecE unconditionally excluded for T3SS).
- [x] 4.2 README updated (default + DSE-off rationale + escape hatch). `docs/design_decisions.md` has no T3SS statement to update.

## 5. Validation (CX3)

- [x] 5.1 CX3 runs (1/4/20 genome, 2026-06-23) confirm it: Salmonella T3SS = 4 substrates (sane, not a flagellar blowup), `T3SS / DLP / window` row present, `T3SS / DSE` empty/n.a. across all genomes. Fleet T3SS counts 0-12/genome, 61 pooled over 20. All steps succeeded.
- [x] 5.2 Validation recorded in NOTES.md (CX3 RESULT 2026-06-23). Independent multi-genome enrichment-background bug found (null_mean ~7x low in batched runs) — logged separately in NOTES; does NOT affect T3SS detection/DSE-exclusion, single-genome enrichment, or substrate calls.
