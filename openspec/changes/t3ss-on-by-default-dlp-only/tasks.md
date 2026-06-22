## 1. Flip the default (all entry points)

- [ ] 1.1 Drop `T3SS` from the default `excluded_systems` in `src/ssign_app/core/runner.py` (`PipelineConfig`), `src/ssign_app/cli.py` (flag `default=` + help text), `src/ssign_app/Home.py` (both the multiselect default and the `excluded_systems=` fallback), and `src/ssign_app/scripts/system_filtering.py` (CLI `--excluded-systems` default). New default: `Flagellum, Tad`.
- [ ] 1.2 Grep the repo for the literal `Flagellum.*Tad.*T3SS` / `"Flagellum,Tad,T3SS"` to confirm no other entry point keeps the old default (cross-file-drift guard).

## 2. DeepSecE never trusted for T3SS

- [ ] 2.1 In `src/ssign_app/scripts/cross_validate_predictions.py` `_dse_flag`, make the flag unconditional: `t3ss_flagged = (dse_type == "T3SS")` (drop the `and not has_t3ss`).
- [ ] 2.2 Remove `has_t3ss` from `_dse_flag` (and its caller threading) once unconditional; grep for `_genome_has_t3ss` / `has_t3ss` and delete `_genome_has_t3ss` if it is now unused. Update the module docstring (the "DSE T3SS reliability guard" note) to say DSE T3SS is always flagged.

## 3. Tests

- [ ] 3.1 cross_validate: assert `dse_T3SS_flagged=True` and DSE excluded from the evidence count when DeepSecE predicts T3SS in a genome **with** a validated T3SS (the case the old conditional let through).
- [ ] 3.2 cross_validate: same assertion in a genome **without** a validated T3SS (unchanged behaviour); and a non-T3SS DSE type still counts.
- [ ] 3.3 config/CLI: assert the default `excluded_systems` is `["Flagellum", "Tad"]` (no T3SS) across `PipelineConfig`, the CLI parser, and `system_filtering` default.
- [ ] 3.4 enrichment: on a small fixture genome with a validated T3SS, assert the stats contain a `T3SS` row with `tool=DLP, mode=window` and no `T3SS` row with `tool=DSE`.
- [ ] 3.5 Run the full unit suite (`pytest tests/unit/ -q`); all green.

## 4. Docs

- [ ] 4.1 `CLAUDE.md`: update the `excluded_systems` default and Key Parameters; reframe critical-bug-fix #4 from "T3SS excluded by default" to "T3SS detected by default; DeepSecE excluded for T3SS (DLP-only)".
- [ ] 4.2 `README` and `docs/design_decisions.md`: update any statement that T3SS is excluded by default; document the `--excluded-systems Flagellum,Tad,T3SS` escape hatch and the DSE-off-for-T3SS rationale.

## 5. Validation (CX3)

- [ ] 5.1 Re-run Salmonella LT2 (SPI-1/SPI-2) on CX3 with `--enrichment-stats`; confirm a sane T3SS substrate count (no 1,808-style flagellar blowup), a `T3SS / DLP / window` enrichment row, and no `T3SS / DSE` row.
- [ ] 5.2 Record the validation outcome in NOTES and mark task #68 done.
