## 1. Runner: resolve BLAST NR (dir→prefix)

- [x] 1.1 In `PipelineConfig.__post_init__` (runner.py DB-resolution block), resolve `blastp_db` when unset: locate the NR directory via `SSIGN_BLAST_NR` + marker root (the `SSIGN_BLAST_NR` manifest entry / `resolve_path`), then set `blastp_db = <dir>/nr`. Do NOT override an explicit `blastp_db`.
- [x] 1.2 Log the resolution (`Resolved blastp_db -> <path>`) consistent with the other DB resolution log lines.
- [x] 1.3 Confirm the existing "BLASTp requires a local database" error still fires when NR is genuinely absent (no regression to runner.py:2099).

## 2. CX3 scripts: full-tier wiring

- [x] 2.1 `submit_batched_overnight.sh`: add `--tier <base|extended|full>` flag (default `extended`); thread `TIER=<value>` into the `-v` string and echo it in the summary.
- [x] 2.2 `run_batched_multi.pbs`: read `: "${TIER:=extended}"`; replace the hardcoded `--tier extended` with `--tier "$TIER"`.
- [x] 2.3 `run_batched_multi.pbs`: add guarded exports next to pfam/pdb70 — `SSIGN_HHSUITE_UNICLUST` from `$DBROOT/hhsuite/uniref30` and `SSIGN_BLAST_NR` from `$DBROOT/blast_nr` (each only when the dir exists); echo them in the DB summary line.
- [x] 2.4 `--dry-run` a `--tier full` submit locally; confirm `TIER=full` appears in the `-v` string and the extended default is unchanged when `--tier` is omitted.

## 3. Tests

- [x] 3.1 Unit test: `blastp_db` unset + `SSIGN_BLAST_NR` pointing at a temp dir containing an `nr.pdb` file → `config.blastp_db` resolves to `<dir>/nr`.
- [x] 3.2 Unit test: explicit `blastp_db` (or `--blastp-db`) is preserved even when `SSIGN_BLAST_NR` is set.
- [x] 3.3 Unit test: `blastp_db` unset + no `SSIGN_BLAST_NR` + no marker NR → `blastp_db` stays empty (error path intact).
- [x] 3.4 Run `pytest tests/unit/ -q`; suite green.

## 4. Docs

- [x] 4.1 Update the install / CX3 how-to (and the cx3-submit skill if needed) with the `--tier full` submit, the NR + UniRef30 prerequisites, and the `--walltime 24h+` runtime caveat (blastp/hhsuite not size-modelled).
- [x] 4.2 Note the full-tier database env vars (`SSIGN_BLAST_NR`, `SSIGN_HHSUITE_UNICLUST`) in the DB-resolution docs so a non-CX3 operator can wire them by hand.

## 5. CX3 validation (user-driven, gates the large panel)

- [ ] 5.1 On CX3: `fetch_databases.sh --tier full` onto `$EPHEMERAL`; `ssign doctor --tier full` reports NR + UniRef30 present.
- [ ] 5.2 Smoke test: `submit_batched_overnight.sh --tier full --tutorial-all` (or 2 genomes); confirm BLASTp + HH-suite steps complete (`N/N steps succeeded`) and record real per-tool wallclock.
- [ ] 5.3 Append the measured blastp/hhsuite wallclock to `memory/calibration/runs.jsonl` and size the panel `--walltime` from it.
