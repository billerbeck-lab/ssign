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

- [x] 5.1 On CX3: `fetch_databases.sh --tier full` onto `$EPHEMERAL` + Bakta-full; `ssign doctor --tier full` = 8/8 databases OK (Bakta db, Swiss-Prot, UniRef30 all present). 2026-07-08.
- [x] 5.2 Smoke test: `submit_batched_overnight.sh --tier full` (2 xantho genomes). **3 runs**: smoke-test-1 (3238415) found the HH-suite + BLASTp bugs; smoke-test-2 (3252429) validated those fixes, found the EggNOG OOM; smoke-test-3 (3256040, 2026-07-09) **PASSED clean** — HH-suite + BLASTp (Swiss-Prot) + EggNOG all `-> OK`, 16/16 both genomes. All four full-tier fixes confirmed.
- [ ] 5.3 Append the measured blastp/hhsuite wallclock to `memory/calibration/runs.jsonl` and size the panel `--walltime` from it. **Blocked on clean per-tool timings**: the smoke-test-3 `ssign.run.log` collapsed the live `elapsed` fields (terminal redraws in the tee), so per-tool HH-suite/BLASTp wallclocks need the PBS `.o` log or a re-timed run. Panel `--walltime 24:00:00` is a safe interim.

## 6. Post-smoke-test robustness fixes (defects the first full-tier smoke test surfaced)

The first full-tier smoke test (run 3238415) exposed three real defects a new-lab / non-CX3
user would also hit. NR was additionally replaced by Swiss-Prot as the full-tier BLASTp default
(curated, minutes vs NR's hours on a substrate set), so tasks 1.x/2.3/3.x now target
`SSIGN_BLAST_SWISSPROT` + `<dir>/swissprot`, not NR.

- [x] 6.1 HH-suite DB dir→ffindex-prefix: `run_hhsuite.py` `_resolve_ffindex_prefix` converts the resolved DB dir to the ffindex prefix hhblits/hhsearch expect (was IsADirectoryError). Commit 98b2f5f.
- [x] 6.2 BLASTp default NR→Swiss-Prot: fetch/manifest/runner/PBS/docs target Swiss-Prot; NR opt-in via `--blastp-db`. Commits 74829d6, e62bc60.
- [x] 6.3 HH-suite not on PATH in PBS runs: `run_batched_multi.pbs` + `run_k12_validation.pbs` relied on the submitter's `~/.bashrc`, which does NOT export `hhblits`. Both now put HH-suite on PATH self-containedly (keep any existing copy, else auto-detect the conda env shipping `hhblits`).
- [x] 6.4 Bakta full/light `db*` glob collision (fetch skip + runner resolve): `fetch_databases.sh` keys the skip-guard on the exact variant subdir (`db/` full vs `db-light/` light) so extended→full re-downloads; `dependency_manifest.py` resolve_path returns `os.path.dirname(max(matches))` so `db/` (full) deterministically beats `db-light/`. Commit 99fac0b.
- [x] 6.5 EggNOG DIAMOND OOM in the parallel annotation group (surfaced by full-tier smoke-test-2, run 3252429): EggNOG autodetected `--dbmem`/`--block_size`/`--index_chunks` off whole-node RAM while sharing the node with 5 other tools → overcommit → DIAMOND OOM-killed. New `parallel_share_ram_gb()` (RAM analogue of `parallel_share_cpus()`) divides by `SSIGN_PARALLEL_GROUP_SIZE`; the three autodetect fns use it. Commit 1e20f66.
