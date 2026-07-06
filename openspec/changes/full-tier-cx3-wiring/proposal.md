## Why

A `--tier full` run adds two annotation tools, BLASTp-vs-NR and HH-suite-vs-UniRef30, on the substrate set. Today a full run silently fails: the runner never reads the NR database env var, so BLASTp aborts with "BLASTp requires a local database", and the CX3 batch wrapper hardcodes `--tier extended` and never exports UniRef30, so HH-suite's MSA step cannot run. The full tier is unusable end-to-end until these are closed, which blocks the planned NR/HH-suite characterization of discovered substrates.

## What Changes

- **Runner resolves the BLAST NR database from its environment variable** like every other database. Currently `blastp_db` is the one database the runner's resolution block omits, so the manifest's `SSIGN_BLAST_NR` entry is dead. After this change, when `blastp_db` is unset and NR is discoverable (via `SSIGN_BLAST_NR` or the fetch-tool marker root), the runner resolves the `-db` prefix automatically. BLAST+ needs the DB name prefix (`.../blast_nr/nr`) whereas the resolver returns the containing directory (`.../blast_nr/`); the runner bridges directory to prefix by appending the `nr` basename. `--blastp-db` still overrides.
- **CX3 batch submit gains a `--tier` flag.** `scripts/cx3/submit_batched_overnight.sh` passes a `TIER` env var through to the PBS script, which runs `ssign run --tier "$TIER"` (default `extended`, so existing extended submits are unchanged).
- **The PBS wrapper exports the full-tier database paths.** `scripts/cx3/run_batched_multi.pbs` exports `SSIGN_HHSUITE_UNICLUST` (UniRef30) and `SSIGN_BLAST_NR` (NR) from `$EPHEMERAL/ssign-databases` when present, alongside the pfam/pdb70 exports it already sets, so HH-suite and BLASTp find their databases with no per-run flags.
- **Docs note the full-tier runtime caveat.** BLASTp and HH-suite are not in the runtime effort model, so size-aware timeouts do not scale them (fixed caps: BLASTp 4h floor; HH-suite 1h + 30min per protein). Full-tier panels need `--walltime 24h+`.

No change to substrate-calling, enrichment, or figure logic.

## Capabilities

### New Capabilities
- `full-tier-db-resolution`: A full-tier run resolves its BLAST NR and HH-suite UniRef30 databases automatically from the standard env vars / fetch marker root, so BLASTp and HH-suite run without the operator passing explicit `--blastp-db` / `--hhsuite-uniclust-db` flags.

### Modified Capabilities
<!-- None: no existing spec's requirements change. -->

## Impact

- `src/ssign_app/core/runner.py` — DB-resolution block (`PipelineConfig.__post_init__`, ~lines 399-448): add `blastp_db` resolution from `SSIGN_BLAST_NR` with dir→prefix handling.
- `scripts/cx3/submit_batched_overnight.sh` — add `--tier` flag, thread `TIER` into the `-v` string.
- `scripts/cx3/run_batched_multi.pbs` — accept `TIER` (default `extended`); export `SSIGN_HHSUITE_UNICLUST` + `SSIGN_BLAST_NR`.
- `docs/` (install / CX3 how-to) — full-tier database prerequisites + walltime caveat.
- Depends on `fetch_databases.sh --tier full` having placed NR (`blast_nr/`) and UniRef30 (`hhsuite/uniref30/`) on `$EPHEMERAL`. No new Python dependencies.
- Validation: a 2-4 tutorial-genome smoke test at `--tier full` on CX3 confirms the wiring and measures real BLASTp/HH-suite wallclock before any large panel.
