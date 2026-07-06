## Context

`--tier full` enables `blastp` (vs NR) and `hhsuite` (vs UniRef30) on top of the extended tier (see `constants.py` `_TIER_ENABLED`). Three independent gaps make a full run fail end-to-end today:

1. `PipelineConfig.__post_init__` (runner.py ~399-448) resolves every optional database from its env var + marker root, **except** `blastp_db`. The manifest defines `SSIGN_BLAST_NR` (subpath `blast_nr`, sentinel `nr.pdb`) but nothing reads it, so `blastp_db` stays empty and the step returns the definitive "BLASTp requires a local database" error (runner.py:2099).
2. `scripts/cx3/run_batched_multi.pbs` hardcodes `--tier extended` (line 89) and exports only pfam + pdb70, never UniRef30.
3. `scripts/cx3/submit_batched_overnight.sh` has no way to select a tier.

The extended tier is the everyday path and must stay byte-for-byte unchanged; this change only adds a full-tier route.

## Goals / Non-Goals

**Goals:**
- A full-tier run resolves NR + UniRef30 with no explicit `--*-db` flags, matching how bakta/eggnog/IPS/ECOD/hhsuite-pfam already resolve.
- A one-flag CX3 submit (`--tier full`) that exports the two extra databases.
- Preserve the exact extended-tier behaviour and the existing "no NR → clear error" path.

**Non-Goals:**
- Adding blastp/hhsuite to the runtime effort model (they keep fixed timeout caps; the walltime caveat is documented, not fixed here).
- Changing substrate-calling, enrichment, or figure logic.
- Fetching the databases (that is `fetch_databases.sh --tier full`, already implemented).
- Graceful HH-suite degradation when UniRef30 is absent (still out of scope, per the `_EXTENDED_ADDS` comment).

## Decisions

**D1 — Resolve `blastp_db` in the runner's Step-A env passthrough, and bridge directory→prefix.**
BLAST+ `-db` needs the database *name prefix* (`.../blast_nr/nr`), but the shared `DatabasePath.resolve_path` returns the *directory* that holds the sentinel (`.../blast_nr/`), and `ssign doctor` treats `SSIGN_BLAST_NR` as that directory. To keep one consistent meaning for the env var (a directory, doctor-compatible) while giving BLAST+ what it needs, the runner resolves the NR directory the normal way and then appends the `nr` basename to form `blastp_db`. `nr` is hardcoded because NR is the only database `fetch_databases.sh --tier full` installs under `blast_nr/`; a user wanting Swiss-Prot or a custom DB passes `--blastp-db` explicitly (which still wins).

- *Alternative A (rejected):* make `SSIGN_BLAST_NR` carry the prefix and pass it verbatim. Rejected: it diverges from doctor's directory semantics and from every other DB env var, and breaks doctor's sentinel glob (`$SSIGN_BLAST_NR/nr.pdb`).
- *Alternative B (rejected):* append `/nr` inside `run_blastp.py`. Rejected: scatters the NR-specific knowledge into the tool wrapper; the runner is where all DB resolution already lives.

**D2 — `--tier` is a submit flag threaded through a `TIER` env var, defaulting to `extended`.**
`submit_batched_overnight.sh` adds `--tier <base|extended|full>` and puts `TIER=<value>` in the PBS `-v` string; `run_batched_multi.pbs` reads `: "${TIER:=extended}"` and runs `ssign run --tier "$TIER"`. Defaulting to `extended` means every existing submit is unchanged, and the hardcoded `--tier extended` line is replaced by the variable.

- *Alternative (rejected):* keep hardcoded tier and pass `--tier full` via `SSIGN_EXTRA_ARGS`. Rejected: relies on argparse last-wins and packs space-separated paths into the PBS `-v` value, exactly the paste-fragility the cx3-submit skill warns against.

**D3 — PBS exports NR + UniRef30 conditionally, like the other DBs.**
Add two guarded exports next to the existing pfam/pdb70 block: `SSIGN_HHSUITE_UNICLUST` from `$DBROOT/hhsuite/uniref30` and `SSIGN_BLAST_NR` from `$DBROOT/blast_nr`, each only when the directory exists. `SSIGN_HHSUITE_UNICLUST` is already read by the runner; `SSIGN_BLAST_NR` becomes live via D1.

## Risks / Trade-offs

- **HH-suite walltime blowout on large panels** → HH-suite runs per-substrate (thread pool, cpu_per_job=2) with a 1h + 30min per-protein cap and is not size-modelled. A 74-genome panel (~1000-1200 substrates) is the tail. Mitigation: validate on a 2-4 genome smoke test first, measure real per-protein time, and size `--walltime 24h+`; document the caveat.
- **Hardcoded `nr` basename** → wrong if a future full-tier ships a differently-named BLAST DB. Mitigation: `--blastp-db` override always wins; the assumption is documented and matches the current `fetch_databases.sh` layout (sentinel `nr.pdb`).
- **`$EPHEMERAL` purge** → NR/UniRef30 can vanish after ~30 days. Mitigation: the exports are guarded on directory existence, and `ssign doctor --tier full` reports what survived before submitting.

## Migration Plan

1. Land the runner + script changes with unit tests (NR resolution: env-set, explicit-override, absent-error).
2. On CX3: `fetch_databases.sh --tier full` onto `$EPHEMERAL`; verify with `ssign doctor --tier full`.
3. Smoke test: `submit_batched_overnight.sh --tier full --tutorial-all` (or 2 genomes), confirm BLASTp + HH-suite complete and record wallclock.
4. Only then run a full-tier panel with `--walltime` sized from the smoke test.

Rollback: revert the change; extended-tier submits are unaffected at every step because the default tier stays `extended`.

## Open Questions

- None blocking. The first full-tier genome target (Xanthobacter panel vs benchmark corpus) is a run-time choice, not a design dependency.
