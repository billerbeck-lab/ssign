# CLI reference

Complete flag list for the `ssign` command. For "how do I X" recipes, see
[`how-to/run.md`](../how-to/run.md).

`ssign` is a command-line tool with three subcommands (running it with no
subcommand prints a short usage hint):

```bash
ssign run input.gbff --outdir <dir>  # run the pipeline
ssign doctor --tier <tier>           # verify the install
ssign fetch-databases --tier <tier>  # download reference databases
```

Most flags below are for `ssign run`; `doctor` and `fetch-databases` are
covered under [Other subcommands](#other-subcommands).

Most boolean flags use `argparse.BooleanOptionalAction`, so each `--<flag>`
accepts a `--no-<flag>` inverse (e.g. `--skip-blastp` and `--no-skip-blastp`).
The exceptions are `--resume`, `--skip-localization-gate`, `--imports-only`,
`--dry-run`, and `--version`: these are plain on/off switches with no `--no-`
form.

## Top-level flags

| Flag | Type | Default | Description |
|---|---|---|---|
| `--version` | bool | `false` | Print the ssign version and exit. |

## `ssign run` essentials

| Flag | Type | Default | Description |
|---|---|---|---|
| `INPUT_PATH` | path (positional) | required | Path to the input genome (GenBank `.gbff`/`.gbk`, GFF3 `.gff`, or FASTA). |
| `--outdir` | path | `./results` | Output directory. |
| `--sample-id` | str | input stem | Prefix for output filenames. Defaults to the input filename's stem. |
| `--original-filename` | str | `""` | Original filename when `INPUT_PATH` is a temp upload (informational). |
| `--resume` | bool | `false` | Skip steps that already succeeded in a previous run (reads `<outdir>/.ssign/<sample-id>_progress.json`). |
| `--tier` | choice | unset → `extended` | Install tier the run targets (`base`/`extended`/`full`); sets each tool's default on/off state to what that tier ships. Unset uses what `fetch-databases` recorded, else `extended`. |
| `--combined-summary` | bool | `true` | (multi-genome) Write a top-level `combined_results.csv` aggregating every genome's substrates with a `source_genome` column. |

## SS detection (MacSyFinder)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--wholeness-threshold` | float | `0.8` | Minimum MacSyFinder system completeness. |
| `--excluded-systems` | list | `Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P` | Space-separated SS types to exclude from substrate calling. T3SS is NOT excluded by default. |
| `--macsyfinder-db-type` | choice | `ordered_replicon` | MacSyFinder `--db-type`. Choices: `ordered_replicon`, `unordered`. |
| `--cpu-per-genome` | int | CPU count | CPUs available to per-genome subtools (passed as `-w` to MacSyFinder, `-num_threads` to BLAST, etc.). |

## Prediction thresholds

| Flag | Type | Default | Description |
|---|---|---|---|
| `--conf-threshold` | float | `0.8` | DeepLocPro extracellular probability minimum. |
| `--proximity-window` | int | `3` | +/-N genes per SS component for proximity neighbourhood. |
| `--required-fraction-correct` | float | `0.8` | Fraction of SS components that must be correctly localized for the system to pass. |
| `--deepsece-min-prob` | float | `0.8` | DeepSecE minimum probability to call a protein secreted. |
| `--signalp-min-prob` | float | `0.5` | SignalP minimum probability for a signal peptide. |
| `--dlp-confidence-threshold` | float | `0.8` | Minimum DeepLocPro max-probability for an SS-machinery component to count in the localization-correctness gate. Components below this are excluded from both sides of `fraction_correct`. Distinct from `--conf-threshold` (which gates extracellular calls). |
| `--skip-localization-gate` | bool | `false` | Disable the literature-derived localization-correctness gate (debug / ad-hoc). |

## Enrichment analysis

| Flag | Type | Default | Description |
|---|---|---|---|
| `--enrichment-stats` | bool | `false` | Per-SS-type circular-shift enrichment test: emits fold (enrichment) + permutation p + BH q per system type for DeepLocPro, DeepSecE, and SignalP, plus enrichment figures. Forces whole-genome DeepLocPro + DeepSecE + SignalP (local; the rotation null needs every gene's positivity in gene order), ~13 min/genome. Pool across genomes for statistical power. |

## ORF prediction and annotation

| Flag | Type | Default | Description |
|---|---|---|---|
| `--use-input-annotations` | bool | `false` | Trust input GenBank annotations and skip Bakta re-annotation. |
| `--run-bakta` | bool | `true` | Run Bakta on FASTA contigs input (default on). GenBank input is governed by `--use-input-annotations` instead. |
| `--bakta-db` | path | `""` | Bakta database directory (required when `--run-bakta`). |
| `--bakta-threads` | int | `0` | Threads passed to Bakta. `0` = match `--cpu-per-genome` (the cgroup-allocated count). |

## DTU prediction tools (DeepLocPro and SignalP)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--deeplocpro-mode` | choice | auto (local when a local install is detected) | `local` (canonical, DTU academic licence required) or `remote` (opt-in fallback: DTU webserver, no licence needed but depends on DTU hosting the service). Unset = auto: local if `deeplocpro` is on `PATH`/`--deeplocpro-path`/`$SSIGN_DEEPLOCPRO_PATH`, otherwise ssign stops with install instructions (it does not auto-submit to the webserver). |
| `--deeplocpro-path` | path | `""` | Path to local DeepLocPro install. Empty falls back to `deeplocpro` on `PATH`. |
| `--signalp-mode` | choice | auto (local when a local install is detected) | `local` (canonical; obtain SignalP 6.0 from the DTU portal, DTU does not redistribute it) or `remote` (opt-in fallback: DTU webserver, no licence needed but depends on DTU hosting the service). Unset = auto: local if `signalp6` is on `PATH`/`--signalp-path`/`$SSIGN_SIGNALP_PATH`, otherwise ssign stops with install instructions (it does not auto-submit to the webserver). |
| `--signalp-path` | path | `""` | Path to local SignalP 6 install. Empty falls back to `signalp6` on `PATH`. |
| `--skip-deeplocpro` | bool | `false` | Skip the DeepLocPro step (overrides the `--tier` default). |
| `--skip-signalp` | bool | `false` | Skip the SignalP step. |
| `--skip-deepsece` | bool | `false` | Skip the DeepSecE step. |
| `--dlp-whole-genome` | bool | `false` | Run DeepLocPro on every protein, not just the SS neighbourhood. |
| `--dse-whole-genome` | bool | `false` | Run DeepSecE on every protein, not just the SS neighbourhood. |
| `--sp-whole-genome` | bool | `false` | Run SignalP on every protein, not just the SS neighbourhood. |

## BLASTp

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-blastp` | bool | `true` (off), `false` at `--tier full` | Skip BLASTp. BLASTp is enabled only at `--tier full`; off at base/extended (the default tier is extended). |
| `--blastp-db` | path | `""` | Path to BLAST database (NR or Swiss-Prot). |
| `--blastp-exclude-taxid` | str | `""` | Comma-separated taxids to exclude (e.g. the query organism). |
| `--blastp-min-pident` | float | `80.0` | BLASTp percent-identity floor. |
| `--blastp-min-qcov` | float | `80.0` | BLASTp query-coverage floor. |
| `--blastp-evalue` | float | `1e-5` | BLASTp e-value threshold. |

## HH-suite

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-hhsuite` | bool | `true` (off), `false` at `--tier full` | Skip HH-suite. HH-suite is enabled only at `--tier full`; off at base/extended (the default tier is extended). Its hhblits step needs UniRef30, which only ships with the full tier. |
| `--hhsuite-pfam-db` | path | `""` | HH-suite Pfam database. Falls back to `$SSIGN_HHSUITE_PFAM`. |
| `--hhsuite-pdb70-db` | path | `""` | HH-suite PDB70 database. Falls back to `$SSIGN_HHSUITE_PDB70`. |
| `--hhsuite-uniclust-db` | path | `""` | UniClust / UniRef30 database. Falls back to `$SSIGN_HHSUITE_UNICLUST`. |
| `--hhsuite-min-prob` | float | (constants) | HH-suite probability floor. Defaults to `ssign_lib.constants.HHSUITE_MIN_PROB`. |

## InterProScan

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-interproscan` | bool | `false` (on), `true` at `--tier base` | Skip InterProScan. InterProScan is on at the extended (default) and full tiers; off only at `--tier base`. |
| `--interproscan-db` | path | `""` | InterProScan install directory. |
| `--interproscan-min-evalue` | float | `1e-5` | InterProScan e-value threshold. |

## pLM-BLAST

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-plmblast` | bool | `false` | Skip pLM-BLAST. On by default at the extended tier (ECOD30, ~40 min/genome). |
| `--plmblast-db` | path | `""` | Path to ECOD pLM-BLAST database (ECOD30 default; ECOD50/70/90 also supported). |
| `--plmblast-cpc` | int | `90` | pLM-BLAST cosine percentile cutoff (Kamiński 2023 default). Drop to 70-80 for more permissive matching on short proteins (<200 aa), at the cost of longer search wallclock. |

## EggNOG-mapper

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-eggnog` | bool | `false` (on), `true` at `--tier base` | Skip EggNOG-mapper. EggNOG is on at the extended (default) and full tiers; off only at `--tier base`. Pass `--skip-eggnog` to disable it. |
| `--eggnog-db` | path | `""` | EggNOG database directory. |
| `--eggnog-dbmem` | bool | auto | Load `eggnog.db` into RAM (`--dbmem`, ~44 GB resident). Default auto: on only when the job's RAM share is >= 50 GB, else the on-disk SQLite is memory-mapped. Force with `--eggnog-dbmem` / `--no-eggnog-dbmem`. |

## Miscellaneous annotation

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-annotation` | bool | `false` | Skip every annotation tool at once (BLASTp, HH-suite, InterProScan, pLM-BLAST, EggNOG, ProtParam). Predictions and substrate calls still run. A per-tool `--no-skip-<tool>` (e.g. `--no-skip-eggnog`) overrides this to keep that one tool on. |
| `--skip-protparam` | bool | `false` | Skip the ProtParam physicochemical-property step. |
| `--t5ass-annotate-whole` | bool | `false` | For T5aSS (classical autotransporter) substrates, run EggNOG / BLASTp / pLM-BLAST / HHsuite / ProtParam a second time on the full protein and emit `t5ass_whole_*` columns alongside the default passenger-only annotations. Lets you compare functional (passenger) vs structural (β-barrel-dominated whole-AT) hits side by side. InterProScan unchanged (already domain-aware). See `docs/explanation/design_decisions.md` § 4.3. |
| `--filter-dse-type-mismatch` | bool | `true` | Drop DSE-only substrates whose predicted SS type does not match the nearby MacSyFinder system. |
| `--ortholog-min-pident` | float | `40.0` | Ortholog grouping percent-identity floor. |
| `--ortholog-min-qcov` | float | `70.0` | Ortholog grouping query-coverage floor. |

## Runtime & diagnostics

| Flag | Type | Default | Description |
|---|---|---|---|
| `--scratch-dir` | path | `""` (auto) | Directory for scratch/temp files (tool working dirs). Auto-resolves: keep `$TMPDIR` if it has adequate free space, else a dir under `--outdir`. Set this when running in a container whose `/tmp` is a small tmpfs (avoids Bakta "No space left on device"). |
| `--monitor-resources` | bool | `true` | Write `<outdir>/runtime_data/{step_timings,resource_samples}.csv` during a run. |
| `--monitor-interval-s` | float | `5.0` | Sampling interval for `resource_samples.csv` (seconds). |

## Figures

| Flag | Type | Default | Description |
|---|---|---|---|
| `--dpi` | int | `300` | Figure resolution. |
| `--fig-ss-comp` | bool | `true` | Render figure `01` (secreted proteins by SS type; per genome, stacked by SS type for a group). |
| `--fig-physicochemical` | bool | `true` | Render figure `02` (size & physicochemical properties by SS type; length + ProtParam when present). |
| `--fig-func-summary` | bool | `true` | Render figures `03`-`06` (functional categories by SS type: COG/KEGG/EggNOG/consensus). |

## Other subcommands

Besides `ssign run`, the CLI has two install helpers.

### `ssign doctor`

Verify the install: Python packages, external binaries, databases, model weights.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--tier` | choice | `extended` | Install tier to verify against (`base`/`extended`/`full`). |
| `--imports-only` | bool | `false` | Only check Python imports; skip binaries / DBs / weights (used by CI). |
| `--data-root` | path | `~/.ssign` | Root for databases + models. `SSIGN_*` env vars override per-DB paths. |

### `ssign fetch-databases`

Download the reference databases for a tier (wraps `scripts/fetch_databases.sh`).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--tier` | choice | required (must be passed) | Which tier's databases to download (base ~4 GB, extended ~100 GB, full ~500 GB). |
| `--target` | path | `~/.ssign/databases` | Destination directory. |
| `--dry-run` | bool | `false` | Print what would be downloaded without downloading anything. |
