# CLI reference

Complete flag list for the `ssign` command. For "how do I X" recipes, see
[`how-to/configure.md`](../how-to/configure.md).

`ssign` has two modes:

```bash
ssign                                # launch the Streamlit GUI
ssign run input.gbff --outdir <dir>  # run the pipeline non-interactively
```

Boolean flags use `argparse.BooleanOptionalAction`, so each `--<flag>`
accepts a `--no-<flag>` inverse (e.g. `--skip-blastp` and `--no-skip-blastp`).

## Top-level flags

| Flag | Type | Default | Description |
|---|---|---|---|
| `--version` | bool | `false` | Print the ssign version and exit. |
| `--no-browser` | bool | `false` | (GUI mode) Start the server without opening a browser. |
| `--port` | int | `8501` | (GUI mode) Port for the Streamlit server. Falls forward to the next free port if `8501` is in use. |

## `ssign run` essentials

| Flag | Type | Default | Description |
|---|---|---|---|
| `INPUT_PATH` | path (positional) | required | Path to the input genome (GenBank `.gbff`/`.gbk`, GFF3 `.gff`, or FASTA). |
| `--outdir` | path | `./results` | Output directory. |
| `--sample-id` | str | input stem | Prefix for output filenames. Defaults to the input filename's stem. |
| `--original-filename` | str | `""` | Original filename when `INPUT_PATH` is a temp upload (informational). |
| `--resume` | bool | `false` | Skip steps that already succeeded in a previous run (reads `<outdir>/.ssign/<sample-id>_progress.json`). |

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

## ORF prediction and annotation

| Flag | Type | Default | Description |
|---|---|---|---|
| `--use-input-annotations` | bool | `false` | Trust input GenBank annotations and skip Bakta re-annotation. |
| `--run-bakta` | bool | `false` | Run Bakta on FASTA input or to re-annotate GenBank. |
| `--bakta-db` | path | `""` | Bakta database directory (required when `--run-bakta`). |
| `--bakta-threads` | int | `4` | Threads passed to Bakta. |

## DTU prediction tools (DeepLocPro and SignalP)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--deeplocpro-mode` | choice | `local` | `local` (canonical, DTU academic licence required) or `remote` (opt-in fallback: DTU webserver, no licence needed but depends on DTU hosting the service). |
| `--deeplocpro-path` | path | `""` | Path to local DeepLocPro install. Empty falls back to `deeplocpro` on `PATH`. |
| `--signalp-mode` | choice | `local` | `local` (canonical; obtain SignalP 6.0 from the DTU portal — DTU does not redistribute it) or `remote` (opt-in fallback: DTU webserver, no licence needed but depends on DTU hosting the service). |
| `--signalp-path` | path | `""` | Path to local SignalP 6 install. Empty falls back to `signalp6` on `PATH`. |
| `--skip-signalp` | bool | `false` | Skip the SignalP step. |
| `--skip-deepsece` | bool | `false` | Skip the DeepSecE step. |
| `--dlp-whole-genome` | bool | `false` | Run DeepLocPro on every protein, not just the SS neighbourhood. |
| `--dse-whole-genome` | bool | `false` | Run DeepSecE on every protein, not just the SS neighbourhood. |
| `--sp-whole-genome` | bool | `false` | Run SignalP on every protein, not just the SS neighbourhood. |

## BLASTp

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-blastp` | bool | `false` | Skip BLASTp. |
| `--blastp-db` | path | `""` | Path to BLAST database (NR or Swiss-Prot). |
| `--blastp-exclude-taxid` | str | `""` | Comma-separated taxids to exclude (e.g. the query organism). |
| `--blastp-min-pident` | float | `80.0` | BLASTp percent-identity floor. |
| `--blastp-min-qcov` | float | `80.0` | BLASTp query-coverage floor. |
| `--blastp-evalue` | float | `1e-5` | BLASTp e-value threshold. |

## HH-suite

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-hhsuite` | bool | `true` | Skip HH-suite (off by default; needs large databases). |
| `--hhsuite-pfam-db` | path | `""` | HH-suite Pfam database. Falls back to `$SSIGN_HHSUITE_PFAM`. |
| `--hhsuite-pdb70-db` | path | `""` | HH-suite PDB70 database. Falls back to `$SSIGN_HHSUITE_PDB70`. |
| `--hhsuite-uniclust-db` | path | `""` | UniClust / UniRef30 database. Falls back to `$SSIGN_HHSUITE_UNICLUST`. |
| `--hhsuite-min-prob` | float | (constants) | HH-suite probability floor. Defaults to `ssign_lib.constants.HHSUITE_MIN_PROB`. |

## InterProScan

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-interproscan` | bool | `false` | Skip InterProScan. |
| `--interproscan-db` | path | `""` | InterProScan install directory. |
| `--interproscan-min-evalue` | float | `1e-5` | InterProScan e-value threshold. |

## pLM-BLAST

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-plmblast` | bool | `false` | Skip pLM-BLAST. On by default at the extended tier (ECOD30, ~40 min/genome). |
| `--plmblast-db` | path | `""` | Path to ECOD pLM-BLAST database (ECOD30 default; ECOD50/70/90 also supported). |

## EggNOG-mapper

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-eggnog` | bool | `true` | Skip EggNOG-mapper (off by default). |
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

## Figures

| Flag | Type | Default | Description |
|---|---|---|---|
| `--dpi` | int | `300` | Figure resolution. |
| `--fig-category` | bool | `true` | Render the functional-category figure. |
| `--fig-ss-comp` | bool | `true` | Render the SS-component composition figure. |
| `--fig-tool-heatmap` | bool | `true` | Render the tool-coverage heatmap. |
| `--fig-substrate-count` | bool | `true` | Render the per-SS substrate-count figure. |
| `--fig-func-summary` | bool | `true` | Render the functional-summary figure. |
