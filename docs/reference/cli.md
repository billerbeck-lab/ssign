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
| `--excluded-systems` | list | `Flagellum Tad T3SS` | Space-separated SS types to exclude. |
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
| `--deeplocpro-mode` | choice | `remote` | `remote` (DTU webserver, no licence needed) or `local` (DTU academic licence required). |
| `--deeplocpro-path` | path | `""` | Path to local DeepLocPro install (required when `--deeplocpro-mode local`). |
| `--signalp-mode` | choice | `remote` | `remote` (DTU webserver) or `local`. DTU does not redistribute SignalP 6.0; users obtain it from the DTU portal. |
| `--signalp-path` | path | `""` | Path to local SignalP 6 install (required when `--signalp-mode local`). |
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
| `--skip-plmblast` | bool | `true` | Skip pLM-BLAST (off by default). |
| `--plmblast-db` | path | `""` | Path to ECOD70 pLM-BLAST database. |

## EggNOG-mapper

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-eggnog` | bool | `true` | Skip EggNOG-mapper (off by default). |
| `--eggnog-db` | path | `""` | EggNOG database directory. |

## PLM-Effector

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-plm-effector` | bool | `true` | Skip PLM-Effector (off by default; GPU-heavy). |
| `--plm-effector-weights-dir` | path | `""` | Directory containing PLM-Effector weights and the ProtT5 cache. |
| `--plm-effector-types` | list | `T1SE T2SE T3SE T4SE T6SE` | Secretion-system effector types to predict. |

## Miscellaneous annotation

| Flag | Type | Default | Description |
|---|---|---|---|
| `--skip-protparam` | bool | `false` | Skip the ProtParam physicochemical-property step. |
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
