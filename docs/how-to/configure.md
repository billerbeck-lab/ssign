# Configuring ssign

Recipes for the most common "how do I X" questions. For the full flag
reference, see [`reference/cli.md`](../reference/cli.md). For environment
variables, see [`reference/env_vars.md`](../reference/env_vars.md).

## Launching ssign

### GUI

```bash
pip install ssign
ssign                        # opens the Streamlit GUI in your browser
ssign --port 8502            # custom port
ssign --no-browser           # do not auto-open the browser
```

### CLI

```bash
ssign run input.gbff --outdir results
ssign run input.gbff --outdir results --skip-hhsuite --skip-plmblast
ssign run input.gbff --outdir results --resume   # pick up after a failed step
```

The CLI is the right interface for HPC and batch use. The GUI is for
single-genome interactive runs on a desktop.

## Input formats

| Format        | Extensions              | Notes |
|---------------|-------------------------|-------|
| GenBank       | `.gbff`, `.gbk`, `.gb`  | Recommended. Carries protein sequences and gene order; ssign re-annotates with Bakta by default unless `--use-input-annotations` is set. |
| FASTA contigs | `.fasta`, `.fna`, `.fa` | Bakta (or pyrodigal as fallback) predicts ORFs. |

GFF3 input requires a paired FASTA and is not yet wired through the
top-level `ssign run` interface. `extract_proteins.py` supports it
internally; v1.x will expose a paired-input flag.

## Skipping annotation tools

Any annotation tool can be skipped on the command line:

```bash
ssign run input.gbff --outdir results \
    --skip-signalp \
    --skip-hhsuite \
    --skip-plmblast
```

The pipeline still runs everything else and records "skipped" for the
matching columns in the output.

## Filtering or including secretion-system types

By default, Flagellum, Tad, and T3SS are excluded from substrate
identification. Override with a different list (space-separated):

```bash
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad T3SS
```

To include T3SS (e.g. when running on an organism known to carry it):

```bash
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad
```

T3SS is excluded by default because DeepSecE produces a high false-positive
rate for it (mostly flagellar misclassifications). The full rationale,
including the 74-genome benchmark behind the decision, is in
[`explanation/design_decisions.md`](../explanation/design_decisions.md) § 2.1
and § 3.3.

## Switching SignalP or DeepLocPro between local and webserver

ssign is offline-first: the canonical (default) path is a local install of
SignalP and DeepLocPro. If the binaries are on `PATH`, no flags are needed:

```bash
ssign run input.gbff --outdir results
```

If they live outside `PATH`, point ssign at them:

```bash
ssign run input.gbff --outdir results \
    --signalp-path /path/to/signalp6/bin \
    --deeplocpro-path /path/to/deeplocpro
```

If you do not have a DTU licence, opt into the DTU webserver fallback
(no licence needed on your part, internet required):

```bash
ssign run input.gbff --outdir results \
    --signalp-mode remote \
    --deeplocpro-mode remote
```

See [`how-to/install.md`](install.md) for the SignalP / DeepLocPro local
install steps.
