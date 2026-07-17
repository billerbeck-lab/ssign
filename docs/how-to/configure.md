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

`ssign` will be available from PyPI at the v1.0.0 release; until then,
install from source (git clone + `pip install -e .`).

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

GFF3 works via `ssign run` when a same-stem nucleotide FASTA sits
alongside it (e.g. `x.gff` + `x.fna`).

### `--use-input-annotations`: leave it off for most use cases

The flag tells ssign to trust the input GenBank's `/product` qualifiers
and skip Bakta re-annotation. Don't pass it unless you have a specific
reason to: most public GenBank records (including NCBI's K-12 reference)
lack the Pfam-domain annotations MacSyFinder needs to detect secretion
systems, so the substrate recall drops sharply.

On the K-12 reference genome:

| Configuration            | Substrates detected |
|--------------------------|---------------------|
| Default (Bakta re-annotation) | 17 |
| `--use-input-annotations`     | 8  |

Same drop applies to GFF3-derived inputs (NCBI's GFF3 also doesn't carry
Pfam IDs in the format ssign needs). Only use `--use-input-annotations`
when you've manually curated a GenBank record with Pfam-tagged CDS
qualifiers and want them preserved.

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

By default, Flagellum, Tad, and the type-IV pili / uptake appendages
(T4aP, T4bP, MSH, ComM, Archaeal-T4P) are excluded from substrate
identification. T3SS is **not** excluded (it is detected and
substrate-called by default). Override with a different list
(space-separated):

```bash
# the default set, written out explicitly:
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P
```

To also exclude T3SS (restore the pre-2026 behaviour), append it:

```bash
ssign run input.gbff --outdir results \
    --excluded-systems Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P T3SS
```

T3SS is included by default, but DeepSecE is never trusted for it (it
produces a high false-positive rate, mostly flagellar misclassifications);
T3SS calls rely on TXSScan detection + DeepLocPro + proximity instead. The
full rationale, including the 74-genome benchmark behind the decision, is in
[`explanation/design_decisions.md`](../explanation/design_decisions.md) § 2.1
and § 3.3.

## Tuning the proximity window

`--proximity-window N` (default `3`) controls how many genes on either
side of each detected SS component are considered candidate substrates.
ssign's proximity-based approach detects substrates that are
**co-located** with their secretion system in the chromosome; widening
the window recovers more candidates near each component but does not
help with substrates encoded far from the secretion machinery.

```bash
ssign run input.gbff --outdir results --proximity-window 5
```

On the K-12 reference: window 3 → 17 substrates, window 5 → 20.

**Compute cost grows roughly linearly** with the window size. The
neighborhood feeds the prediction tools (DeepLocPro, DeepSecE,
SignalP), so a window-5 run does ~1.7x more prediction work than
the default; window-10 is ~3x. Per-substrate annotation tools
(InterProScan, EggNOG, pLM-BLAST) only run on the final filtered set,
so their cost grows more slowly. False-positive count also grows with
the window, substrates further from a system component are more often
co-located by chance than because of secretion biology.

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

If you do not have a DTU licence, you can opt into the DTU webserver
fallback (no licence needed on your part, internet required):

```bash
ssign run input.gbff --outdir results \
    --signalp-mode remote \
    --deeplocpro-mode remote
```

The webserver path is a convenience for users without a DTU licence and
for first-time trials. It depends on DTU continuing to host the service,
which we cannot guarantee long-term, for paper-grade or long-running
cohorts, install SignalP and DeepLocPro locally instead.

See [`how-to/install.md`](install.md) for the SignalP / DeepLocPro local
install steps.
