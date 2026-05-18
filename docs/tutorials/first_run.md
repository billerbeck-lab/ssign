# Example first run

A walkthrough that takes you from a fresh machine to a results table. We
use *Escherichia coli* K-12 substrain MG1655 because it is the canonical
Gram-negative reference: ~4.6 Mb, ~4,300 CDS, secretion systems and
substrates that are well-described in the literature, so you can sanity-
check what ssign reports against what you already know.

About 10-20 minutes on a typical laptop with the canonical local install
of DeepLocPro and SignalP. The "quick taste" variant below skips the
local DTU install and submits to DTU's webserver instead; that route adds
a network-bound wait and is fine for a single-genome trial run.

## Prerequisites

- Linux or macOS
- Python 3.10 or newer (`python3 --version`)
- Local installs of SignalP 6.0 and DeepLocPro, set up per
  [`how-to/install.md`](../how-to/install.md). Both are free with a DTU
  academic licence and run fully offline. If you do not have a DTU
  licence yet, the "Quick taste without DTU tools" section near the end
  shows a webserver variant you can use instead.
- ~10 GB of free disk: ~3 GB for the install (PyTorch + ESM) and ~7 GB
  downloaded on first DeepSecE run for the cached ESM-1b language model

> **Running on an HPC cluster?** Don't run from a login node or a
> JupyterHub session — both are typically throttled to ~1 CPU and DeepSecE
> will take 60-90 min instead of seconds. Submit a proper compute job (and
> request a GPU if you can; DeepSecE is ~100x faster on GPU). See
> [`how-to/run_on_hpc.md`](../how-to/run_on_hpc.md) for templates.

## 1. Install ssign

```bash
python3 -m venv ~/.ssign-env
source ~/.ssign-env/bin/activate
pip install ssign
ssign --version
```

The last line prints the installed version (e.g. `ssign 0.9.0`). If you
see `command not found`, your venv is not active; rerun the `source` line.

## 2. Download E. coli K-12

The reference assembly is GCF_000005845.2 (ASM584v2). Fetch the GenBank
file from NCBI's RefSeq FTP:

```bash
mkdir -p ~/ssign-tutorial && cd ~/ssign-tutorial

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
mv GCF_000005845.2_ASM584v2_genomic.gbff ecoli_k12.gbff

ls -lh ecoli_k12.gbff
```

The file is ~13 MB.

## 3. Run ssign

The base tier (`pip install ssign` with no extras) covers secretion-system
detection, secreted-protein prediction, proximity analysis, and reporting.
With SignalP and DeepLocPro installed locally per `install.md`, the
canonical run command just skips the annotation tools that need their
own databases (BLASTp, HH-suite, EggNOG, InterProScan, pLM-BLAST):

```bash
ssign run ecoli_k12.gbff --outdir ecoli_results \
    --use-input-annotations \
    --skip-blastp \
    --skip-eggnog \
    --skip-hhsuite \
    --skip-interproscan \
    --skip-plmblast
```

What that command says, in plain English:

- `ecoli_k12.gbff` is the input.
- Write outputs to `./ecoli_results/`.
- Trust the input GenBank's existing annotations rather than re-annotating
  with Bakta. (Bakta is not installed in the base tier.)
- SignalP and DeepLocPro default to local mode, so ssign will look for
  the binaries on `PATH` (or whichever directory you gave to
  `--signalp-path` / `--deeplocpro-path`).
- Skip the five annotation tools that need extra databases. Everything
  else (DeepLocPro, SignalP, DeepSecE predictions, MacSyFinder detection,
  proximity analysis, ProtParam, reporting) runs as normal.

ssign will print step-by-step progress with percentages. DeepSecE is
usually the slowest step on a CPU; DeepLocPro and SignalP take a few
minutes each on a modern laptop when run locally.

### Quick taste without DTU tools

If you have not installed SignalP and DeepLocPro locally yet and just
want to see ssign produce output, swap in the DTU webserver. No licence
required on your part, but the run will take longer (network round-trip
plus DTU queueing) and depends on DTU continuing to host the service.
The webserver is fine for a one-off trial; for paper-grade or
multi-genome work, install locally per `install.md` so your pipeline is
not at the mercy of a third-party service.

```bash
ssign run ecoli_k12.gbff --outdir ecoli_results \
    --use-input-annotations \
    --signalp-mode remote \
    --deeplocpro-mode remote \
    --skip-blastp \
    --skip-eggnog \
    --skip-hhsuite \
    --skip-interproscan \
    --skip-plmblast
```

In this variant the DLP and SignalP steps submit batches to DTU and
poll for results; expect ~30 minutes total wall time, mostly waiting.

## 4. Look at the output

Once the run finishes, `ecoli_results/` contains:

```
ecoli_results/
├── ecoli_k12_results.csv       Main results
├── ecoli_k12_results_raw.csv   All columns, no filtering
├── ecoli_k12_summary.txt       Plain-text report
├── figures/
│   └── ecoli_k12/*.png         Five summary figures
└── .ssign/
    └── ecoli_k12_progress.json Resume manifest
```

Open `ecoli_k12_results.csv` in a spreadsheet, or use the flat
`ecoli_k12_results_raw.csv` from pandas (the chunked CSV interleaves
sections with different column sets, so a plain `read_csv` doesn't
work):

```python
import pandas as pd
df = pd.read_csv("ecoli_results/ecoli_k12_results_raw.csv")
print(df[["locus_tag", "predicted_localization", "nearby_ss_types", "gbff_annotation"]].head(20))
```

The CSV is organised in up to three chunks, each prefixed with a `#`
header. Empty chunks are omitted:

1. `# Secreted Proteins` — one row per predicted substrate, with
   annotation columns from every tool that ran.
2. `# Secretion Systems (with secreted proteins)` — the SS instances
   whose neighbourhoods contained at least one secreted protein.
3. `# Secretion Systems (other)` — SS instances detected but without
   high-confidence substrates.

For column-by-column meaning, see [`reference/output_files.md`](../reference/output_files.md).

### Sanity check: what should we see in K-12?

K-12 has well-described secretion machinery you can use to gut-check the
output:

- **Type 2 secretion (T2SS):** the *gsp* operon. ssign detects a complete
  T2SS in the secretion-systems chunk, with the *gspC–gspM* + *gspO*
  components and at least one substrate (chitinase ChiA / b3338).
- **Type 5a secretion (T5aSS):** the canonical *E. coli* autotransporters
  — adhesin Ag43 (*flu* / b2000), YfaL, YpjA, YcgV, YhjY. ssign reports
  these as both detected systems and their own substrates ("T5SS-self").
- **Type 6 secretion (T6SS)** is *not* detected at default settings — the
  K-12 MG1655 T6SS operon is H-NS-suppressed and degenerate, so it falls
  below MacSyFinder's wholeness threshold. Lowering
  `--wholeness-threshold 0.5` will surface it as an incomplete system.
- **Flagella** are excluded by default (they are not a substrate-export
  system in the relevant sense). You will not see them in the chunked
  CSV. They are still listed in the raw CSV.

`ecoli_k12_summary.txt` is the same information laid out for reading: a
substrate count, per-SS breakdown, and the enrichment-test results.

## 5. What to try next

- **Add an annotation tool.** Re-run with `--skip-interproscan` removed
  (after [installing InterProScan](../how-to/install.md#interproscan))
  and you get InterPro domain calls in the substrate rows.
- **Run a different organism.** The same command works on any Gram-negative
  GenBank or FASTA. *Pseudomonas aeruginosa* PAO1 (GCF_000006765.1) and
  *Vibrio cholerae* N16961 (GCF_000006745.1) are common follow-ups.
- **Run a cohort.** ssign is built to run over many genomes; the CLI is
  designed for shell loops or HPC arrays. See
  [`how-to/run_on_hpc.md`](../how-to/run_on_hpc.md) for a SLURM/PBS
  template.
- **Tune what counts as a substrate.** See
  [`how-to/configure.md`](../how-to/configure.md) for the common knobs
  (proximity window, system completeness threshold, which SS types to
  include or exclude).
