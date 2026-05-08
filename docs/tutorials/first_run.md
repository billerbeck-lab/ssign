# Your first ssign run: E. coli K-12

A walkthrough that takes you from a fresh machine to a results table. We
use *Escherichia coli* K-12 substrain MG1655 because it is the canonical
Gram-negative reference: ~4.6 Mb, ~4,300 CDS, secretion systems and
substrates that are well-described in the literature, so you can sanity-
check what ssign reports against what you already know.

About 30 minutes, mostly waiting for DTU's webserver.

## Prerequisites

- Linux or macOS
- Python 3.10 or newer (`python3 --version`)
- An internet connection (this tutorial uses the DTU webserver fallback for
  DeepLocPro and SignalP, so no DTU licence is required on your part)

No GPU, no DTU licence, no system packages required. For a fully offline
run with local DTU tools (the canonical ssign mode), see
[`how-to/install.md`](../how-to/install.md).

## 1. Install ssign

```bash
python3 -m venv ~/.ssign-env
source ~/.ssign-env/bin/activate
pip install ssign
ssign --version
```

The last line should print something like `ssign 1.0.0`. If you see
`command not found`, your venv is not active; rerun the `source` line.

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
For a tutorial run we opt into the DTU webserver fallback for DLP and
SignalP (no licence needed) and skip the annotation tools that need their
own databases (BLASTp, HH-suite, EggNOG, InterProScan, pLM-BLAST):

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

What that command says, in plain English:

- `ecoli_k12.gbff` is the input.
- Write outputs to `./ecoli_results/`.
- Trust the input GenBank's existing annotations rather than re-annotating
  with Bakta. (Bakta is not installed in the base tier.)
- Use the DTU webserver fallback for SignalP and DeepLocPro, so we don't
  need a DTU academic licence for this first run.
- Skip the five annotation tools that need extra databases. Everything
  else (DeepLocPro, SignalP, PLM-Effector predictions, MacSyFinder
  detection, proximity analysis, ProtParam, reporting) runs as normal.

ssign will print step-by-step progress with percentages. The two slowest
steps are DeepLocPro and SignalP, both submitted to DTU's webserver and
queued behind other users (a few minutes each, usually).

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

Open `ecoli_k12_results.csv` in a spreadsheet or pandas:

```python
import pandas as pd
df = pd.read_csv("ecoli_results/ecoli_k12_results.csv", skiprows=1, comment="#")
print(df[["locus_tag", "predicted_localization", "nearby_ss_types", "gbff_annotation"]].head(20))
```

The CSV is organised in three blocks separated by `# Section` headers:

1. **Secreted Proteins**: one row per predicted substrate, with annotation
   columns from every tool that ran.
2. **Secretion Systems (with secreted proteins)**: the SS instances whose
   neighbourhoods contained at least one secreted protein.
3. **Other Secretion Systems**: SS instances detected but without
   high-confidence substrates.

For column-by-column meaning, see [`reference/output_files.md`](../reference/output_files.md).

### Sanity check: what should we see in K-12?

K-12 has well-described secretion machinery you can use to gut-check the
output:

- **Type 2 secretion (T2SS):** the *gsp* operon. ssign should detect a
  T2SS in section 2 of the CSV.
- **Type 6 secretion (T6SS):** the H-NS-suppressed *evf* / *vas* loci.
  ssign should detect a T6SS in section 2 or 3.
- **Flagella** are excluded by default (they are not a substrate-export
  system in the relevant sense). You will not see them in the chunked
  CSV. They are still listed in the raw CSV.

`ecoli_k12_summary.txt` is the same information laid out for reading: a
substrate count, per-SS breakdown, and the enrichment-test results.

## 5. What to try next

- **Add an annotation tool.** Re-run with `--skip-interproscan` removed
  (after [installing InterProScan](../how-to/install.md#interproscan-optional))
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
