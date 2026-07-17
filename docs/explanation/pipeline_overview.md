# How the pipeline works

A narrative walkthrough of what ssign does during a run, in the order it
happens. For per-decision rationale and citations, see
[`design_decisions.md`](design_decisions.md). For the column-by-column
output reference, see [`reference/output_files.md`](../reference/output_files.md).

ssign is command-line only: `ssign run <genome> --outdir <out>`. It runs in
six phases, each producing an intermediate output the next phase consumes.
(Internally these expand into a ~16-step DAG, roughly one step per tool.)

```
input → proteins + gene order
      → secretion systems detected
      → secretion candidates predicted
      → substrates filtered by proximity + voting
      → optional functional annotation
      → integrated CSV + report + figures
```

## Phase 1: input processing

ssign accepts annotated GenBank, paired GFF3, or raw FASTA contigs.

For GenBank input, it reads the protein translations, locus tags, and each
protein's contig coordinates straight from the file. By default it also
re-annotates with Bakta, because incoming GenBank annotations have unknown
provenance (old Prokka, recent PGAP, manual curation, a private pipeline); one
consistent caller across a cohort is what makes downstream consensus voting
comparable. Curated-GenBank users can opt out with `--use-input-annotations`.

For FASTA input, Bakta (or pyrodigal as a fallback) calls ORFs from scratch.

This phase produces three core artefacts: a protein FASTA, a gene-info table
(one row per CDS with locus, contig, start, end, product), and a gene-order
table (the same proteins sorted by chromosome position).

## Phase 2: secretion-system detection

ssign hands the protein FASTA and gene order to **MacSyFinder v2** with the
**TXSScan** profile bundle. MacSyFinder scans for HMM hits to secretion-system
component genes (T1SS ABC, T2SS Gsp, T3SS Yop, T4SS VirB, T5SS
autotransporters, T6SS VAS, plus flagellar and Tad-pilus machinery), then
assembles them into complete systems per TXSScan's component requirements.

Two outputs:

- `valid_systems`: one row per system scoring above `--wholeness-threshold`
  (default 0.8 of required components present).
- `ss_components`: one row per individual component HMM hit, pre-grouped by
  which system instance it belongs to.

A genome can carry several systems of the same type (e.g. two T6SSs),
distinguished by their component locations. Flagella, Tad pili, and type-IV
pili are excluded at the next phase, not here (T3SS is not excluded);
detection itself is unconditional.

## Phase 3: secreted-protein prediction

Two independent predictors look at every protein and decide whether it looks
secreted:

- **DeepLocPro** predicts the subcellular localization (extracellular,
  outer-membrane, periplasmic, cytoplasmic, ...).
- **DeepSecE** predicts which secretion-system *type* a protein is a substrate
  of (T1SE / T2SE / T3SE / T4SE / T6SE).

A protein is flagged as a candidate substrate if **either** of these trips.
ssign records `n_prediction_tools_agreeing` (0-2) as a confidence signal, and
`secretion_evidence` lists which tools voted. The "any one" rule is deliberate:
a missed substrate is a missed biological finding, while false positives get
filtered by Phase 4's proximity step. See
[`design_decisions.md` § 3.1](design_decisions.md#31-equal-predictor-rule-dlp--dse-both-trigger).

**SignalP** also runs here but is evidence-only, not a trigger. It detects
classical Sec/Tat signal peptides, which many Gram-negative effectors (T3SS,
T4SS, T6SS, T1SS C-terminal signals) lack by design, so treating it as a
trigger would under-call those classes.

By default DLP, SignalP, and DeepSecE run only on proteins inside the SS
neighbourhood (Phase 4) to save compute. Three flags (`--dlp-whole-genome`,
`--sp-whole-genome`, `--dse-whole-genome`) run them across the whole genome
instead, which cohort-wide enrichment analysis needs.

ssign is offline-first: the canonical path uses local DLP and SignalP installs.
Users without a DTU academic licence can opt into the DTU webserver fallback
(`--deeplocpro-mode remote`, `--signalp-mode remote`; no licence needed,
internet required), but its long-term availability depends on DTU, so local
installs are the durable choice for publication and cohort work.

## Phase 4: substrate identification

This phase combines Phase 2's detected systems with Phase 3's secreted-looking
proteins and decides which are real candidate substrates. Two filters run in
series:

1. **Per-component proximity.** For each individual SS component detected in
   Phase 2, ssign takes a window of `--proximity-window` genes (default ±3)
   around that component on the *same contig*. The union across all components
   of one system instance is that system's "neighbourhood". Candidate
   substrates from Phase 3 inside the neighbourhood pass; candidates outside
   any neighbourhood are dropped. Using a per-component window rather than a
   span across the system's full footprint is load-bearing: an early bug used
   the system-wide span and produced ~26 false positives, because
   secretion-system genes can span tens of kb. See
   [`design_decisions.md` § 5.1](design_decisions.md#51-per-ss-component-window-not-full-system-span).

2. **Localization quorum.** For each candidate system, the fraction of detected
   components correctly localized (membrane-spanning OM components on the OM,
   etc.) must exceed `--required-fraction-correct` (default 0.8). Systems whose
   components look mislocalized are dropped from substrate calls.

T5SS handling is special: the autotransporter (T5aSS), two-partner secretion
(T5bSS), and chaperone-usher (T5cSS) classes export their own passenger domain
rather than a separate substrate. ssign's T5SS handler calls each component as
its own substrate when the per-component biology fits (T5aSS_PF03797 and
T5cSS_PF03895 pass on extracellular OR outer-membrane; the T5bSS translocator
pore is OM-only).

Output: `substrates_filtered` (per-genome) is the load-bearing CSV from this
phase.

## Phase 5: optional functional annotation

Six tools run independently against the substrate set. Each is opt-in via a
`--skip-*` flag, so a base-tier user without large databases can finish at
Phase 4 with no annotations.

| Tool | What it adds |
|---|---|
| **BLASTp** | Best hit description against Swiss-Prot (default) or NR (opt-in). |
| **HH-suite** | Best Pfam domain and best PDB structural homolog via HMM-vs-HMM. |
| **EggNOG-mapper** | Orthology, COG/KEGG categories, GO terms. |
| **InterProScan** | Protein domain calls across Pfam, SMART, PROSITE, etc. |
| **pLM-BLAST** | Embedding-based remote-homology search against ECOD30. |
| **ProtParam** | Physicochemical properties (GRAVY, MW, pI, charge, aromaticity). |

Each writes a per-tool CSV to the work directory. None depend on each other, so
they run in parallel when cores allow.

## Phase 6: integration and reporting

ssign reads every Phase 5 output and merges them into a single "integrated"
CSV keyed on `locus_tag` + `sample_id`. From that CSV it produces:

- **`<sample-id>_results.csv`**: chunked output (secreted proteins with merged
  annotation columns, then secretion systems with associated substrates, then
  other detected systems).
- **`<sample-id>_results_raw.csv`**: every column from the integrated CSV,
  unfiltered.
- **Annotation consensus voting.** Each substrate is classified into one of 27
  broad functional categories (Apparatus-associated, Protease/Peptidase,
  Pore-forming toxin, Adhesin, Autotransporter passenger, ...) by
  keyword-matching tool descriptions. The most-supported category becomes
  `broad_annotation`, with tool names listed as evidence. See
  [`design_decisions.md` § 4.1](design_decisions.md#41-27-category-broad-functional-voting).
- **Enrichment testing** (opt-in, `--enrichment-stats`). A per-SS-type
  circular-shift permutation test asks whether secreted-predicted proteins
  cluster around each SS type's components more than a genome-structure-
  preserving null (the exact circular rotations of the gene-ordered positivity
  vector), emitting fold + permutation p + BH q per system type.
  Output: `<sample-id>_enrichment_stats.tsv` (see
  [`output_files.md`](../reference/output_files.md)).
- **Up to six summary figures**, the curated `01`-`06` set: `01` secreted
  proteins by SS type, `02` size & physicochemical properties, and `03`-`06`
  functional categories (COG, KEGG, EggNOG, curated consensus). Toggle with the
  `--fig-*` flags.
- **HTML and text reports** that bring it together in human-readable form.

## Ortholog grouping

Alongside integration, ssign groups a genome's substrate proteins into ortholog
families. The step is substrate-scoped: it runs an all-vs-all BLASTp on that
genome's filtered substrate set (not on whole genomes), keeps hits above the
identity and coverage cutoffs, and clusters them single-linkage with Union-Find.
Each substrate gains `ortholog_group`, `og_n_members`, and `og_mean_pident`
columns in the integrated CSV, and the groups are written to
`<sample-id>_ortholog_groups.csv`.

It runs automatically, but only when NCBI BLAST+ (`blastp` and `makeblastdb`)
is on PATH. Without BLAST+ the step soft-skips: every substrate becomes its own
singleton group and the run continues.

A multi-genome run adds a second, cross-genome pass on top of the per-genome
grouping: it pools every genome's substrates, runs one all-vs-all BLASTp over the
combined set, and clusters that. It writes `cross_genome_orthologs.csv` (each
substrate with its `sample_id` and cross-genome group) and
`cross_genome_ortholog_groups.csv` (per group: `n_members`, `n_genomes`,
`genomes`, `members`, `mean_pident`), merges an `xg_ortholog_group` column into
every genome's integrated CSV, and draws `figures/07_cross_genome_orthologs.png`:
a two-panel conservation figure (group-size distribution + the most-conserved
groups, stacked by secretion-system type and labelled by function). Same BLAST+
requirement as the per-genome pass.
