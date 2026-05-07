# How the pipeline works

A narrative walkthrough of what ssign does during a run, in the order it
happens. For per-decision rationale and citations, see
[`design_decisions.md`](design_decisions.md). For the column-by-column
output reference, see [`reference/output_files.md`](../reference/output_files.md).

ssign runs in six phases. Each phase produces an intermediate output that
the next phase consumes. The whole flow:

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

For GenBank input, it extracts the protein translations, the locus tags,
and each protein's contig coordinates straight from the file. By default
it also re-annotates the input with Bakta because incoming GenBank
annotations have unknown provenance: a file from "GenBank" might have
been annotated by old Prokka, recent PGAP, manual curation, or a private
pipeline. Re-annotating with one consistent caller across an entire
cohort is what makes downstream consensus voting comparable. Users with
carefully curated GenBanks can opt out with `--use-input-annotations`.

For FASTA input, Bakta (or pyrodigal as a fallback) calls ORFs from
scratch.

This phase produces three core artefacts: a protein FASTA, a gene-info
table (one row per CDS with locus, contig, start, end, product), and a
gene-order table (the same proteins sorted by chromosome position).

## Phase 2: secretion-system detection

ssign hands the protein FASTA and the gene order to **MacSyFinder v2**
with the **TXSScan** profile bundle. MacSyFinder scans for HMM hits to
secretion-system component genes (T1SS ABC, T2SS Gsp, T3SS Yop, T4SS
VirB, T5SS autotransporters, T6SS VAS, plus flagellar and Tad-pilus
machinery), then tries to assemble those hits into complete systems
based on TXSScan's per-system component requirements.

Two outputs:

- `valid_systems`: one row per system that scored above the
  `--wholeness-threshold` (default 0.8 of required components present).
- `ss_components`: one row per individual component HMM hit,
  pre-grouped by which system instance it belongs to.

A genome can carry several systems of the same type (e.g. two T6SSs),
and they get distinguished by their component locations. Flagella, Tad
pili, and T3SS are excluded by default at the next phase, not here;
detection itself is unconditional.

## Phase 3: secreted-protein prediction

Three independent predictors look at every protein and decide whether
it looks secreted:

- **DeepLocPro** predicts the subcellular localization (extracellular,
  outer-membrane, periplasmic, cytoplasmic, ...).
- **DeepSecE** predicts which secretion-system *type* a protein is a
  substrate of (T1SE / T2SE / T3SE / T4SE / T6SE).
- **PLM-Effector** uses a stack of protein language models (ESM-1b,
  ESM-2, ProtT5) plus XGBoost to classify per-SS-type secretion.

A protein is flagged as a candidate substrate if **any one** of these
three trips. ssign records `n_prediction_tools_agreeing` (0-3) as a
confidence signal, and `secretion_evidence` lists which tools voted.
The "any one" rule is deliberate: false negatives from any single tool
are expensive (a missed substrate is a missed biological finding),
while false positives get filtered by Phase 4's proximity step. See
[`design_decisions.md` § 3.1](design_decisions.md#31-equal-predictor-rule-dlp--dse--plm-e-all-trigger).

**SignalP** also runs in this phase but does not contribute to the
trigger count. It detects classical Sec/Tat signal peptides, which
many Gram-negative effectors (T3SS, T4SS, T6SS, T1SS C-terminal
signals) lack by design. Treating it as a trigger would under-call
those effector classes; recording it as evidence preserves the
information value without biasing the call.

By default DLP, SignalP, and DeepSecE only run on proteins inside the
SS neighbourhood (Phase 4) to save compute. Three flags
(`--dlp-whole-genome`, `--sp-whole-genome`, `--dse-whole-genome`) run
them across the entire genome instead, which is needed for cohort-wide
enrichment analysis.

ssign is offline-first: the canonical execution path uses local DLP and
SignalP installs. Users without a DTU academic licence can opt into the
DTU webserver fallback with `--deeplocpro-mode remote` and
`--signalp-mode remote` (no licence needed on the user's part, internet
required).

## Phase 4: substrate identification

This phase combines Phase 2's "we found a secretion system" with Phase 3's
"this protein looks secreted" and decides which proteins are real
candidate substrates.

Two filters run in series:

1. **Per-component proximity.** For each individual SS component
   detected in Phase 2, ssign takes a window of `--proximity-window`
   genes (default ±3) around that component on the *same contig*. The
   union across all components of one system instance is the
   "neighbourhood" of that system. Candidate substrates from Phase 3
   that fall inside the neighbourhood pass; candidates outside any
   neighbourhood are dropped. This per-component window (rather than a
   span across the system's full footprint) is load-bearing: an early
   bug used the system-wide span and produced ~26 false positives
   because secretion-system genes can span tens of kb. See
   [`design_decisions.md` § 5.1](design_decisions.md#51-per-ss-component-window-not-full-system-span).

2. **Localization quorum.** For each candidate system, the fraction of
   detected components correctly localized (membrane-spanning OM
   components on the OM, etc.) must exceed
   `--required-fraction-correct` (default 0.8). Systems whose
   components look mislocalized are dropped from substrate calls.

T5SS handling is special: the autotransporter (T5aSS), two-partner
secretion (T5bSS), and chaperone-usher (T5cSS) classes export their own
passenger domain rather than a separate substrate. ssign's T5SS handler
calls each component as its own substrate when the per-component
biology fits (T5aSS_PF03797 and T5cSS_PF03895 pass on extracellular OR
outer-membrane; T5bSS translocator pore is OM-only).

Output: `substrates_filtered` (per-genome) is the load-bearing CSV from
this phase.

## Phase 5: optional functional annotation

Six tools run independently against the substrate set. Each is opt-in
via a `--skip-*` flag, so a base-tier user without large databases can
finish the pipeline at Phase 4 with no annotations.

| Tool | What it adds |
|---|---|
| **BLASTp** | Best hit description against NR or Swiss-Prot. |
| **HH-suite** | Best Pfam domain match and best PDB structural homolog via HMM-vs-HMM. |
| **EggNOG-mapper** | Orthology, COG/KEGG categories, GO terms. |
| **InterProScan** | Protein domain calls across all Pfam, SMART, PROSITE, etc. |
| **pLM-BLAST** | Embedding-based remote-homology search against ECOD70. |
| **ProtParam** | Physicochemical properties (GRAVY, MW, pI, charge, aromaticity). |

Each writes a per-tool CSV to the work directory. None depend on each
other; with enough cores they could run in parallel (and do, in the
GUI's batch mode).

## Phase 6: integration and reporting

ssign reads every Phase 5 output and merges them into a single
"integrated" CSV keyed on `locus_tag` + `sample_id`. From that CSV it
produces:

- **`<sample-id>_results.csv`**: chunked output (secreted proteins
  with merged annotation columns, then secretion systems with
  associated substrates, then other detected systems).
- **`<sample-id>_results_raw.csv`**: every column from the integrated
  CSV with no filtering.
- **Annotation consensus voting.** Each substrate is classified into one
  of 17 broad functional categories (Adhesin, Autotransporter,
  Protease, ...) by keyword-matching tool descriptions. The
  most-supported category becomes `broad_annotation`; tool names get
  listed as evidence. See
  [`design_decisions.md` § 4.1](design_decisions.md#41-17-category-broad-functional-voting).
- **Enrichment testing.** A Fisher's exact test asks which functional
  categories are over-represented near each SS type, vs the
  whole-genome background. Output: `enrichment_fisher.csv` and a plain-
  text summary.
- **Five summary figures** (functional category, SS-component
  composition, tool-coverage heatmap, substrate-count per SS,
  functional summary). Toggle each with the `--fig-*` flags.
- **HTML and text reports** that bring it together in a human-readable
  form.

## Cross-genome step (multi-genome runs only)

When the GUI processes more than one genome in a single batch, ssign
adds a cross-genome ortholog-grouping step after each genome's Phase 6
finishes. It pools substrates across genomes, runs all-vs-all BLASTp,
and groups orthologs with Union-Find. The group ID is then merged back
into each genome's integrated CSV (`ortholog_group` column) so a user
can see "this T1SS substrate in genome A is the same protein family as
that one in genome B".

This step is opt-in by the GUI's batch flow; CLI users running one
`ssign run` per genome do not get it automatically (the union is across
the run's `--outdir`, so the genomes have to share an outdir to be
grouped).
