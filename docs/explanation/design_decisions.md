# Design decisions and defensible citations

Living notes to seed the paper's Methods section. Each entry records a
decision we've made in `ssign`, the scientific rationale, and the
peer-reviewed citations that back it. When paper-writing starts, pull
from here directly.

Entries are organised by pipeline stage. Each has three parts:

- **Decision** — what we do
- **Rationale / evidence** — why, with peer-reviewed citations
- **Methods-section draft** — drop-in paragraph (edit as needed)

---

## 1. CDS calling

### 1.1 Sole canonical CDS caller: Bakta (via Pyrodigal)

- **Decision:** `ssign` uses Bakta as its sole canonical CDS caller. Bakta
  internally wraps Pyrodigal (a maintained Cython rebinding of
  Prodigal v2.6.3+31b300a) and adds a custom small-ORF (<30 aa) detector.
  No consensus across multiple CDS predictors is performed.

- **Rationale:**
  - **Pyrodigal ≡ Prodigal, but maintained.** Pyrodigal is a
    regression-tested rebinding of Prodigal v2.6.3, not a
    reimplementation. Pyrodigal is actively released
    (v3.7.1, March 2026) while the original Prodigal binary has been
    near-dormant since 2016 — Pyrodigal is the responsible choice for a
    5-year reproducibility target.
  - **Prodigal's accuracy on Gram-negatives is well-characterised.**
    On _E. coli_ K-12 with experimentally verified genes: 100%
    sensitivity on the verified set, 96.5% correct start-site
    identification (Hyatt et al., 2010). In the 15-tool ORForise
    benchmark, Prodigal ranked 1st overall across 12 metrics on 3/6
    model organisms including three Gram-negatives (Dimonaco et al.,
    2022).
  - **Multi-caller consensus has been shown to hurt in this regime.**
    Dimonaco et al. (2022) explicitly tested the union of the five
    best callers across seven organisms and found it improved gene
    detection by only **0.47% on average while substantially increasing
    false-positive CDS**. Their recommendation: "we advise against
    generating aggregated ab initio annotations from multiple tools
    where no existing annotation is available for the genome, as this
    results in poor overall performance."
  - **Inter-tool disagreement is almost entirely start-codon
    placement, not gene presence.** Korandla et al. (2020) (AssessORF)
    showed 88–95% of CDS calls agree between top-tier predictors with
    no clear winner; disagreements are dominated by start-codon choice
    within otherwise-concordant ORFs. Small N-terminal shifts do not
    affect secretion-substrate identification because downstream HMM,
    signal-peptide, and proximity analyses are robust to a few-residue
    start-site difference.
  - **No license-clean alternative adds independent signal.** Prokka
    (Seemann, 2014) wraps the same Prodigal; its author has publicly
    endorsed Bakta as the successor. PGAP (NCBI) bundles GeneMarkS-2+
    under a Georgia Tech non-redistribution clause, making it
    incompatible with a public SHA-pinned Docker image. GeneMarkS-2
    standalone has the same license problem. Pyrodigal-gv is
    giant-virus-tuned, not bacterial.
  - **ML-based CDS callers are not production-ready for a 2031
    reproducibility target.** BALROG (Sommer & Salzberg, 2021) gives
    ~11% false-positive reduction at 100× runtime, but its GitHub
    repository shows near-total maintenance abandonment since early
    2023 (zero releases, unanswered issues, first author moved on).
    Recent protein-language-model pipelines (Bacformer,
    Protein Set Transformer) _consume_ Pyrodigal output rather than
    replacing it — a strong signal that the field still regards
    Prodigal-class callers as canonical.

- **Citations:**
  - [Hyatt et al. (2010). _BMC Bioinformatics_ 11:119](https://doi.org/10.1186/1471-2105-11-119) — Prodigal.
  - [Larralde (2022). _JOSS_ 7(72):4296](https://doi.org/10.21105/joss.04296) — Pyrodigal.
  - [Schwengers et al. (2021). _Microbial Genomics_ 7(11):000685](https://doi.org/10.1099/mgen.0.000685) — Bakta.
  - [Dimonaco et al. (2022). _Bioinformatics_ 38(5):1198–1207](https://doi.org/10.1093/bioinformatics/btab827) — ORForise benchmark.
  - [Korandla et al. (2020). _Bioinformatics_ 36(4):1022–1029](https://doi.org/10.1093/bioinformatics/btz714) — AssessORF.
  - [Yok & Rosen (2011). _BMC Bioinformatics_ 12:20](https://doi.org/10.1186/1471-2105-12-20) — consensus CDS calling on metagenomic reads.
  - [Sommer & Salzberg (2021). _PLOS Comput Biol_ 17(2):e1008727](https://doi.org/10.1371/journal.pcbi.1008727) — BALROG.

- **Methods-section draft:**

  > Coding sequences were predicted by Bakta v1.x (Schwengers et al.,
  > 2021), which internally uses Pyrodigal (Larralde, 2022) — a
  > maintained Cython rebinding of Prodigal v2.6.3 (Hyatt et al., 2010)
  > — plus a custom small-ORF (<30 aa) detector. We did not perform
  > multi-caller consensus voting. The most comprehensive recent
  > benchmark of fifteen prokaryotic gene predictors (Dimonaco et al., 2022) found that union-aggregation of the five top callers
  > improved gene detection by only 0.47% on average while
  > substantially increasing false-positive CDS, and concluded against
  > aggregated ab initio annotation for genomes without pre-existing
  > reference annotations. Inter-tool disagreement in this regime is
  > dominated by start-codon selection within otherwise-concordant ORFs
  > (Korandla et al., 2020), which does not materially affect
  > secretion-substrate identification because downstream HMM,
  > signal-peptide, and proximity analyses are robust to small
  > N-terminal coordinate shifts.

### 1.1.1 EggNOG-mapper is not a CDS caller

- **Decision:** EggNOG-mapper is used strictly as a functional-annotation
  tool on pre-called proteins via `--itype proteins`. It is never used
  as a CDS caller, even though it technically supports two gene-
  prediction modes (`--genepred prodigal` and `--genepred search`).

- **Rationale:** The eggNOG-mapper v2 paper (Cantalapiedra et al., 2021)
  states explicitly: "ORF detection, only available for prokaryotic
  assemblies, is performed using Prodigal." Gene prediction is listed
  as the 4th feature in the paper's abstract — a convenience add-on
  for users who want to start from contigs, not a primary capability.
  The paper makes **no accuracy claims for gene prediction** and does
  **not benchmark it against other callers**. The official GitHub
  README describes eggNOG-mapper as "a tool for fast functional
  annotation of novel sequences." `--genepred prodigal` would be
  exactly equivalent to Bakta's Pyrodigal (same upstream tool).
  `--genepred search` (blastx-based inference against the eggNOG
  database) has no peer-reviewed accuracy data published and is
  homology-biased by construction — genes without eggNOG hits are
  missed entirely and 5' boundaries are unreliable because the mode
  has no explicit start-codon model.

- **Defensive coding:** `run_eggnog.py` passes `--itype proteins`
  explicitly, so EggNOG cannot silently fall into gene-prediction mode
  even if its default ever changes.

- **Citations:**
  - [Cantalapiedra, C. P., et al. (2021). _Molecular Biology and Evolution_ 38(12):5825–5829](https://doi.org/10.1093/molbev/msab293) — eggNOG-mapper v2, Prodigal integration noted as convenience feature.
  - eggnog-mapper GitHub README (accessed 2026-04-24) — describes the tool as a functional annotator.

### 1.2 Re-annotate by default (Phase 3.3)

- **Decision:** When the input is a GenBank file, `ssign` re-runs Bakta
  on the DNA sequences by default and treats Bakta's fresh CDS set as
  canonical. The user's original GenBank annotations are preserved as
  a separate `gbff_annotation` column (mapped to Bakta's CDS by
  reciprocal coordinate overlap) and vote independently in
  annotation-consensus. Users with curated GenBanks can skip the Bakta
  re-run via `--use-input-annotations`.

- **Rationale:** Incoming GenBanks have untrusted provenance — they
  may come from old Prokka runs, private pipelines, manual curation,
  or RefSeq at any vintage. Uniform Bakta re-annotation across a
  cohort is required for reproducible consensus voting; without it,
  "Bakta said X, EggNOG said Y, GBFF said Z" is unsound because the
  three opinions might concern slightly different CDS.

- **Citations:** as above (Bakta, Pyrodigal).

---

## 2. Secretion-system detection

### 2.1 MacSyFinder v2 + TXSScan

- **Decision:** Secretion systems are detected using MacSyFinder v2
  with the TXSScan HMM profile models. Only systems with a
  `wholeness_score ≥ 0.8` are accepted as valid. Flagellum and Tad (and
  the type-IV pili / uptake appendages) are excluded from substrate
  identification by default. T3SS is detected and substrate-called by
  default; DeepSecE is gated for it instead (see §3.3).

- **Rationale:**
  - MacSyFinder v2 is the current community standard for
    secretion-system detection, a full rewrite of the original tool
    focused on better edge-case handling and bundled profile packaging.
  - TXSScan provides curated HMM profiles for T1SS–T6SS plus flagellar
    and Tad-pilus systems.
  - **`wholeness_score ≥ 0.8` threshold**: Systems below 0.8 lack too
    many essential components to be confidently called. The threshold
    matches TXSScan's recommended quorum for Gram-negative genomes.
  - **Excluded systems (default):**
    - _Flagellum:_ not a secretion system in the substrate-export
      sense; its components are not secretion machinery for proteins
      to be functionally secreted. Including it pollutes substrate
      lists with flagellar proteins.
    - _Tad pilus:_ similar reasoning — its cargo is structural, not
      secreted effectors.
    (T3SS is **not** in the default excluded set — it is detected and
    substrate-called by default. The DeepSecE reliability issue below
    (§3.3) is handled by gating DeepSecE for T3SS, not by excluding the
    system.)

- **Citations:**
  - [Neron et al. (2023). _Peer Community Journal_ 3:e28](https://doi.org/10.24072/pcjournal.250) — MacSyFinder v2.
  - [Abby et al. (2016). _Scientific Reports_ 6:23080](https://doi.org/10.1038/srep23080) — TXSScan.

---

## 3. Secreted-protein prediction

### 3.1 Equal-predictor rule: DLP / DSE / PLM-E all trigger

- **Decision (Phase 3.2.b):** DeepLocPro, DeepSecE, and PLM-Effector
  are treated as **equal secretion predictors**. A protein is marked
  as a candidate substrate if **any one** of these three flags
  secretion. The count of agreeing tools is recorded as
  `n_prediction_tools_agreeing` (0-3) as a confidence signal, and
  `secretion_evidence` lists which tools flagged.

- **Rationale:**
  - **DLP** (localisation-based, trained on ~40k bacterial proteins)
    and **DSE** (secretion-type, trained on secretion-system
    effectors) measure different biological signals. A protein can be
    clearly secreted by cellular localisation without being a known
    effector, and vice versa — treating them as equal predictors
    captures both views.
  - **PLM-E** adds a third independent signal based on
    protein-language-model ensemble classification
    (ESM-1b + ESM-2 + ProtT5 + XGBoost stacking per SS type).
  - **Equality with OR-logic** is the current best practice for
    secretion prediction, where false negatives from any single tool
    are expensive (missed substrates) and false positives get filtered
    by downstream proximity analysis.

- **Citations:**
  - [Moreno et al. (2024). _Bioinformatics_ 40(12):btae677](https://doi.org/10.1093/bioinformatics/btae677) — DeepLocPro.
  - [Zhang et al. (2023). _Research_ 6:0258](https://doi.org/10.34133/research.0258) — DeepSecE.
  - [Zheng (2026). _Briefings in Bioinformatics_ 27(2):bbag143](https://doi.org/10.1093/bib/bbag143) — PLM-Effector.

- **Revision (2026-06-17): PLM-Effector demoted to opt-in.** Validation on
  P. aeruginosa PAO1 showed PLM-Effector calls ~25% of the proteome as effectors
  (18% even gated at `max_stacking >= 0.8`), with no reliable enrichment near real
  secretion systems. This is the documented failure mode of effector predictors
  trained/thresholded on balanced sets and applied genome-wide (the model's paper
  reports no genome-scale validation; PAO1 truly has ~4 T3SS effectors). Therefore:
  DLP + DSE are the default equal predictors; **PLM-Effector is OFF by default**
  (opt in with `--no-skip-plm-effector`), and when used its call is **gated at
  `max_stacking >= 0.8`** for consistency with DLP/DSE. It is excluded from the
  enrichment statistical test and retained as an opt-in proximity signal and a
  classifier-training feature. See openspec change
  `enrichment-background-and-plme-default-off`.

- **Revision (2026-06-19): enrichment test is a circular-shift permutation, not a
  binomial.** The opt-in enrichment test asks, per SS type, whether secreted-predicted
  proteins cluster around that type's components more than chance. The old binomial
  assumed each neighborhood protein is independently secreted at the genome background
  rate; secreted genes actually cluster (operons, islands), so that null is too narrow
  and over-states significance (fleet validation: real null SD 14.6 vs 10.7 under
  independence; the binomial gives the same fold but an anti-conservative p). The test
  now rotates each predictor's gene-ordered positivity vector through every circular
  offset (the exact permutation null, computed in one FFT pass; offset 0 is the
  observed) and reports fold + permutation p + BH q. Because it needs every gene's
  positivity in gene order, `--enrichment-stats` forces whole-genome DLP/DSE
  (~13 min/genome). Autotransporters (T5aSS/T5cSS) are tested by self-detection (the
  component is both machinery and substrate), other types by a ±3-gene window. A
  per-type null-distribution figure is emitted. See openspec change
  `enrichment-circular-shift-per-run`.

### 3.2 SignalP is evidence-only, not a trigger

- **Decision (Phase 3.2.b):** SignalP is recorded as a separate
  `signalp_supports_secretion` column but does **not** contribute to
  `is_secreted` or `n_prediction_tools_agreeing`.

- **Rationale:** SignalP detects **Sec/Tat signal peptides**, which
  cover only a subset of secretion pathways. Many Gram-negative
  effectors (T3SS, T4SS, T6SS, T1SS C-terminal signals) lack classical
  signal peptides, so SignalP correctly reports "no signal peptide"
  for them. Treating SignalP as a trigger would under-call these
  effector classes. Keeping it as evidence preserves its information
  value without biasing the secretion call.

- **Citations:**
  - [Teufel et al. (2022). _Nature Biotechnology_ 40(7):1023–1025](https://doi.org/10.1038/s41587-021-01156-3) — SignalP 6.0.

### 3.3 T3SS: DeepSecE not trusted (flagellar-misclassification guard)

- **Decision:** T3SS predictions from DeepSecE are flagged
  (`dse_T3SS_flagged=True`) and excluded from the secretion-evidence
  count **unconditionally** — independent of genome content, because a
  genome-content condition would re-admit the flagellar false positives.
  T3SS itself is **detected and substrate-called by default** (it is
  _not_ in `excluded_systems`); its substrate calls rely on
  MacSyFinder/TXSScan detection + DeepLocPro + proximity, not DeepSecE.
  (Earlier versions excluded T3SS entirely by default; that was changed
  once TXSScan detection proved reliable — openspec `t3ss-detection`.)

- **Why flagellar proteins specifically:** DeepSecE has no flagellum
  output class — its classes are non-secreted, T1SE, T2SE, T3SE, T4SE,
  T6SE (`run_deepsece.py` `PREDICTED_LABELS`). The flagellar export
  apparatus is an evolutionary T3SS homolog, so flagellar proteins share
  sequence features with injectisome effectors and, having no flagellum
  bin to land in, get funneled into T3SE. That is why the false
  positives concentrate in T3SS rather than the other effector types,
  and why the guard is T3SS-specific.

- **Rationale (internal benchmark):** In a 74-genome _Xanthomonas_
  benchmark, MacSyFinder found **0 T3SS systems** while DeepSecE
  predicted **1,808 T3SS substrate candidates** across the same
  cohort — mostly hypothetical proteins and flagellar-protein
  misclassifications. This false-positive rate is too high to accept
  uncritically; gating DSE T3SS calls on MacSyFinder validation
  prevents large-scale over-calling. The reliability issue is T3SS-
  specific; DSE calls for T1SS/T2SS/T4SS/T6SS are not gated.

- **Note for paper:** The 74-genome benchmark is internal. Either
  include these numbers as supplementary data (lab-generated,
  citable as Reid et al., in preparation), or replicate on a
  published _Xanthomonas_ set before submission.

### 3.4 Why PLM-Effector is vendored

- **Decision:** PLM-Effector source code is vendored into
  `src/ssign_app/scripts/plm_effector/` rather than installed as a
  dependency.

- **Rationale:**
  - No PyPI package exists upstream; install is conda-based with
    Python 3.9 + CUDA 11.3 pinned, which conflicts with `ssign`'s
    Python 3.10+ baseline.
  - Upstream code has hardcoded paths to the author's institutional
    machine (`/home/zhengdd/...`) that require edits to run anywhere
    else.
  - Vendoring lets us adapt for Python 3.10+, wrap a clean CLI around
    the entry point, and ship a SHA-pinned Docker image with 5-year
    reproducibility.
  - License: CC-BY 3.0 (not MIT, despite the upstream README badge
    claiming MIT; the actual `LICENSE` file is CC-BY 3.0). CC-BY 3.0
    explicitly permits redistribution with attribution, which we
    preserve via `src/ssign_app/scripts/plm_effector/LICENSE` and a
    citation entry in `CITATION.cff`.

---

## 4. Annotation consensus

### 4.1 17-category broad functional voting

- **Decision:** Each substrate protein is classified into one of 17
  broad functional categories (Adhesin, Autotransporter, Protease,
  Lipase/Esterase, Nuclease, Glycoside hydrolase, Toxin, Transporter,
  Secretion system, Flagellar, Oxidoreductase, Transferase, Chaperone,
  Binding protein, Structural, Regulatory, Hypothetical) by keyword-
  matching tool-specific description strings. Multiple tools vote;
  the most-supported category becomes `broad_annotation`, with tool
  names listed as evidence.

- **Rationale:**
  - 17 categories is a pragmatic compromise between "too broad to be
    useful" (4-5 COG-style superclasses) and "too granular to vote"
    (hundreds of KEGG/GO terms). Keyword matching against descriptions
    is tool-agnostic: any new tool that produces a human-readable
    description string can contribute votes without a schema change.
  - Confidence tiers (High ≥3 tools, Medium =2, Low =1, None =0)
    provide a simple interpretable confidence signal that maps
    directly to the `confidence_tier` column in the final output.

- **Tools voting (8 sources):** Bakta product, EggNOG description,
  BLASTp top hit description, HH-suite top Pfam description, HH-suite
  top PDB description, InterProScan descriptions, pLM-BLAST top ECOD
  hit, original GenBank annotation (GBFF). SignalP is excluded from
  annotation voting because it's a prediction tool, not an
  annotation tool.

### 4.2 Cross-tool field naming

- **Decision:** Bakta and EggNOG output tables use unified field names
  for cross-referenced codes: `ec_numbers`, `cog_ids`, `go_terms`,
  `kegg_ko` (not `kegg_ids`), `pfam_ids`.

- **Rationale:** Both EggNOG's `KEGG_ko` column and Bakta's
  `DbXrefs:KEGG:...` entries carry KEGG Orthology IDs. Using a single
  name (`kegg_ko`) across tools means `annotation_consensus.py` does
  not need per-tool column-name conditionals.

### 4.3 T5aSS passenger-domain annotation (default since v1.0.0)

- **Decision:** For T5aSS (classical autotransporter) substrates with
  a clean domain geometry (PF03797 β-barrel found, passenger length
  ≥ 25 aa), EggNOG / BLASTp / pLM-BLAST / HHsuite / ProtParam read
  the **passenger subsequence** rather than the full protein.
  InterProScan is the lone exception: it keeps reading the full
  protein. Users can opt into a second, full-protein pass via
  `--t5ass-annotate-whole`, which adds side-by-side `t5ass_whole_*`
  columns to the master CSV; the two annotation sets sit alongside
  each other rather than overwriting.

- **Rationale:** Autotransporters have three functional regions
  (N-terminal passenger / α-helical linker / C-terminal β-barrel
  translocator), and the barrel is highly conserved across the family.
  Tools that return a single best hit per query (DIAMOND-based EggNOG,
  blastp's bitscore-sorted top hit, pLM-BLAST's cosine-screened ECOD
  best hit) preferentially match the barrel — for a typical AT, this
  means the headline annotation collapses to "Outer membrane
  autotransporter barrel" (COG3210 or similar) and the functional
  passenger annotation (pertactin / IgA protease / hemoglobin protease
  fold etc.) is hidden. Empirical evidence on K-12 (2026-06-04 run):
  5/5 T5aSS substrates returned β-barrel ECOD families (4KH3, 5KE1,
  4MEE, 3QQ2, 3AEH) when run on the full protein.

  Routing those five tools to the passenger subsequence pushes the
  homology search toward passenger-fold matches. The passenger
  boundary is `seq[sp_end+1 : barrel_start - LINKER_LENGTH]`, where
  `sp_end` is SignalP's cleavage-site position, `barrel_start` is the
  envelope start of PF03797 from the bundled HMM scan, and
  `LINKER_LENGTH = 30` (constants in `scripts/ssign_lib/constants.py`).

  InterProScan stays on the full protein because it is already
  domain-aware: it reports per-domain coordinates and Pfam hits
  separately, so a full-AT input gives the user both the passenger
  Pfam(s) (pertactin etc.) and the barrel Pfam at distinct positions —
  passenger-only input would lose the barrel confirmation without
  gaining any functional information.

- **Fallback to full protein** for T5aSS substrates with:
  - any `t5_quality_flag` set (`barrel_only`, `no_signalp`,
    `no_sec_signal`, `omp_porin_no_at`, `unclassified`) — the
    handler couldn't fix the geometry and we don't trust the slice;
  - `passenger_length < MIN_PASSENGER_FOR_ANNOTATION` (25 aa) —
    pLM-BLAST's cpc-90 pre-screen and blastp's bitscore would not
    find meaningful homology on a residue count this small.

  Non-T5aSS substrates (T5bSS, T5cSS, T1SS effectors, etc.) always
  carry their full sequence — these don't have an autotransporter
  passenger to begin with.

- **Per-row provenance:** the `t5_annotation_source` column in the
  master CSV records `passenger` or `full` for each T5aSS row (empty
  for non-T5aSS). Populated by the side-car TSV from
  `_step_build_passenger_fasta`.

- **`--t5ass-annotate-whole` opt-in second pass:** users who want to
  compare passenger vs full-AT annotations side by side enable this
  flag. The pipeline runs EggNOG / BLASTp / pLM-BLAST / HHsuite /
  ProtParam a second time on the full protein FASTA, restricted to
  clean T5aSS substrates only, and writes the results into
  `t5ass_whole_*`-prefixed columns. v1.0.0 runs this second pass
  sequentially per genome (small T5aSS substrate counts make
  per-tool startup the dominant cost); pooled multi-genome execution
  for the second pass is a future enhancement.

- **Implementation notes:**
  - Per-genome step: `_step_build_passenger_fasta` reads the t5ss
    handler's domains TSV, substitutes passenger sequences into a
    new FASTA keyed by the same locus_tags, and writes a side-car
    TSV recording per-locus routing.
  - Routing dispatch: `PipelineRunner._annotation_input_proteins(tool)`
    returns the passenger-substituted FASTA for the five routed
    tools and the full FASTA for IPS (and any future non-routed
    tool).
  - Multi-genome: `MultiGenomeRunner._pool_segment_d_inputs` builds
    two pooled substrate-only FASTAs (`pooled_substrate_proteins.faa`
    full + `pooled_substrate_proteins_passenger.faa` substituted)
    and the pool runner's annotation steps route via the same helper.

---

## 5. Proximity analysis

### 5.1 Per-SS-component window, not full system span

- **Decision:** Neighborhood proteins are defined as ±3 genes from
  **each individual SS component's gene**, taking the union — not ±3
  from the full system-wide span.

- **Rationale:** Using the system-wide span caused ~26 false positive
  substrates in early benchmarking (internal). SS components can be
  scattered across tens of kb (e.g. T1SS genes split into an operon
  plus an outer-membrane channel located elsewhere), so a "full span"
  window pulls in unrelated neighbouring genes. Per-component ±3 is
  biologically grounded: most T1SS/T4SS substrates are encoded
  immediately adjacent to one of the component genes.

- **Cross-contig guard:** Proximity windows never cross contig
  boundaries. For fragmented assemblies this prevents proteins on
  different contigs from appearing as "neighbours" of an SS component.

### 5.2 T5aSS/T5cSS "hitchhikers"

- **What:** A classical or trimeric autotransporter (T5aSS/T5cSS) is a
  single self-contained polypeptide: its passenger threads through its
  own outer-membrane β-barrel, so the component *is* the substrate (the
  "self" signal). Separately, secreted-predicted proteins are found
  significantly enriched in the ±3-gene window around T5a/c components.
  We call these neighbours **hitchhikers**, on the hypothesis that they
  may piggyback out through the shared T5 pore rather than encoding their
  own machinery.

- **How it surfaces:** Because proximity runs its ±3 window around every
  SS component (§ 5.1), T5a/c neighbours are already substrate-called
  (`substrate_source == "proximity"`, `nearby_ss_types` = the T5 type).
  The enrichment test (`--enrichment-stats`) reports the hitchhiker
  signal as a separate `mode=window` result for T5aSS/T5cSS, alongside
  the `mode=self` autotransporter test (DLP-or-DSE for the hitchhiker
  window, DLP-or-SignalP for self). Figures show the two populations as
  distinct `(self)` / `(hitchhiker)` categories.

- **Status:** The hitchhiker enrichment is descriptive (it quantifies the
  clustering); the pore-sharing mechanism is a hypothesis, not a claim.

---

## 6. Pipeline architecture decisions

### 6.1 Offline-first

- **Decision:** `ssign` v1.0.0 has no external API dependencies. All
  tools run from local binaries and databases; no remote calls.

- **Rationale:** Publication longevity. External services (NCBI remote
  BLAST, EBI InterProScan web, MPI Toolkit HHpred, BioLib-hosted DTU
  tools) get deprecated, rate-limited, or altered over a 5-10 year
  window. A local-only pipeline is the only architecture that can
  reliably reproduce a 2026 analysis in 2031.

### 6.2 Install tiers (base / extended / full)

- **Decision:** Three install tiers differ in which databases
  `scripts/fetch_databases.sh` downloads:
  - **base** (~17 GB): MacSyFinder + DLP + DSE + SignalP + PLM-E +
    Bakta light
  - **extended** (~130 GB): + EggNOG + HH-suite (Pfam + PDB70) +
    InterProScan + pLM-BLAST
  - **full** (~630 GB): + BLAST NR + Bakta full DB + HH-suite UniRef30

- **Rationale:** Lab researchers rarely need BLAST NR (390 GB). Tier-
  aware distribution keeps the minimum-useful install under 20 GB
  while still offering the full reproducibility bundle for users who
  need it.

### 6.3 Nextflow "power mode" deprecated

- **Decision:** The Nextflow DSL2 pipeline (`bin/`, `modules/local/`,
  `workflows/`, `main.nf`, `nextflow.config`) is being removed for
  v1.0.0. The Python `ssign run` CLI + Singularity covers all HPC
  batch use cases.

- **Rationale:** Maintaining two orchestrators (Python runner +
  Nextflow) forced double-bookkeeping of every script between `bin/`
  and `src/ssign_app/scripts/`. Migration guidance for Nextflow users
  will land alongside the dedicated HPC how-to in the Phase 5 docs
  overhaul.

---

## 7. License and distribution

### 7.1 GPL-3.0-or-later

- **Decision:** `ssign` is distributed under GPL-3.0-or-later.

- **Rationale:** The binding constraint is Pyrodigal (imported
  directly in `extract_proteins.py:157`), which is GPL-3. All other
  direct dependencies (biopython, pandas, numpy, streamlit, pyhmmer)
  are permissive (BSD/MIT/Apache) and compatible with GPL-3.
  Switching to Apache 2.0 or MIT would require replacing
  `import pyrodigal` with a subprocess `prodigal` call — not worth the
  engineering for a public research tool where GPL-3 is standard.

---

## 8. Reference fixture genome

- **Decision:** Regression and integration tests use `contig_87` of
  _Xanthobacter tagetidis_ ATCC 700314 (213 kb, 179 CDS). BIMENO_04457
  on this contig is the expected T1SS substrate.

- **Rationale:** T1SS is Teo's research focus. The source assembly is
  fragmented (87 contigs); contig_87 was chosen because it contains
  the full T1SS operon plus the substrate on a single contiguous DNA
  span. The whole contig fits inside 213 kb, which is a reasonable
  git-tracked fixture size.

---

## Consolidated bibliography

For paper Methods section:

- Abby, S. S., et al. (2016). _Scientific Reports_, 6:23080. doi:[10.1038/srep23080](https://doi.org/10.1038/srep23080). — TXSScan.
- Cantalapiedra, C. P., et al. (2021). _Molecular Biology and Evolution_, 38(12):5825–5829. doi:[10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293). — eggNOG-mapper v2.
- Dimonaco, N. J., et al. (2022). _Bioinformatics_, 38(5):1198–1207. doi:[10.1093/bioinformatics/btab827](https://doi.org/10.1093/bioinformatics/btab827). — ORForise benchmark; no-aggregation recommendation.
- Hyatt, D., et al. (2010). _BMC Bioinformatics_, 11:119. doi:[10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119). — Prodigal.
- Korandla, D. R., et al. (2020). _Bioinformatics_, 36(4):1022–1029. doi:[10.1093/bioinformatics/btz714](https://doi.org/10.1093/bioinformatics/btz714). — AssessORF.
- Larralde, M. (2022). _Journal of Open Source Software_, 7(72):4296. doi:[10.21105/joss.04296](https://doi.org/10.21105/joss.04296). — Pyrodigal.
- Lomsadze, A., et al. (2018). _Genome Research_, 28(7):1079–1089. doi:[10.1101/gr.230615.117](https://doi.org/10.1101/gr.230615.117). — GeneMarkS-2.
- Moreno, J., et al. (2024). _Bioinformatics_, 40(12):btae677. doi:[10.1093/bioinformatics/btae677](https://doi.org/10.1093/bioinformatics/btae677). — DeepLocPro.
- Neron, B., et al. (2023). _Peer Community Journal_, 3:e28. doi:[10.24072/pcjournal.250](https://doi.org/10.24072/pcjournal.250). — MacSyFinder v2.
- Schwengers, O., et al. (2021). _Microbial Genomics_, 7(11):000685. doi:[10.1099/mgen.0.000685](https://doi.org/10.1099/mgen.0.000685). — Bakta.
- Seemann, T. (2014). _Bioinformatics_, 30(14):2068–2069. doi:[10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153). — Prokka.
- Sommer, M. J., & Salzberg, S. L. (2021). _PLOS Computational Biology_, 17(2):e1008727. doi:[10.1371/journal.pcbi.1008727](https://doi.org/10.1371/journal.pcbi.1008727). — BALROG.
- Teufel, F., et al. (2022). _Nature Biotechnology_, 40(7):1023–1025. doi:[10.1038/s41587-021-01156-3](https://doi.org/10.1038/s41587-021-01156-3). — SignalP 6.0.
- Yok, N. G., & Rosen, G. L. (2011). _BMC Bioinformatics_, 12:20. doi:[10.1186/1471-2105-12-20](https://doi.org/10.1186/1471-2105-12-20). — Combining gene prediction methods.
- Zhang, Y., et al. (2023). _Research_, 6:0258. doi:[10.34133/research.0258](https://doi.org/10.34133/research.0258). — DeepSecE.
- Zheng, D. (2026). _Briefings in Bioinformatics_, 27(2):bbag143. doi:[10.1093/bib/bbag143](https://doi.org/10.1093/bib/bbag143). — PLM-Effector.

---

_Last updated: with the Phase 3.2.d runner wiring + the Phase 3.3 CDS-calling research. Add new decisions here as they're made._
