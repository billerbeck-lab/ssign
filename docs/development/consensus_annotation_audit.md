# Consensus annotation audit (2026-07-09)

Audit of the custom consensus-annotation logic + keyword word-matching that drives the functional figures. Findings + a proposed v2. **No code changed** — this is for review before we build v2. Closes the open "REVISIT the consensus keyword grouping" flag (NOTES, 2026-06-30).

Graded against: the Xanthobacter 74-genome panel (593 substrates), the benchmark accuracy sheet (55 known secreted proteins with paired tool output + known function), the benchmark ground-truth families (17 effectors), and a focused literature pass on known effector functional classes.

Code under audit: `annotation_consensus.py` (`CATEGORY_PATTERNS`, `classify_description`, `compute_consensus`) and `ssign_lib/functional_vocab.py` (`consensus_bucket`).

---

## TL;DR grade: **C-** (usable signal, systematically biased, not paper-ready)

The pipeline usually produces *a* label, but the label is often wrong or uninformative for the exact proteins the paper cares about (effectors). Three structural problems dominate:

1. **Over-broad rules swamp function.** The Autotransporter, Structural, and Transporter regexes fire on 41% / 34% / 39% of all tool descriptions and together **win 60% of the 593 calls** (Autotransporter 230, Structural 123). Specific function is drowned out.
2. **Effectors get mislabelled as machinery.** The `t[1-9]ss` token in the "Secretion system" rule captures any protein whose annotation names its cognate system, so T6SS amidase effectors (Tae/Tae4) and multiple T3SS effectors (BopA, BipB, YspE, EspZ, CopN) are called *apparatus* and dropped from the functional chart.
3. **The vote is unsound.** `classify_description` returns *every* matching category, so one description splits its vote across unrelated buckets (809 multi-match instances over 593 rows); "Hypothetical" and title-case fallback noise are winnable categories that outvote real calls; correlated tools are counted as independent, inflating confidence.

Accuracy against known truth: **13/55** known secreted proteins land outside the vocabulary, and among those inside there are systematic confusions (RTX toxins → Protease; TPS adhesins → Nuclease; S-layer → Protease).

---

## How it works (for reference)

- `classify_description(desc)` runs 17 case-insensitive regexes and returns **all** that match; if none match, it mints a category from the **first 3 words, title-cased**.
- `compute_consensus({tool: desc})` classifies each tool's description, **frequency-counts categories across tools**, takes `most_common(1)` as `broad_annotation`, and sets confidence by how many tools hit that category (≥3 High, 2 Medium, 1 Low).
- `consensus_bucket()` (figures) passes known categories through, routes `Secretion system`/`Flagellar` → "Apparatus-associated", and dumps everything else (blanks + fallback noise) → "Other".

Five tools feed it: Bakta product, EggNOG, InterProScan, pLM-BLAST/ECOD, BLASTp.

---

## Findings

### A. Voting logic

| # | Problem | Evidence |
|---|---|---|
| A1 | **Vote-splitting.** `classify_description` returns all matches, so a multi-domain description casts several fractional votes across unrelated buckets. | 809 multi-category tool-matches over 593 rows. |
| A2 | **"Hypothetical" is winnable.** Two tools calling a protein hypothetical outvote the one tool with a real function. | 5/593 won "Hypothetical" with a real category present in another tool (`compute_consensus` has no floor semantics). |
| A3 | **Fallback noise is winnable.** The first-3-words fallback mints pseudo-categories ("Twin-Arginine Translocation Pathway", "Six-Bladed Beta-Propeller", "Tungstate Binding") that can win the vote. | 28/593 (5%) `broad_annotation` calls are fallback noise; a **COG2931 RTX-toxin** row was labelled "Six-Bladed Beta-Propeller" (InterProScan noise) instead of Toxin. |
| A4 | **Correlated tools counted as independent.** Bakta/EggNOG/InterProScan often echo the same description; "≥3 agree = High" overstates confidence. | 103/593 rated "High" despite redundant sources. |
| A5 | **Arbitrary tie-break.** Ties resolve by `Counter` insertion order (pattern order), not by specificity or evidence. | structural. |

### B. Regex rules (`CATEGORY_PATTERNS`)

| # | Problem | Evidence |
|---|---|---|
| B1 | **`fli[a-z]` false matches.** Catches "**fli**ppase" and "con**fli**ct" as Flagellar. | probe: `lipid II flippase MurJ → ['Flagellar']`; `conflict → ['Flagellar']`. |
| B2 | **Over-broad Autotransporter.** `barr?el.*domain` + `passenger.*domain` grab generic β-barrel/OM proteins. | fired on 681/1662 descriptions; won 230/593 calls. |
| B3 | **Over-broad Structural.** `outer.*membrane` treats OM *localisation* as a *function*. | fired 573×; won 123/593. |
| B4 | **Category overlap.** `substrate.binding` lives in both Transporter and Binding protein; `acyltransferase` (Lipase) vs `acetyltransferase`/`transferase` (Transferase). | probe: `substrate-binding protein → ['Transporter','Binding protein']`. |
| B5 | **Machinery hijack.** `secretion.*system\|type.*secretion\|t[1-9]ss\|vir[bd]` captures effectors named after their SS and routes them to "Apparatus-associated". | Tae/Tae4 amidase effectors, BopA/BipB/YspE/EspZ/CopN → "Secretion system". |
| B6 | **Dead rules.** Transferase (2 hits) and Flagellar (3, all false) barely fire; the granularity is wrong for the corpus. | rule-usage counts. |
| B7 | **"first match wins" comment is inaccurate** — `classify_description` returns all matches; order only affects the tie-break. | code read. |

### C. In-the-wild distribution (Xantho 593)

- Autotransporter 230 · Structural 123 · Nuclease 37 · Regulatory 35 · Glycoside hydrolase 33 · Protease 28 · blank 27 · Toxin 23 · … · fallback-noise 28.
- Figure buckets: Other 55, Apparatus-associated 7, and the two over-broad buckets (Autotransporter 230 + Structural 123) = **60% of everything**.
- Confidence: High 103 / Medium 319 / Low 144 / no_hits 27.

### D. Accuracy vs known truth (55 gold proteins)

Systematic, not random, errors:

- **RTX toxins → Protease**: `lktA` (Mannheimia leukotoxin), `rtxA` → Protease (should be Toxin).
- **TPS adhesins → Nuclease**: `hxuA`, `hmw1A`, `hmw2A`, `cdrA` (two-partner-secretion adhesins/hemophores) → Nuclease.
- **S-layer → Protease**: `rsaA` (Caulobacter S-layer) → Protease.
- **Phospholipase → Transferase**: `plaA` → Transferase.
- **13/55 land outside the vocabulary** (Other/fallback/blank), including several T3SS/T6SS effectors.

Correct where the annotation is unambiguous: serralysins → Protease; trimeric-AT adhesins (`yadA`,`nadA`,`sadA`) → Adhesin; classic autotransporters (`espP`,`flu`,`icsA`) → Autotransporter; `TseL` → Lipase.

### E. Coverage gaps vs the literature

The current vocabulary has no home (or only a too-generic one) for the effector activities that define each secretion system:

- **T6SS**: peptidoglycan **amidase/endopeptidase** (Tae) — not in the Protease regex → "Other" or hijacked to machinery; **muramidase/glucosaminidase** (Tge); **NAD(P)+ glycohydrolase**; Tle phospholipases land in Lipase (OK). ([Hernandez 2020](https://onlinelibrary.wiley.com/doi/10.1111/cmi.13241); [TseP eLife 2024](https://elifesciences.org/reviewed-preprints/101125))
- **T3SS**: **GAP/GEF**, **phosphothreonine lyase**, **ADP-ribosyltransferase**, **glycosyltransferase** (arginine-GlcNAc, e.g. NleB), **acetyltransferase**, **kinase**, **ubiquitin E3 ligase / deubiquitinase** — all collapse to generic "Transferase" or "Other". ([FEMS Microbiol Rev 2011](https://academic.oup.com/femsre/article/35/6/1100/524381); [Trends Microbiol 2021](https://www.cell.com/trends/microbiology/fulltext/S0966-842X(21)00262-6))
- **T1SS/T5SS**: **RTX toxin**, **hemophore/heme-binding**, **S-layer**, SPATE proteases — no dedicated category. ([ASM Bacterial Secretion Systems](https://journals.asm.org/doi/10.1128/microbiolspec.vmbf-0012-2015); [SecretEPDB](https://www.nature.com/articles/srep41031))

---

## Proposed v2 (for review — not implemented)

Three changes: fix the vote, fix the rules, and adopt a literature-grounded taxonomy.

### 1. Voting logic

- **One vote per tool.** `classify_description` returns the **single most-specific** category for a description (specific enzymatic/effector classes rank above generic localisation/structure), not every match. Kills vote-splitting (A1) and the Structural/Transporter swamp (B2/B3).
- **Hypothetical/Other are floors, never competitors.** A protein is "Hypothetical" only if *no* tool yields a functional category; it can never outvote one (A2).
- **Drop the title-case fallback entirely.** No match → `Unclassified`, never a mintable winning label (A3).
- **Deduplicate correlated tools.** Collapse identical/near-identical descriptions to one vote before counting; report `n_distinct_tool_calls` so "High" means genuine cross-source agreement (A4).
- **Specificity tie-break**, not insertion order (A5).

### 2. Machinery detection (separate from function)

- Remove `t[1-9]ss` / "type … secretion" from the *functional* matcher (B5). Detect apparatus components by **component identity** (VgrG/Hcp/PAAR/Tss*, Gsp*, Vir*B/D, Sec/Tat translocon, flagellar Fli/Flg/Flh gene families with word boundaries) — a curated machinery list, not a substring that eats effectors.

### 3. Literature-grounded category set

Replace the 17 ad-hoc buckets with a two-level scheme (fine class → figure bucket):

- **Hydrolases**: Protease/peptidase · Peptidoglycan hydrolase (amidase/muramidase, Tae/Tge) · Lipase/phospholipase (Tle) · Glycoside hydrolase / carbohydrate-active · Nuclease · Phosphatase.
- **Transferases (split out the effector activities)**: Glycosyltransferase · ADP-ribosyltransferase · Acetyl/acyltransferase · Kinase · Methyltransferase.
- **Host-signalling effectors**: GTPase modulator (GAP/GEF) · Phosphothreonine lyase · Ubiquitin-pathway (E3 ligase / deubiquitinase) · Effector (mimicry/other).
- **Toxins/membrane-active**: RTX toxin · Pore-forming/hemolysin · Other toxin.
- **Surface/binding**: Adhesin (fimbrial/trimeric-AT/TPS) · Hemophore / heme+metal acquisition · S-layer · Lectin/CBM.
- **Redox**: Oxidoreductase.
- **Machinery** (→ Apparatus-associated): secretion/uptake apparatus · Flagellar.
- **Low-information**: Transporter/channel · Regulatory · Chaperone · Unclassified · Hypothetical (floor).

### 4. Concrete regex fixes (carry into v2 patterns)

- `fli[a-z]` → `\bfl[ai][A-Z]\b`-style flagellar-gene match with word boundaries (B1).
- Autotransporter requires `autotransporter|passenger domain`, not bare `barrel domain` (B2).
- Drop bare `outer.*membrane` as a *functional* category; OM localisation is metadata, not function (B3).
- De-overlap `substrate.binding` (→ transport/solute-binding only) and the transferase family (B4).
- Add: `amidase`, `muramidase`, `glucosaminidase`, `ADP-ribosyl`, `glycosyltransferase`, `phosphothreonine lyase`, `deubiquitinase|ubiquitin ligase`, `RTX|repeat.in.toxin|calcium-binding repeat`, `hemophore|heme|TonB-dependent|siderophore`, `S-layer`, `TPS|two-partner|hemagglutinin` (E).

### Suggested validation for v2

Re-run this audit's grading script against v2 and require: (a) no single category > ~25% of calls, (b) 0 effectors routed to Apparatus by keyword, (c) ≥ the current correct count on the 55-protein accuracy sheet with fewer outside-vocab, (d) fallback-noise = 0.

---

## Proposed v2 — prototyped + measured (2026-07-09)

The v2 above was prototyped (`docs/development/v2lib.py`) and graded head-to-head with v1 on the Xantho 593, the 55-gold accuracy sheet, and the full-tier smoke-test-3 24-substrate set (`docs/development/regrade.py`). One finding reshaped it.

### The deepest finding: two of the tools report *structure*, not *function*

Feeding structural/fold annotations into a functional keyword matcher is the dominant error source. Evidence:
- **pLM-BLAST/ECOD** returns fold codes ("β-helix", "α+β") — pure structure.
- **HHpred_PDB** returns PDB titles that are misleading: the same GDSL lipase whose Pfam says `Lipase_GDSL_lke` gets PDB `Thermolabile hemolysin` (a structural homolog). Matching "hemolysin" → toxin is wrong.
- **InterProScan** mixes real domains (`Peptidase_S6`) with fold-superfamily names (`pectate lyase C-like domain superfamily`, `Haem peroxidase domain superfamily`) that trigger enzyme keywords for proteins with no such activity (autotransporter β-helix → "Glycoside hydrolase"; VasX toxin → "Oxidoreductase").
- **EggNOG, Bakta product, BLASTp/Swiss-Prot, GBFF** are the genuinely functional annotators. (Note: BLASTp/Swiss-Prot returned *nothing* on the environmental Roseixanthobacter substrates — too distant from Swiss-Prot's model organisms; relevant caveat for the paper.)

### v2 mechanism (as prototyped)

1. **One weighted vote per tool**, deduped by identical string. Tool credibility tiers: **3** = BLASTp/Swiss-Prot, EggNOG, Bakta, GBFF · **2** = HHpred_Pfam, InterProScan · **1** = pLM-BLAST/ECOD, HHpred_PDB.
2. **Fold-name down-weight**: any description that IS a structural fold name ("… superfamily", "…-like domain", "β-helix/barrel/propeller", "TIM barrel") is capped to weight 1, so a fold can't outvote a real functional call regardless of source tool.
3. **Specificity-ordered single category per description** (first match wins on a most-specific → generic list); no title-case fallback (→ `Unclassified`); `Hypothetical`/`Unclassified` are floors that never beat a functional call.
4. **Machinery by component identity** (VgrG/Hcp/PAAR/Tss/Gsp/Vir/Sec/Tat/flagellar with boundaries) + translocator context ("haemolysin *secretion* protein", ShlB/FhaC/HecB/TpsB) → `Apparatus-associated`, which **competes by weight** (so a translocator the good tools agree on wins over a weak structural mislabel).
5. **Taxonomy** = the two-level lit-grounded set above, minus a standalone "RTX toxin" bucket (RTX is a T1SS secretion motif spanning proteases/cytolysins/cyclases, not a function — serralysins → Protease, CyaA/cytolysins → Pore-forming), plus **Beta-lactamase** and **Phage/mobile element**.

### Measured v1 → v2

| | v1 | v2 |
|---|---|---|
| Title-case fallback noise (Xantho 593) | 37 | **0** |
| Effectors keyword-hijacked to machinery | several (Tae amidase, T3SS effectors) | **0** (Tae → Peptidoglycan hydrolase; translocators → Apparatus by evidence) |
| Structure-vs-function false matches | e.g. 71 autotransporter→GH, hemolysin-fold→toxin | fixed (GDSL lipase stays Lipase; flu→Adhesin; icsA→Autotransporter) |
| 55-gold outside-vocabulary | 13/55 | **8/55**, all now *honest* Unclassifieds (no tool has a function) |
| Correct on annotatable gold | mixed (RTX→Protease errors, adhesins→Nuclease) | serralysins→Protease, SPATE→Protease, Map→GTPase modulator, hasA→Hemophore, cdrA/flu→Adhesin, cya→Pore-forming |

### Residual errors are garbage-in (the one lever left)

A few stay wrong because the annotation tools themselves have no correct call: `Tde1` → Phage (its His-Me nuclease domain is genuinely annotated phage-associated), `VasX` → Oxidoreductase (only a "haem peroxidase domain" fold hit exists), `apxIIIA` → Protease (pore-former carrying a peptidase domain). No keyword logic fixes these; they need a **small curated effector-family override** (known gene/family → category), the recommended follow-on to v2.

Prototype: `docs/development/v2lib.py` (matcher) + `docs/development/regrade.py` (grading). Implementation tracked as the v2 task.

---

## Sources

- [Hernandez et al. 2020, T6SS effector proteins (Cell Microbiol)](https://onlinelibrary.wiley.com/doi/10.1111/cmi.13241)
- [TseP chimeric amidase/lysozyme effectors (eLife 2024)](https://elifesciences.org/reviewed-preprints/101125)
- [Activity, delivery, diversity of T6SS effectors (Mol Microbiol)](https://onlinelibrary.wiley.com/doi/10.1111/mmi.14648)
- [Galán/Dean, T3SS effector domains & motifs (FEMS Microbiol Rev 2011)](https://academic.oup.com/femsre/article/35/6/1100/524381)
- [T3SS effector network hypothesis (Trends Microbiol 2021)](https://www.cell.com/trends/microbiology/fulltext/S0966-842X(21)00262-6)
- [Bacterial Secretion Systems overview (ASM Microbiol Spectr)](https://journals.asm.org/doi/10.1128/microbiolspec.vmbf-0012-2015)
- [SecretEPDB effector knowledgebase (Sci Rep)](https://www.nature.com/articles/srep41031)

Grading script: `docs/development/grade_consensus.py` (Xantho + benchmark; reads the local `~/Desktop/cx3_runs` panel + `validation_sweeps/benchmark` sheets — adjust paths to re-run).
