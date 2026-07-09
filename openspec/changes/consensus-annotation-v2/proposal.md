## Why

The functional-annotation consensus that drives the paper's functional figures was audited (`docs/development/consensus_annotation_audit.md`) and graded **C-**. It systematically mislabels the effectors the paper cares about: 60% of calls are won by three over-broad rules, effectors named after their secretion system are hijacked into the "machinery" bucket, structural fold names (from pLM-BLAST/ECOD and InterProScan superfamilies) are matched as if they were functions, and a title-case fallback mints noise categories that outvote real calls. A prototyped v2 fixes these on the audit's own benchmark data; this change lands it in production.

## What Changes

- **Tool-credibility weighting**: one deduped, weighted vote per annotation tool. Tier-1 functional-name tools (BLASTp/Swiss-Prot, EggNOG, Bakta, GBFF) outweigh Tier-2 domain tools (HHpred_Pfam, InterProScan) which outweigh Tier-3 structure-only tools (pLM-BLAST/ECOD, HHpred_PDB).
- **Fold-name down-weight**: a description that IS a structural fold/superfamily name is capped to weight 1, so structure can never outvote function.
- **One specificity-ranked category per description** (most-specific → generic, first match wins). **BREAKING** to output values: removes the first-3-words title-case fallback (unmatched → `Unclassified`); `Hypothetical`/`Unclassified` become floors that never beat a functional call.
- **Machinery by component identity + translocator context** (VgrG/Hcp/PAAR/Tss/Gsp/Vir/Sec/Tat/flagellar with word boundaries; ShlB/FhaC/HecB/TpsB) → `Apparatus-associated`, competing by weight. Removes the `t[1-9]ss` keyword that hijacked effectors.
- **Literature-grounded taxonomy**: drop the standalone `RTX toxin` bucket (RTX is a secretion motif, not a function); add Peptidoglycan hydrolase, ADP-ribosyltransferase, Glycosyltransferase, Phosphothreonine lyase, Ubiquitin-pathway, GTPase modulator, Beta-lactamase, Hemophore/metal-uptake, S-layer, Phage/mobile-element.
- `consensus_bucket` keeps routing machinery → Apparatus and unknown → Other, tracking the new category set.
- **Out of scope** (separate follow-on): a curated effector-family override for garbage-in residuals where the tools themselves carry no correct functional call (Tde1, VasX, apxIIIA).

## Capabilities

### New Capabilities
- `annotation-consensus`: the tool-weighted functional-consensus classifier that assigns each substrate a broad functional category from its multi-tool annotation descriptions, plus the figure-side bucket collapse.

### Modified Capabilities
<!-- none: run-figures consumes the consensus but its requirement (show functional categories) is unchanged; only the vocabulary the consensus emits changes. -->

## Impact

- Code: `src/ssign_app/scripts/annotation_consensus.py` (rewrite `CATEGORY_PATTERNS`, `classify_description`, `compute_consensus`), `src/ssign_app/scripts/ssign_lib/functional_vocab.py` (`consensus_bucket`, category set). `integrate_annotations.py` already passes the 8 tool→column names the weighting keys on (no change needed there).
- Output: the `broad_annotation` / `broad_consensus_annotation` / figure-bucket values change vocabulary (new categories; no more title-case noise; `Unclassified` instead of invented labels). Downstream figures (`run-figures`) pick up the new categories automatically via `consensus_bucket`.
- Tests: `tests/unit/test_annotation_consensus.py`, `tests/unit/test_functional_vocab.py`.
- Validation: `docs/development/regrade.py` against the Xantho 593 + 55-gold + full-tier-24 fixtures.
- No new runtime dependencies; no pipeline-step or DB changes.
