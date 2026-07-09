## 1. Port the v2 matcher into annotation_consensus.py

- [x] 1.1 Replace `CATEGORY_PATTERNS` with the v2 ordered (most-specific → generic) taxonomy from `docs/development/v2lib.py`, including the new effector categories and the dropped RTX bucket.
- [x] 1.2 Rewrite `classify_description` to return a single category (first match wins) or None — no title-case fallback; add the `__APPARATUS__` translocator-context rule and the component-identity machinery detector (`_MACH`).
- [x] 1.3 Add the `_IS_FOLD` structural-fold detector.

## 2. Weighted voting in compute_consensus

- [x] 2.1 Add the `TOOL_WEIGHT` tier map (BLASTp/EggNOG/Bakta/GBFF=3, HHpred_Pfam/InterProScan=2, pLM-BLAST/HHpred_PDB=1; default 2).
- [x] 2.2 Rewrite `compute_consensus`: dedup identical descriptions keeping max weight; cap fold-name descriptions to weight 1; weighted vote where functional categories AND Apparatus compete, `Hypothetical` is a floor, no functional/apparatus → `Unclassified`.
- [x] 2.3 Preserve the output dict keys (`broad_annotation`, `broad_consensus_annotation`, `detailed_annotation`, `n_tools_agreeing`, `concordance_ratio`, `confidence_tier`, etc.) so `integrate_annotations` and the CSV schema are unchanged; confidence recomputed from the winner's weighted support (High ≥5, Medium ≥3, Low else; `no_hits` sentinel kept).
- [x] 2.4 Confirm the tool→description keys `compute_consensus` receives (from `integrate_annotations.TOOL_HIT_COLUMNS`) match the `TOOL_WEIGHT` names — all 8 match, no adjustment needed.

## 3. Figure-side vocabulary

- [x] 3.1 Update `ssign_lib/functional_vocab.py` `consensus_bucket` + `_KNOWN_CATEGORIES` to the v2 category set (imports `CATEGORY_NAMES`); `Apparatus-associated` passes through, `Unclassified`/blank/unknown → `Other`.
- [x] 3.2 Grep for any hard-coded old category names in figures/report/docs — none found (the one `generate_figures.py` match is the SS-type axis label, unrelated).

## 4. Tests

- [x] 4.1 Rewrite `tests/unit/test_annotation_consensus.py`: weighted vote, fold down-weight, single-category classify, floors, machinery-by-identity (no `t[1-9]ss` hijack), RTX-by-activity, output schema.
- [x] 4.2 Update `tests/unit/test_functional_vocab.py` for the new `consensus_bucket` category set.
- [x] 4.3 `pytest tests/unit/ -q` green — **1449 passed**.

## 5. Golden fixture + validation

- [x] 5.1 Regenerate the end-to-end golden fixture (`tests/fixtures/golden/`): only the GBFF-driven consensus columns of the one substrate row change (`Autotransporter` → `Autotransporter passenger`); regenerated surgically (verified equal to `compute_consensus`) since the golden run has no local DeepLocPro here. Documented in `REGENERATE.md`.
- [x] 5.2 Run `docs/development/regrade.py` — **GATE PASSED**: port fidelity 672/672 rows identical to the prototype; Xantho top category 39% (<42%); v2 fallback-noise 0 (v1 had 37); full-tier-24 fold-name false calls fixed (Thermolabile-hemolysin→Lipase, ShlB→Apparatus); 55-gold outside-vocab 13→3.
- [x] 5.3 Open Question resolved (design.md): autotransporters classify by **passenger function** (espP→Protease, flu/cdrA→Adhesin), falling to `Autotransporter passenger` only when unspecific (icsA). Matches the regenerated golden.
