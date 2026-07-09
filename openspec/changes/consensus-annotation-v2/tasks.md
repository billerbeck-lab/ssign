## 1. Port the v2 matcher into annotation_consensus.py

- [ ] 1.1 Replace `CATEGORY_PATTERNS` with the v2 ordered (most-specific → generic) taxonomy from `docs/development/v2lib.py`, including the new effector categories and the dropped RTX bucket.
- [ ] 1.2 Rewrite `classify_description` to return a single category (first match wins) or None — no title-case fallback; add the `__APPARATUS__` translocator-context rule and the component-identity machinery detector (`_MACH`).
- [ ] 1.3 Add the `_IS_FOLD` structural-fold detector.

## 2. Weighted voting in compute_consensus

- [ ] 2.1 Add the `TOOL_WEIGHT` tier map (BLASTp/EggNOG/Bakta/GBFF=3, HHpred_Pfam/InterProScan=2, pLM-BLAST/HHpred_PDB=1; default 2).
- [ ] 2.2 Rewrite `compute_consensus`: dedup identical descriptions keeping max weight; cap fold-name descriptions to weight 1; weighted vote where functional categories AND Apparatus compete, `Hypothetical` is a floor, no functional/apparatus → `Unclassified`.
- [ ] 2.3 Preserve the output dict keys (`broad_annotation`, `broad_consensus_annotation`, `detailed_annotation`, `n_tools_agreeing`, `concordance_ratio`, `confidence_tier`, etc.) so `integrate_annotations` and the CSV schema are unchanged; recompute confidence from the winning category's weighted support.
- [ ] 2.4 Confirm the tool→description keys `compute_consensus` receives (from `integrate_annotations.TOOL_HIT_COLUMNS`) match the `TOOL_WEIGHT` names; adjust the map if any name differs.

## 3. Figure-side vocabulary

- [ ] 3.1 Update `ssign_lib/functional_vocab.py` `consensus_bucket` + `_KNOWN_CATEGORIES` to the v2 category set; keep machinery → `Apparatus-associated`, unknown/blank → `Other`.
- [ ] 3.2 Grep for any hard-coded old category names in figures/report/docs and update.

## 4. Tests

- [ ] 4.1 Rewrite `tests/unit/test_annotation_consensus.py`: weighted vote, fold down-weight, single-category classify, floors, machinery-by-identity (no `t[1-9]ss` hijack), RTX-by-activity.
- [ ] 4.2 Update `tests/unit/test_functional_vocab.py` for the new `consensus_bucket` category set.
- [ ] 4.3 `pytest tests/unit/ -q` green.

## 5. Golden fixture + validation

- [ ] 5.1 Regenerate the end-to-end golden fixture (`tests/fixtures/golden/`) since `broad_annotation` values change; note the vocabulary change in `REGENERATE.md`.
- [ ] 5.2 Run `docs/development/regrade.py`; confirm fallback-noise = 0, 0 keyword-hijacked effectors, ≤ v1 outside-vocab on the 55-gold with in-vocab calls correct, no single category > ~40% on Xantho. Record the numbers in the change.
- [ ] 5.3 Resolve the design Open Question (autotransporter passenger-function vs structural class) against the regenerated figures before finalizing.
