## Context

`compute_consensus` frequency-counts *every* category matched by *every* tool's free-text description and takes the mode. The audit (`docs/development/consensus_annotation_audit.md`) showed this fails because (a) `classify_description` returns all matches, splitting votes; (b) two of the eight input tools (pLM-BLAST/ECOD, HHpred_PDB) report structural *fold* names, not functions, and a third (InterProScan) mixes fold-superfamily names in; (c) machinery and effectors are conflated via a `t[1-9]ss` keyword; (d) unmatched descriptions mint title-case noise categories that can win the vote. The v2 design was prototyped (`docs/development/v2lib.py`) and measured against v1 on the Xantho 593-substrate panel, a 55-protein gold accuracy sheet, and the full-tier 24-substrate smoke test.

## Goals / Non-Goals

**Goals:**
- Never let a structural fold name outvote a functional annotation.
- Never let machinery/translocator components be labelled with a cargo function, and never let effectors be labelled as machinery.
- Emit an honest `Unclassified` instead of an invented label when no tool has a functional call.
- Cover the effector activity classes that define each secretion system (per the literature), at a useful granularity.
- Keep the pipeline interface unchanged (same input columns, same output column names).

**Non-Goals:**
- A curated gene/family override for garbage-in cases where the tools carry no correct call (Tde1, VasX, apxIIIA). That is a separate follow-on change.
- Changing which tools run, the DB tiers, or any pipeline step.
- Re-scoring or re-deriving secretion-system membership (unchanged).

## Decisions

1. **Weighted vote over frequency count.** Each tool casts one deduped vote; tool weight encodes credibility (functional-name tools 3, domain tools 2, structure-only tools 1). *Alternative considered:* strict precedence (Tier-1 wins outright). Rejected: weighting degrades gracefully when a Tier-1 tool is silent (common on environmental genomes where Swiss-Prot BLASTp returns nothing) while still letting two Tier-2 agreements beat one Tier-3.
2. **Fold-name down-weight as a cross-cutting cap.** A description matching the structural-fold regex ("… superfamily", "…-like domain", "beta-helix/barrel/propeller", "TIM barrel") is capped to weight 1 regardless of source tool. *Alternative:* drop ECOD/PDB entirely. Rejected: they are the only signal for some otherwise-unannotated proteins; capping keeps them as tie-breakers without letting a fold name dominate.
3. **Single specificity-ranked category per description.** Patterns ordered most-specific → generic; first match wins. *Alternative:* keep multi-match but weight by specificity. Rejected as more complex for no measured gain; single-match eliminated the 809 vote-splits directly.
4. **Machinery is evidence-based and competes by weight.** Component-identity regexes (with word boundaries) + translocator context route to `Apparatus-associated`, which competes in the weighted vote (so a translocator the good tools agree on wins) rather than being a last-resort floor. Only `Hypothetical`/`Unclassified` are floors.
5. **Taxonomy from the literature, not ad-hoc.** Two-level set grounded in the T3/T6/T1/T5 effector reviews cited in the audit; RTX dropped as a category (it is a secretion motif spanning multiple functions).

## Risks / Trade-offs

- **Output vocabulary changes** (breaking for anyone parsing `broad_annotation` strings) → the values were never a stable API; figures consume them via `consensus_bucket`, which is updated in lockstep. Document in the change; the golden end-to-end fixture is regenerated.
- **`Unclassified` share rises** (v1 always emitted *a* label) → intended: on Xantho the rise is genuinely-unannotated OM β-barrels / empty-annotation rows. Honest > wrong for a functional figure.
- **Tool weights are hand-set, not learned** → they are justified per-tool by the audit evidence and validated on three datasets; revisit only if a tool's annotation character changes.
- **Residual garbage-in errors remain** → explicitly out of scope; the curated-override follow-on addresses them.

## Migration Plan

- Rewrite in place; no data migration. Regenerate the golden end-to-end fixture (`tests/fixtures/golden/`) since `broad_annotation` values change; document the vocabulary change in the fixture's `REGENERATE.md`.
- Rollback = revert the change; no persisted state.
- Validation gate before merge: `docs/development/regrade.py` must show fallback-noise 0, 0 keyword-hijacked effectors, ≤ v1 outside-vocab on the 55-gold with in-vocab calls correct, and no single category > ~40% on Xantho.

## Open Questions

- Do the run figures want autotransporters split by passenger *function* (espP → Protease) vs the *structural* class? The prototype reports passenger function; confirm this is the desired figure behaviour when regenerating the golden outputs.
