# Golden-output reference

This directory holds **frozen reference outputs** for the end-to-end pipeline
regression test (`tests/integration/test_pipeline_e2e_golden.py`). Each file
is the byte-for-byte output of running ssign on the bundled minimal T5aSS
fixture (`tests/fixtures/Xanthobacter_T5aSS_minimal.gbff`, 20 kb / 9 CDS).

## Why this exists

Unit tests catch logic bugs in individual functions. Property-based tests
catch invariants. Neither catches **integration regressions**: a code change
that is locally correct but quietly alters a downstream column, drops a row,
or shifts a probability by 0.001. For a tool whose output is going into a
publication, that class of regression is exactly what we cannot afford.

The golden-output test runs the full pipeline, then diffs every produced TSV
against the reference here. Any diff fails CI with a clear "this changed"
message. Maintainers approving the change re-generate the references and
commit the new state.

## Layout

```
tests/fixtures/golden/
├── README.md              ← this file
├── REGENERATE.md          ← step-by-step regeneration instructions
└── t5ass_minimal/         ← reference outputs for the T5aSS fixture
    ├── substrates_filtered.tsv
    ├── substrates_unfiltered.tsv
    ├── ss_components.tsv
    ├── valid_systems.tsv
    ├── predictions.tsv
    └── master.csv         ← integrate_annotations final output
```

The `t5ass_minimal/` directory is populated; see plan task 4.i.

## Diff strategy

Some columns are non-deterministic across runs (timestamps, full file paths,
ProtParam float precision artefacts, run UUIDs). The test normalises those
before diffing. Non-determinism beyond the normalised set is treated as a
real regression; investigate before re-freezing.

## Updating the references

Don't edit these files by hand. The regeneration command in `REGENERATE.md`
is the only sanctioned way to update them. Each update should include in
its commit message:

1. Why the diff is intentional (link to the code change that caused it).
2. A summary of which columns / rows changed and by how much.
3. Confirmation that the change is biologically reasonable (not just code-clean).
