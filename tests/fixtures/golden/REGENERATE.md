# Regenerating golden-output references

When an intentional pipeline change alters one or more outputs, follow these
steps to refresh the frozen reference files.

## 1. Confirm the diff is intentional

Run the golden-output test and read the diff carefully:

```bash
SSIGN_RUN_FULL_PIPELINE=1 \
SSIGN_DEEPLOCPRO_PATH=/path/to/dir/containing/deeplocpro \
pytest -m integration tests/integration/test_pipeline_e2e_golden.py -v
```

If the diff is unexpected, **stop and investigate** — it may be a real
regression. The whole point of golden outputs is to make accidental
behaviour drift visible.

## 2. Run the pipeline on the minimal T5aSS fixture

```bash
SSIGN_RUN_FULL_PIPELINE=1 \
SSIGN_DEEPLOCPRO_PATH=/path/to/dir/containing/deeplocpro \
SSIGN_GOLDEN_REGEN_DIR=/tmp/ssign_golden_regen \
pytest -m integration tests/integration/test_pipeline_e2e_golden.py::test_pipeline_e2e_matches_golden_outputs -v
```

When `SSIGN_GOLDEN_REGEN_DIR` is set, the test runs the pipeline into that
directory instead of using a tempdir, so the produced files survive after
the test exits and you can copy them into the reference directory.

## 3. Copy the produced outputs into the reference directory

The test runs the pipeline with `use_input_annotations=True` (skip Bakta),
DeepLocPro local (offline), and every other prediction/annotation tool
disabled. The user-facing outputs land in `$SSIGN_GOLDEN_REGEN_DIR`
(`$RUN` below) and the intermediate TSVs land in a separate work directory
(`$WORK`). On a diff failure the test prints both paths under
`runner.work_dir=` and `outdir=` — copy them into shell vars and run:

```bash
DEST=tests/fixtures/golden/t5ass_minimal
RUN=/tmp/ssign_golden_regen
WORK=...   # runner.work_dir from the pytest output

# User-facing outputs (outdir):
cp $RUN/t5ass_minimal_results.csv      $DEST/
cp $RUN/t5ass_minimal_results_raw.csv  $DEST/
cp $RUN/t5ass_minimal_summary.txt      $DEST/

# Intermediate work-dir TSV/CSVs that capture each pipeline phase:
for name in gene_info gene_order valid_systems ss_components \
            deeplocpro predictions substrates substrates_filtered \
            t5ss_substrates integrated enrichment_fisher; do
    src=$(ls $WORK/t5ass_minimal_${name}.* 2>/dev/null | head -1)
    [ -n "$src" ] && cp "$src" $DEST/
done
```

## 4. Re-run the test to confirm green

```bash
SSIGN_RUN_FULL_PIPELINE=1 \
SSIGN_DEEPLOCPRO_PATH=/path/to/dir \
pytest -m integration tests/integration/test_pipeline_e2e_golden.py -v
```

## 5. Commit with a clear rationale

The commit message must include:

- **Why** the diff is intentional (link the upstream code change).
- **What** changed (column names? row counts? probability values?).
- **Biological sanity check** — does the new output still call the known
  T5aSS substrate (BIMENO_04457) correctly?

A diff in the golden-output reference without that justification is grounds
for reviewer rejection.
