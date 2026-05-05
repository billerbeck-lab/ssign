# Regenerating golden-output references

When an intentional pipeline change alters one or more outputs, follow these
steps to refresh the frozen reference files.

## 1. Confirm the diff is intentional

Run the golden-output test and read the diff carefully:

```bash
SSIGN_RUN_FULL_PIPELINE=1 pytest -m integration \
    tests/integration/test_pipeline_e2e_golden.py -v
```

If the diff is unexpected, **stop and investigate** — it may be a real
regression. The whole point of golden outputs is to make accidental
behaviour drift visible.

## 2. Run the pipeline on the minimal T5aSS fixture

```bash
mkdir -p /tmp/ssign_golden_regen
ssign run \
    tests/fixtures/Xanthobacter_T5aSS_minimal.gbff \
    --outdir /tmp/ssign_golden_regen \
    --sample-id t5ass_minimal \
    --use-input-annotations
```

`--use-input-annotations` skips Bakta re-annotation (the input GenBank
already has CDS annotations) so the test remains fast and reproducible
without a Bakta install.

## 3. Copy the produced outputs into the reference directory

```bash
cp /tmp/ssign_golden_regen/substrates_filtered.tsv  tests/fixtures/golden/t5ass_minimal/
cp /tmp/ssign_golden_regen/substrates_unfiltered.tsv tests/fixtures/golden/t5ass_minimal/
cp /tmp/ssign_golden_regen/ss_components.tsv        tests/fixtures/golden/t5ass_minimal/
cp /tmp/ssign_golden_regen/valid_systems.tsv        tests/fixtures/golden/t5ass_minimal/
cp /tmp/ssign_golden_regen/predictions.tsv          tests/fixtures/golden/t5ass_minimal/
cp /tmp/ssign_golden_regen/master.csv               tests/fixtures/golden/t5ass_minimal/
```

## 4. Re-run the test to confirm green

```bash
SSIGN_RUN_FULL_PIPELINE=1 pytest -m integration \
    tests/integration/test_pipeline_e2e_golden.py -v
```

## 5. Commit with a clear rationale

The commit message must include:

- **Why** the diff is intentional (link the upstream code change).
- **What** changed (column names? row counts? probability values?).
- **Biological sanity check** — does the new output still call the known
  T5aSS substrate (BIMENO_04457) correctly?

A diff in the golden-output reference without that justification is grounds
for reviewer rejection.
