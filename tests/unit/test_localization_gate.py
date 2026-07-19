"""Unit tests for the literature-derived localization-correctness gate.

Covers the helper in ``ssign_lib.localization_gate``: TSV loading,
per-system scoring, and ss_type aggregation. The shipped TSV under
``ssign_lib/data/`` is also smoke-tested for shape and DLP-vocab
validity since the gate is only as good as the table behind it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ssign_app.scripts.ssign_lib.localization_gate import (
    SystemVerdict,
    aggregate_failed_ss_types,
    evaluate_system,
    load_component_localizations,
)


def _write_table(tmp_dir, rows: list[str]) -> str:
    header = "ss_type\ttxsscan_gene_name\tprotein_alias\tacceptable_localizations\tcitation_pmid\tcitation_doi\tnotes\n"
    path = Path(tmp_dir) / "rules.tsv"
    path.write_text("# header comment\n" + header + "\n".join(rows) + "\n")
    return str(path)


def _row(loc: str, prob: float = 0.9) -> dict[str, str]:
    """Build a fake cross_validate row with the prob in one of the 5 classes.

    Column names match the production cross_validate_predictions.py
    output exactly. Note the asymmetric ``dlp_`` prefix on the
    extracellular column (historical: see cross_validate FIELDNAMES).
    Building rows with the real names means that if a future schema
    rename drops the prefix from one but not the other, this test
    file fails loud instead of silently scoring with 0-probability.
    """
    base = {
        "dlp_extracellular_prob": "0",
        "periplasmic_prob": "0",
        "outer_membrane_prob": "0",
        "cytoplasmic_prob": "0",
        "cytoplasmic_membrane_prob": "0",
        "predicted_localization": loc,
    }
    col_for_label = {
        "Extracellular": "dlp_extracellular_prob",
        "Periplasmic": "periplasmic_prob",
        "Outer Membrane": "outer_membrane_prob",
        "Cytoplasmic": "cytoplasmic_prob",
        "Cytoplasmic Membrane": "cytoplasmic_membrane_prob",
    }
    base[col_for_label[loc]] = str(prob)
    return base


class TestLoadComponentLocalizations:
    def test_parses_basic_rows(self, tmp_dir):
        path = _write_table(
            tmp_dir,
            [
                "T1SS\tT1SS_abc\tHlyB\tCytoplasmic Membrane\t1\t10/1\tnote",
                "T1SS\tT1SS_omf\tTolC\tOuter Membrane\t2\t10/2\tnote",
            ],
        )
        rules = load_component_localizations(path)
        assert rules[("T1SS", "T1SS_abc")] == frozenset({"Cytoplasmic Membrane"})
        assert rules[("T1SS", "T1SS_omf")] == frozenset({"Outer Membrane"})

    def test_multi_compartment_pipe_splits(self, tmp_dir):
        path = _write_table(
            tmp_dir,
            ["T2SS\tT2SS_gspG\tGspG\tCytoplasmic Membrane|Periplasmic\t1\t10/1\tnote"],
        )
        rules = load_component_localizations(path)
        assert rules[("T2SS", "T2SS_gspG")] == frozenset({"Cytoplasmic Membrane", "Periplasmic"})

    def test_empty_acceptable_is_sentinel_for_skip(self, tmp_dir):
        path = _write_table(
            tmp_dir,
            ["pT4SSi\tT4SS_I_traQ\tDotE\t\t\t\tno evidence found"],
        )
        rules = load_component_localizations(path)
        assert rules[("pT4SSi", "T4SS_I_traQ")] == frozenset()

    def test_skips_comments_and_blank_lines(self, tmp_dir):
        path = Path(tmp_dir) / "rules.tsv"
        path.write_text(
            "# header notes\n"
            "# more notes\n"
            "\n"
            "ss_type\ttxsscan_gene_name\tprotein_alias\tacceptable_localizations\tcitation_pmid\tcitation_doi\tnotes\n"
            "T1SS\tT1SS_abc\tHlyB\tCytoplasmic Membrane\t1\t10/1\tnote\n"
        )
        rules = load_component_localizations(str(path))
        assert rules == {("T1SS", "T1SS_abc"): frozenset({"Cytoplasmic Membrane"})}


class TestEvaluateSystem:
    def setup_method(self) -> None:
        self.rules = {
            ("T1SS", "T1SS_abc"): frozenset({"Cytoplasmic Membrane"}),
            ("T1SS", "T1SS_mfp"): frozenset({"Cytoplasmic Membrane"}),
            ("T1SS", "T1SS_omf"): frozenset({"Outer Membrane"}),
        }

    def test_all_components_correct_passes(self):
        preds = {
            "LOC1": _row("Cytoplasmic Membrane"),
            "LOC2": _row("Cytoplasmic Membrane"),
            "LOC3": _row("Outer Membrane"),
        }
        components = [("LOC1", "T1SS_abc"), ("LOC2", "T1SS_mfp"), ("LOC3", "T1SS_omf")]
        v = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert v.n_scored == 3
        assert v.n_correct == 3
        assert v.fraction_correct == 1.0
        assert v.passed

    def test_two_of_three_correct_passes_at_default_threshold(self):
        # 2/3 = 0.667 < 0.8, but >= 0.6 (a more permissive threshold)
        preds = {
            "LOC1": _row("Cytoplasmic Membrane"),
            "LOC2": _row("Cytoplasmic Membrane"),
            "LOC3": _row("Extracellular"),  # wrong
        }
        components = [("LOC1", "T1SS_abc"), ("LOC2", "T1SS_mfp"), ("LOC3", "T1SS_omf")]
        v_strict = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert not v_strict.passed
        v_permissive = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.6)
        assert v_permissive.passed

    def test_low_confidence_components_excluded_from_denominator(self):
        # LOC1 is wrong but DLP confidence is low — should be skipped entirely.
        # The system then has 2/2 confidently-correct components → passes.
        preds = {
            "LOC1": _row("Extracellular", prob=0.5),  # low confidence, excluded
            "LOC2": _row("Cytoplasmic Membrane"),
            "LOC3": _row("Outer Membrane"),
        }
        components = [("LOC1", "T1SS_abc"), ("LOC2", "T1SS_mfp"), ("LOC3", "T1SS_omf")]
        v = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert v.n_scored == 2  # LOC1 excluded
        assert v.n_correct == 2
        assert v.passed

    def test_no_rule_components_skipped(self):
        # Component with no rule in the table → skipped entirely
        preds = {
            "LOC1": _row("Cytoplasmic Membrane"),
            "LOC2": _row("Extracellular"),
            "LOC3": _row("Outer Membrane"),
        }
        components = [
            ("LOC1", "T1SS_abc"),
            ("LOC2", "T1SS_unknown"),  # no rule
            ("LOC3", "T1SS_omf"),
        ]
        v = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert v.n_scored == 2  # LOC2 skipped (no rule)
        assert v.n_correct == 2

    def test_empty_rule_set_skips_component(self):
        rules = {
            ("T1SS", "T1SS_abc"): frozenset({"Cytoplasmic Membrane"}),
            ("T1SS", "T1SS_mystery"): frozenset(),  # explicit sentinel
        }
        preds = {"LOC1": _row("Cytoplasmic Membrane"), "LOC2": _row("Cytoplasmic")}
        components = [("LOC1", "T1SS_abc"), ("LOC2", "T1SS_mystery")]
        v = evaluate_system("T1SS_1", "T1SS", components, preds, rules, 0.8, 0.8)
        assert v.n_scored == 1
        assert v.n_correct == 1
        assert v.passed

    def test_zero_scored_components_fails_open(self):
        # All components have no rule → n_scored == 0 → pass (don't drop on no evidence)
        preds = {"LOC1": _row("Cytoplasmic"), "LOC2": _row("Cytoplasmic")}
        components = [("LOC1", "unknown_gene_1"), ("LOC2", "unknown_gene_2")]
        v = evaluate_system("X_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert v.n_scored == 0
        assert v.passed  # fail-open

    def test_missing_prediction_skips_component(self):
        # Component is in ss_components but not in DLP output (e.g. extract failed)
        preds = {"LOC1": _row("Cytoplasmic Membrane")}
        components = [("LOC1", "T1SS_abc"), ("LOC2", "T1SS_mfp")]
        v = evaluate_system("T1SS_1", "T1SS", components, preds, self.rules, 0.8, 0.8)
        assert v.n_scored == 1  # only LOC1 had a prediction

    def test_multi_compartment_acceptable_match(self):
        rules = {("T2SS", "T2SS_gspG"): frozenset({"Cytoplasmic Membrane", "Periplasmic"})}
        preds = {"LOC1": _row("Periplasmic"), "LOC2": _row("Cytoplasmic Membrane")}
        for locus in ("LOC1", "LOC2"):
            v = evaluate_system("T2SS_1", "T2SS", [(locus, "T2SS_gspG")], preds, rules, 0.8, 0.8)
            assert v.n_correct == 1


class TestAggregateFailedSsTypes:
    def test_all_systems_pass_no_failed_types(self):
        verdicts = [
            SystemVerdict("T1SS_1", "T1SS", 3, 3, 1.0, True),
            SystemVerdict("T2SS_1", "T2SS", 5, 5, 1.0, True),
        ]
        assert aggregate_failed_ss_types(verdicts) == set()

    def test_all_systems_of_type_fail_marks_type_failed(self):
        verdicts = [
            SystemVerdict("T1SS_1", "T1SS", 3, 0, 0.0, False),
            SystemVerdict("T1SS_2", "T1SS", 3, 1, 0.33, False),
            SystemVerdict("T2SS_1", "T2SS", 5, 5, 1.0, True),
        ]
        assert aggregate_failed_ss_types(verdicts) == {"T1SS"}

    def test_partial_failure_within_type_does_not_fail_the_type(self):
        # 2 of 3 T1SS systems failed, but 1 passed → keep T1SS substrates
        verdicts = [
            SystemVerdict("T1SS_1", "T1SS", 3, 0, 0.0, False),
            SystemVerdict("T1SS_2", "T1SS", 3, 1, 0.33, False),
            SystemVerdict("T1SS_3", "T1SS", 3, 3, 1.0, True),  # one passes
        ]
        assert aggregate_failed_ss_types(verdicts) == set()


class TestSchemaCompatibility:
    """Pin the gate to the cross_validate predictions TSV column names.

    Three independent simplify-review agents (2026-06-06) flagged a
    silent bug where the gate read ``extracellular_prob`` but
    cross_validate writes ``dlp_extracellular_prob``. The bug was
    masked because the test fixture used the same wrong name. Tests
    here use the real production FIELDNAMES so a future rename in
    either side fails the suite loud.
    """

    def test_gate_reads_real_cross_validate_extracellular_column(self):
        # Build a row exactly as cross_validate_predictions.py would emit it
        # for a protein DLP called Extracellular at 0.95 confidence.
        from ssign_app.scripts.cross_validate_predictions import FIELDNAMES

        assert "dlp_extracellular_prob" in FIELDNAMES, (
            "cross_validate dropped dlp_extracellular_prob; update _PROB_COLS in "
            "ssign_lib/localization_gate.py to match."
        )

        rules = {("T5aSS", "T5aSS_PF03797"): frozenset({"Extracellular", "Outer Membrane"})}
        # Real-schema row: extracellular at 0.95, others 0. If the gate were
        # reading the wrong column name, this would score as low-confidence
        # (max_prob ~ 0) and the component would be excluded → n_scored == 0
        # → fail-open. With the right column, it gets scored and matches.
        row = {f: "" for f in FIELDNAMES}
        row["predicted_localization"] = "Extracellular"
        row["dlp_extracellular_prob"] = "0.95"
        verdict = evaluate_system("T5aSS_1", "T5aSS", [("LOC1", "T5aSS_PF03797")], {"LOC1": row}, rules, 0.8, 0.8)
        assert verdict.n_scored == 1, "Gate missed the extracellular column — column-name drift?"
        assert verdict.n_correct == 1
        assert verdict.passed


class TestShippedTable:
    """Smoke tests against the real component_localizations.tsv that ships."""

    @pytest.fixture
    def rules(self) -> dict:
        return load_component_localizations()

    def test_file_loads_with_at_least_one_rule(self, rules):
        assert len(rules) > 50, f"Shipped table parsed only {len(rules)} rules; expected ~90"

    def test_all_acceptable_values_use_dlp_vocab(self, rules):
        # Every non-empty acceptable_localizations entry must be one of DLP's
        # 5 class strings — anything else means a typo or a stale label.
        dlp_vocab = {"Cytoplasmic", "Cytoplasmic Membrane", "Periplasmic", "Outer Membrane", "Extracellular"}
        offenders = []
        for (ss_type, gene), accept in rules.items():
            for label in accept:
                if label not in dlp_vocab:
                    offenders.append(f"({ss_type}, {gene}): {label!r}")
        assert not offenders, "Non-DLP-vocab labels in shipped table:\n  " + "\n  ".join(offenders)

    def test_some_rows_use_empty_sentinel(self, rules):
        # The TSV is supposed to carry "no evidence found" rows with empty
        # acceptable_localizations. If that breaks (e.g. someone normalises
        # the file and drops empties), the gate stops being honest about
        # unknowns and starts mis-scoring them.
        assert any(not accept for accept in rules.values()), "No empty-sentinel rows; check #58a notes"
