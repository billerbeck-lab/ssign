"""Tests for t5ss_handler.py.

T5SS components are self-secreting — the component IS its own substrate.
Pin: only ss_type values starting with "T5" are picked up; non-T5
components are ignored; predictions are looked up by locus_tag and merged
in; the substrate output tags everything with tool="T5SS-self"; the domain
output marks every T5 component as "T5SS-component".
"""

import os
import sys

import pytest


from _helpers import (  # noqa: E402
    PREDICTIONS_FIELDS,
    SS_COMPONENT_FIELDS,
    make_prediction_row,
    make_ss_component_row,
    read_tsv_rows,
    run_script_main,
    write_tsv,
)
from t5ss_handler import main as t5ss_main  # noqa: E402


def _run_t5ss(monkeypatch, tmp_dir, components, predictions):
    ss_path = write_tsv(
        os.path.join(tmp_dir, "ss_components.tsv"),
        SS_COMPONENT_FIELDS,
        components,
    )
    pred_path = write_tsv(
        os.path.join(tmp_dir, "predictions.tsv"),
        PREDICTIONS_FIELDS,
        predictions,
    )
    out_substrates = os.path.join(tmp_dir, "t5ss_substrates.tsv")
    out_domains = os.path.join(tmp_dir, "t5ss_domains.tsv")
    run_script_main(
        monkeypatch,
        t5ss_main,
        [
            "t5ss_handler",
            "--ss-components",
            ss_path,
            "--predictions",
            pred_path,
            "--sample",
            "test_sample",
            "--out-substrates",
            out_substrates,
            "--out-domains",
            out_domains,
        ],
    )
    return read_tsv_rows(out_substrates), read_tsv_rows(out_domains)


@pytest.mark.parametrize("ss_type", ["T5aSS", "T5bSS", "T5cSS"])
def test_every_t5_subtype_recognised(monkeypatch, tmp_dir, ss_type):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", ss_type)],
        [make_prediction_row("AT_1", dlp_ext=0.9)],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    assert substrates[0]["nearby_ss_types"] == ss_type


def test_non_t5_component_ignored(monkeypatch, tmp_dir):
    # T2SS, T1SS, T6SS components must NOT appear as T5SS-self substrates.
    substrates, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [
            make_ss_component_row("T1_1", "T1SS"),
            make_ss_component_row("T2_1", "T2SS"),
            make_ss_component_row("T5_1", "T5aSS"),
            make_ss_component_row("T6_1", "T6SSi"),
        ],
        [make_prediction_row(locus, dlp_ext=0.9) for locus in ["T1_1", "T2_1", "T5_1", "T6_1"]],
    )
    assert {r["locus_tag"] for r in substrates} == {"T5_1"}
    assert {r["locus_tag"] for r in domains} == {"T5_1"}


def test_substrate_tool_is_t5ss_self(monkeypatch, tmp_dir):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [make_prediction_row("AT_1", dlp_ext=0.9)],
    )
    assert substrates[0]["tool"] == "T5SS-self"


def test_prediction_fields_merged_into_substrate(monkeypatch, tmp_dir):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [
            make_prediction_row(
                "AT_1",
                dlp_ext=0.95,
                dse_type="T1SS",
                dse_prob=0.7,
                signalp_pred="OTHER",
                signalp_prob=0.05,
                product="autotransporter passenger",
            )
        ],
    )
    s = substrates[0]
    assert s["dlp_extracellular_prob"] == "0.95"
    assert s["dse_ss_type"] == "T1SS"
    assert s["dse_max_prob"] == "0.7"
    assert s["product"] == "autotransporter passenger"


def test_t5ss_kept_even_when_dlp_low(monkeypatch, tmp_dir):
    # T5aSS median outer-membrane probability is ~0.47; the wrapper does NOT
    # apply DLP thresholds — components are always emitted as substrates.
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [make_prediction_row("AT_1", dlp_ext=0.05)],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    assert substrates[0]["dlp_extracellular_prob"] == "0.05"


def test_missing_prediction_yields_empty_fields(monkeypatch, tmp_dir):
    # A T5 component with no row in predictions.tsv still appears as a
    # substrate; merge fields fall back to empty strings (or 0.0 for dlp).
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    s = substrates[0]
    assert s["dlp_extracellular_prob"] in ("0", "0.0")
    assert s["dse_ss_type"] == ""
    assert s["product"] == ""


def test_domain_output_marks_every_t5_component(monkeypatch, tmp_dir):
    _, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [
            make_ss_component_row("T5a_1", "T5aSS"),
            make_ss_component_row("T5b_1", "T5bSS"),
        ],
        [],
    )
    by_locus = {r["locus_tag"]: r for r in domains}
    assert by_locus["T5a_1"]["domain_group"] == "T5SS-component"
    assert by_locus["T5b_1"]["domain_group"] == "T5SS-component"
    assert by_locus["T5a_1"]["ss_type"] == "T5aSS"


def test_no_t5_components_writes_empty_outputs(monkeypatch, tmp_dir):
    substrates, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("T2_1", "T2SS")],
        [make_prediction_row("T2_1", dlp_ext=0.9)],
    )
    assert substrates == []
    assert domains == []


def test_invalid_dlp_prob_falls_back_to_zero(monkeypatch, tmp_dir):
    # If the predictions row has a non-numeric extracellular_prob, the
    # wrapper logs and falls back to 0.0 rather than crashing.
    bad_pred = make_prediction_row("AT_1")
    bad_pred["dlp_extracellular_prob"] = "not-a-number"
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [bad_pred],
    )
    assert substrates[0]["dlp_extracellular_prob"] in ("0", "0.0")
