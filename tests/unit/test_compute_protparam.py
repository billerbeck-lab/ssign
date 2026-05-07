"""Tests for compute_protparam.py.

Wraps BioPython's ProteinAnalysis to compute MW / pI / GRAVY / instability /
aromaticity / charge@pH7 per substrate. Defensive on input: filters out
non-standard amino acids, skips short sequences (<10 aa) and proteins
missing from the FASTA, swallows per-protein ProtParam exceptions.
"""

import csv
import os
import sys

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta


from _helpers import run_script_main, write_tsv  # noqa: E402
from compute_protparam import main as protparam_main  # noqa: E402

SUBSTRATE_FIELDS = ["locus_tag"]
EXPECTED_OUTPUT_FIELDS = [
    "locus_tag",
    "mw_da",
    "isoelectric_point",
    "gravy",
    "instability_index",
    "aromaticity",
    "charge_ph7",
]


def _run(monkeypatch, tmp_dir, substrate_ids, sequences):
    substrates = write_tsv(
        os.path.join(tmp_dir, "substrates.tsv"),
        SUBSTRATE_FIELDS,
        [{"locus_tag": pid} for pid in substrate_ids],
    )
    proteins = os.path.join(tmp_dir, "proteins.faa")
    write_fasta(sequences, proteins)
    out = os.path.join(tmp_dir, "protparam.csv")
    run_script_main(
        monkeypatch,
        protparam_main,
        [
            "compute_protparam",
            "--substrates",
            substrates,
            "--proteins",
            proteins,
            "--sample",
            "test_sample",
            "--output",
            out,
        ],
    )
    with open(out) as f:
        return list(csv.DictReader(f)), out


# A 30-residue real-ish protein (well above the 10-aa floor)
_REAL_SEQ = "MKTLLLTLLCAFSVAQAVDLPTQEPALGK"


# ---------------------------------------------------------------------------
# Happy path + output schema
# ---------------------------------------------------------------------------


def test_all_six_properties_emitted(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ})
    assert len(rows) == 1
    assert set(rows[0].keys()) == set(EXPECTED_OUTPUT_FIELDS)
    # Values are numeric, populated
    for field in EXPECTED_OUTPUT_FIELDS:
        if field == "locus_tag":
            continue
        assert rows[0][field] != ""
        float(rows[0][field])


def test_output_is_csv_not_tsv(monkeypatch, tmp_dir):
    # Comma-separated output (the script uses default csv.DictWriter delimiter)
    _, out_path = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ})
    with open(out_path) as f:
        first_line = f.readline()
    assert "," in first_line
    assert "\t" not in first_line


# ---------------------------------------------------------------------------
# Non-standard amino-acid handling
# ---------------------------------------------------------------------------


def test_non_standard_aa_stripped(monkeypatch, tmp_dir):
    # Sequence with X (unknown), U (selenocysteine), * (stop) interspersed —
    # ProtParam still works once the stripped sequence has ≥10 residues.
    seq_with_garbage = "MK*TLLLTXXLLCAFSVAQAVDLPTUQEPALGK"
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": seq_with_garbage})
    assert len(rows) == 1


def test_sequence_under_10_clean_aa_skipped(monkeypatch, tmp_dir):
    # 8 standard residues + 6 X's = 8 clean AAs, below the floor → skipped silently.
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": "MKTLLLLLXXXXXX"})
    assert rows == []


def test_lowercase_sequence_handled(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ.lower()})
    assert len(rows) == 1


# ---------------------------------------------------------------------------
# Missing / empty edge cases
# ---------------------------------------------------------------------------


def test_substrate_id_missing_from_fasta_skipped(monkeypatch, tmp_dir):
    # Substrate references "ORPHAN" which has no entry in the proteins FASTA.
    # The wrapper logs and skips rather than crashing.
    rows, _ = _run(monkeypatch, tmp_dir, ["ORPHAN"], {"P1": _REAL_SEQ})
    assert rows == []


def test_empty_substrates_writes_header_only(monkeypatch, tmp_dir):
    _, out = _run(monkeypatch, tmp_dir, [], {"P1": _REAL_SEQ})
    with open(out) as f:
        lines = f.readlines()
    assert len(lines) == 1
    assert lines[0].strip().split(",") == EXPECTED_OUTPUT_FIELDS


def test_partial_substrate_overlap_writes_only_matched(monkeypatch, tmp_dir):
    # Two substrates requested, only one is in the FASTA → only one row emitted
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        ["P1", "P2"],
        {"P1": _REAL_SEQ},
    )
    assert {r["locus_tag"] for r in rows} == {"P1"}


# ---------------------------------------------------------------------------
# Output value sanity (BioPython invariants)
# ---------------------------------------------------------------------------


def test_mw_positive_and_proportional_to_length(monkeypatch, tmp_dir):
    short = _REAL_SEQ
    long_seq = _REAL_SEQ * 4
    rows, _ = _run(
        monkeypatch,
        tmp_dir,
        ["SHORT", "LONG"],
        {"SHORT": short, "LONG": long_seq},
    )
    by_locus = {r["locus_tag"]: r for r in rows}
    short_mw = float(by_locus["SHORT"]["mw_da"])
    long_mw = float(by_locus["LONG"]["mw_da"])
    assert short_mw > 0
    assert long_mw > short_mw  # Longer sequence ⇒ heavier protein


def test_pi_in_biological_range(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ})
    pi = float(rows[0]["isoelectric_point"])
    # Biological pI bounds — anything outside [2, 12] would indicate parser corruption
    assert 2.0 <= pi <= 12.0


def test_aromaticity_in_unit_range(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ})
    aromaticity = float(rows[0]["aromaticity"])
    # Aromaticity is a fraction (F+W+Y / total) — must be in [0, 1]
    assert 0.0 <= aromaticity <= 1.0


def test_output_precision_two_decimals_for_mw(monkeypatch, tmp_dir):
    rows, _ = _run(monkeypatch, tmp_dir, ["P1"], {"P1": _REAL_SEQ})
    # MW is rounded to 2 decimals; pI/GRAVY/etc to 4. Pin the format.
    mw_str = rows[0]["mw_da"]
    if "." in mw_str:
        decimals = mw_str.split(".")[1]
        assert len(decimals) <= 2
