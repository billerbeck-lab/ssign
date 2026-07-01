"""Unit tests for the curated figure generator (headless, Agg backend)."""

import os

import pandas as pd
from _helpers import run_script_main

from ssign_app.scripts.generate_figures import main as gen_main

_TYPES = ["T1SS", "T2SS", "T3SS", "T5aSS", "T6SS"]

# Functional figures (04-07), by-SS-type only, fixed numbers per source.
_FUNCTIONAL = [
    "04_cog_category_by_sstype.png",
    "05_kegg_function_by_sstype.png",
    "06_eggnog_description_by_sstype.png",
    "07_consensus_function_by_sstype.png",
]
_PER_GENOME = [
    "01_secreted_by_genome.png",
    "02_autotransporter_self_detection.png",
    "03_physicochemical.png",
    *_FUNCTIONAL,
]
_PHYSCHEM_COLS = ["aa_length", "gravy", "mw_da", "isoelectric_point", "instability_index", "aromaticity", "charge_ph7"]


def _make_df(n=30, genome="genome_a"):
    rows = []
    for i in range(n):
        t = _TYPES[i % len(_TYPES)]
        rows.append(
            {
                "locus_tag": f"{genome}_{i:03d}",
                "sample_id": genome,
                "tool": ["DLP", "DSE", "DLP,DSE", "T5SS-self"][i % 4],
                "nearby_ss_types": t,
                "dlp_extracellular_prob": round(0.3 + 0.6 * ((i % 7) / 6), 3),
                "outer_membrane_prob": round(0.2 + 0.7 * ((i % 5) / 4), 3),
                "signalp_prediction": "SP(Sec/SPI)" if i % 3 == 0 else "OTHER",
                # Size + physicochemical (ProtParam) columns -> figure 03.
                "aa_length": 120 + 30 * (i % 20),
                "gravy": round(-0.5 + 0.1 * (i % 10), 3),
                "mw_da": 20000 + 2500 * (i % 20),
                "isoelectric_point": round(4.5 + 0.3 * (i % 15), 2),
                "instability_index": round(20 + 2 * (i % 15), 1),
                "aromaticity": round(0.05 + 0.005 * (i % 12), 3),
                "charge_ph7": round(-15 + 2 * (i % 16), 1),
                # Functional sources (real value formats).
                "cog_category": ["U", "MU", "S", ""][i % 4],
                "kegg_ko": ["ko:K10953;ko:K12516", "ko:K15125", "-", ""][i % 4],
                "eggnog_description": ["Type V secretory pathway adhesin", "Hemolysin activator", "-", ""][i % 4],
                "broad_annotation": ["Toxin", "Protease", "Adhesin", "Hypothetical"][i % 4],
            }
        )
    return pd.DataFrame(rows)


def _pngs(d):
    return sorted(f for f in os.listdir(d) if f.endswith(".png")) if os.path.isdir(d) else []


def _run(monkeypatch, df, out, extra=()):
    csv = out.parent / "integrated.csv"
    df.to_csv(csv, index=False)
    run_script_main(
        monkeypatch,
        gen_main,
        ["generate_figures.py", "--master-csvs", str(csv), "--outdir", str(out), "--dpi", "80", *extra],
    )
    return _pngs(out)


def test_per_genome_emits_numbered_set(monkeypatch, tmp_path):
    produced = _run(monkeypatch, _make_df(), tmp_path / "figs")
    assert produced == _PER_GENOME
    for f in produced:
        assert (tmp_path / "figs" / f).stat().st_size > 0
    assert not any(f.startswith("fig") or f.startswith("P") for f in produced)


def test_physicochemical_skipped_when_all_size_columns_absent(monkeypatch, tmp_path):
    df = _make_df().drop(columns=_PHYSCHEM_COLS)
    produced = _run(monkeypatch, df, tmp_path / "figs")
    assert "03_physicochemical.png" not in produced
    assert "01_secreted_by_genome.png" in produced
    for f in _FUNCTIONAL:
        assert f in produced


def test_physicochemical_emitted_with_only_length(monkeypatch, tmp_path):
    # Length lives in the physicochemical figure, so length alone still emits it.
    df = _make_df().drop(columns=[c for c in _PHYSCHEM_COLS if c != "aa_length"])
    produced = _run(monkeypatch, df, tmp_path / "figs")
    assert "03_physicochemical.png" in produced


def test_autotransporter_skipped_without_om_column(monkeypatch, tmp_path):
    df = _make_df().drop(columns=["outer_membrane_prob"])
    produced = _run(monkeypatch, df, tmp_path / "figs")
    assert "02_autotransporter_self_detection.png" not in produced


def test_pooled_mode_two_genomes(monkeypatch, tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    _make_df(genome="genome_a").to_csv(a, index=False)
    _make_df(genome="genome_b").to_csv(b, index=False)
    out = tmp_path / "pool"
    run_script_main(
        monkeypatch,
        gen_main,
        [
            "generate_figures.py",
            "--master-csvs",
            str(a),
            str(b),
            "--outdir",
            str(out),
            "--mode",
            "pooled",
            "--dpi",
            "80",
        ],
    )
    # The curated set over all genomes combined, named 0N_pooled_*. No P0N_ figures.
    expected = sorted(f"{f[:2]}_pooled_{f[3:]}" for f in _PER_GENOME)
    assert _pngs(out) == expected


def test_pooled_mode_single_genome_noop(monkeypatch, tmp_path):
    csv = tmp_path / "one.csv"
    _make_df(genome="only").to_csv(csv, index=False)
    out = tmp_path / "pool1"
    run_script_main(
        monkeypatch,
        gen_main,
        ["generate_figures.py", "--master-csvs", str(csv), "--outdir", str(out), "--mode", "pooled", "--dpi", "80"],
    )
    assert _pngs(out) == []


def test_toggle_skips_named_figures(monkeypatch, tmp_path):
    produced = _run(monkeypatch, _make_df(), tmp_path / "figs", extra=["--no-physicochemical", "--no-func-summary"])
    assert "03_physicochemical.png" not in produced
    assert not any(f in produced for f in _FUNCTIONAL)
    assert "01_secreted_by_genome.png" in produced
    assert "02_autotransporter_self_detection.png" in produced


def test_no_functional_source_columns_skips_block(monkeypatch, tmp_path):
    df = _make_df().drop(columns=["cog_category", "kegg_ko", "eggnog_description", "broad_annotation"])
    produced = _run(monkeypatch, df, tmp_path / "figs")
    assert not any(f in produced for f in _FUNCTIONAL)
    assert "01_secreted_by_genome.png" in produced


def test_functional_numbers_are_fixed_per_source(monkeypatch, tmp_path):
    # Dropping COG removes only 04; KEGG/EggNOG/consensus keep their fixed numbers.
    df = _make_df().drop(columns=["cog_category"])
    produced = _run(monkeypatch, df, tmp_path / "figs")
    assert "04_cog_category_by_sstype.png" not in produced
    assert "05_kegg_function_by_sstype.png" in produced
    assert "07_consensus_function_by_sstype.png" in produced
