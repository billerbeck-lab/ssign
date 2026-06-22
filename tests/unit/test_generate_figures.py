"""Unit tests for the curated figure generator (headless, Agg backend)."""

import os

import pandas as pd
from _helpers import run_script_main

from ssign_app.scripts.generate_figures import main as gen_main

_TYPES = ["T1SS", "T2SS", "T3SS", "T5aSS", "T6SS"]
_PER_GENOME = [
    "01_substrates_per_type.png",
    "02_secretion_evidence.png",
    "03_localization_confidence.png",
    "04_signalp_by_type.png",
    "05_tool_coverage.png",
    "06_protein_length.png",
    "07_functional_categories.png",
]


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
                "signalp_prediction": "SP(Sec/SPI)" if i % 3 == 0 else "OTHER",
                "aa_length": 120 + 30 * (i % 20),
                "broad_annotation": ["toxin", "protease", "adhesin", "hypothetical protein"][i % 4],
                "blastp_hit_description": "some hit" if i % 2 == 0 else None,
            }
        )
    return pd.DataFrame(rows)


def _pngs(d):
    return sorted(f for f in os.listdir(d) if f.endswith(".png")) if os.path.isdir(d) else []


def test_per_genome_emits_numbered_set(monkeypatch, tmp_path):
    csv = tmp_path / "integrated.csv"
    _make_df().to_csv(csv, index=False)
    out = tmp_path / "figs"
    run_script_main(
        monkeypatch, gen_main, ["generate_figures.py", "--master-csvs", str(csv), "--outdir", str(out), "--dpi", "80"]
    )

    produced = _pngs(out)
    assert produced == _PER_GENOME
    for f in produced:
        assert (out / f).stat().st_size > 0
    # No legacy figN_ names.
    assert not any(f.startswith("fig") for f in produced)


def test_missing_column_skips_only_its_figure(monkeypatch, tmp_path):
    df = _make_df().drop(columns=["signalp_prediction"])  # only fig 04 depends on it
    csv = tmp_path / "integrated.csv"
    df.to_csv(csv, index=False)
    out = tmp_path / "figs"
    run_script_main(
        monkeypatch, gen_main, ["generate_figures.py", "--master-csvs", str(csv), "--outdir", str(out), "--dpi", "80"]
    )

    produced = _pngs(out)
    assert "04_signalp_by_type.png" not in produced
    # The others still render (05 has blastp coverage independent of signalp).
    for f in [
        "01_substrates_per_type.png",
        "02_secretion_evidence.png",
        "05_tool_coverage.png",
        "07_functional_categories.png",
    ]:
        assert f in produced


def test_physicochemical_off_by_default(monkeypatch, tmp_path):
    csv = tmp_path / "integrated.csv"
    _make_df().to_csv(csv, index=False)
    out = tmp_path / "figs"
    run_script_main(
        monkeypatch, gen_main, ["generate_figures.py", "--master-csvs", str(csv), "--outdir", str(out), "--dpi", "80"]
    )
    assert not any("physicochemical" in f for f in _pngs(out))


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
    produced = _pngs(out)
    assert produced == ["P01_substrates_per_genome.png", "P02_sstype_by_genome.png", "P03_evidence_basis.png"]


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


def test_toggle_skips_named_figure(monkeypatch, tmp_path):
    csv = tmp_path / "integrated.csv"
    _make_df().to_csv(csv, index=False)
    out = tmp_path / "figs"
    run_script_main(
        monkeypatch,
        gen_main,
        [
            "generate_figures.py",
            "--master-csvs",
            str(csv),
            "--outdir",
            str(out),
            "--dpi",
            "80",
            "--no-length",
            "--no-signalp",
        ],
    )
    produced = _pngs(out)
    assert "06_protein_length.png" not in produced
    assert "04_signalp_by_type.png" not in produced
    assert "01_substrates_per_type.png" in produced
