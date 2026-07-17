"""Cross-genome ortholog wiring + conservation figure.

Two surfaces are covered without needing BLAST+ installed:

* ``run_cross_genome_orthologs`` (runner) — the all-vs-all BLASTp call is
  monkeypatched to emit the canned ``run_ortholog_grouping`` output, so the test
  exercises the pooling, the sample_id / n_genomes / genomes augmentation, and
  the xg_ortholog_group merge-back.
* ``run_cross_genome_ortholog_figure.render`` — pure matplotlib, exercised on a
  synthetic groups CSV + integrated CSVs.
"""

import csv
import os

import pandas as pd
import run_cross_genome_ortholog_figure as figmod

from ssign_app.core import runner as runner_mod


def _write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def _write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _integrated(path, sample_id, rows):
    """Minimal integrated CSV: locus_tag + sample_id + annotation + ss type."""
    _write_csv(
        path,
        ["locus_tag", "sample_id", "nearby_ss_types", "broad_consensus_annotation"],
        [
            {"locus_tag": lt, "sample_id": sample_id, "nearby_ss_types": ss, "broad_consensus_annotation": ann}
            for lt, ss, ann in rows
        ],
    )


def _fake_ortholog_grouping(script, args, timeout=None):
    """Stand in for run_ortholog_grouping.py: read the (prefixed-id) combined
    FASTA the runner wrote and single-linkage-cluster it by identical sequence,
    emitting the per-protein + groups CSVs with the same prefixed ids."""
    argd = dict(zip(args[::2], args[1::2]))
    seqs: dict[str, str] = {}
    cur = None
    with open(argd["--substrates-fasta"]) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = ""
            elif cur:
                seqs[cur] += line
    by_seq: dict[str, list[str]] = {}
    for pid, s in seqs.items():
        by_seq.setdefault(s, []).append(pid)
    groups = sorted(by_seq.values(), key=lambda g: (-len(g), g))
    per_rows, group_rows = [], []
    for i, members in enumerate(groups, 1):
        og = f"OG_{i:03d}"
        for pid in members:
            per_rows.append(
                {"locus_tag": pid, "ortholog_group": og, "og_n_members": len(members), "og_mean_pident": 100.0}
            )
        group_rows.append(
            {
                "ortholog_group": og,
                "n_members": len(members),
                "members": ";".join(sorted(members)),
                "mean_pident": 100.0,
            }
        )
    _write_csv(argd["--output"], ["locus_tag", "ortholog_group", "og_n_members", "og_mean_pident"], per_rows)
    _write_csv(argd["--output-groups"], ["ortholog_group", "n_members", "members", "mean_pident"], group_rows)
    return 0, "", ""


# --------------------------------------------------------------------------- #
# run_cross_genome_orthologs (runner)
# --------------------------------------------------------------------------- #
def test_run_cross_genome_orthologs_augments_and_merges(tmp_dir, monkeypatch):
    a_int = os.path.join(tmp_dir, "A_integrated.csv")
    b_int = os.path.join(tmp_dir, "B_integrated.csv")
    a_faa = os.path.join(tmp_dir, "A.faa")
    b_faa = os.path.join(tmp_dir, "B.faa")
    _integrated(a_int, "A", [("a1", "T2SS", "Protease"), ("a2", "T1SS", "Adhesin")])
    _integrated(b_int, "B", [("b1", "T2SS", "Protease")])
    _write_fasta(a_faa, {"a1": "MKAA", "a2": "MKBB", "zz": "MKZZ"})
    _write_fasta(b_faa, {"b1": "MKAA"})

    monkeypatch.setattr(runner_mod, "run_script", _fake_ortholog_grouping)

    res = runner_mod.run_cross_genome_orthologs(
        [("A", a_int, a_faa), ("B", b_int, b_faa)],
        output_dir=tmp_dir,
    )

    assert res["n_proteins"] == 3  # a1, a2, b1 (zz is not a substrate)
    assert res["n_groups"] == 2
    assert res["genomes_updated"] == 2

    # per-protein CSV gained sample_id
    og = pd.read_csv(res["orthologs_csv"]).set_index("locus_tag")
    assert og.loc["a1", "sample_id"] == "A"
    assert og.loc["b1", "sample_id"] == "B"

    # group CSV gained n_genomes + genomes (unique genomes, paralogs deduped)
    grp = pd.read_csv(res["ortholog_groups_csv"]).set_index("ortholog_group")
    assert grp.loc["OG_001", "n_genomes"] == 2
    assert grp.loc["OG_001", "genomes"] == "A;B"
    assert grp.loc["OG_002", "n_genomes"] == 1
    assert grp.loc["OG_002", "genomes"] == "A"

    # xg_ortholog_group merged into each integrated CSV
    a_df = pd.read_csv(a_int).set_index("locus_tag")
    assert a_df.loc["a1", "xg_ortholog_group"] == "OG_001"
    assert a_df.loc["a2", "xg_ortholog_group"] == "OG_002"
    b_df = pd.read_csv(b_int).set_index("locus_tag")
    assert b_df.loc["b1", "xg_ortholog_group"] == "OG_001"


def test_shared_locus_tag_across_genomes_counts_both(tmp_dir, monkeypatch):
    """A locus_tag identical in two genomes (e.g. a shared RefSeq WP_ accession)
    must be kept as two records and counted as 2 genomes, not deduped to 1."""
    a_int = os.path.join(tmp_dir, "A_integrated.csv")
    b_int = os.path.join(tmp_dir, "B_integrated.csv")
    a_faa = os.path.join(tmp_dir, "A.faa")
    b_faa = os.path.join(tmp_dir, "B.faa")
    # WP_0001 is the SAME id in both genomes with the SAME sequence (conserved).
    _integrated(a_int, "A", [("WP_0001", "T2SS", "Protease"), ("a_uniq", "T1SS", "Adhesin")])
    _integrated(b_int, "B", [("WP_0001", "T2SS", "Protease"), ("b_uniq", "T1SS", "Adhesin")])
    _write_fasta(a_faa, {"WP_0001": "MKAA", "a_uniq": "MKCC"})
    _write_fasta(b_faa, {"WP_0001": "MKAA", "b_uniq": "MKDD"})

    monkeypatch.setattr(runner_mod, "run_script", _fake_ortholog_grouping)
    res = runner_mod.run_cross_genome_orthologs([("A", a_int, a_faa), ("B", b_int, b_faa)], output_dir=tmp_dir)

    assert res["n_proteins"] == 4  # both WP_0001 copies kept
    grp = pd.read_csv(res["ortholog_groups_csv"])
    shared = grp[grp["n_genomes"] == 2]
    assert len(shared) == 1
    assert shared.iloc[0]["genomes"] == "A;B"

    # Each genome's own WP_0001 row gets the shared group merged in.
    a_df = pd.read_csv(a_int).set_index("locus_tag")
    b_df = pd.read_csv(b_int).set_index("locus_tag")
    assert a_df.loc["WP_0001", "xg_ortholog_group"] == shared.iloc[0]["ortholog_group"]
    assert b_df.loc["WP_0001", "xg_ortholog_group"] == shared.iloc[0]["ortholog_group"]


def test_run_cross_genome_orthologs_skips_below_two(tmp_dir, monkeypatch):
    a_int = os.path.join(tmp_dir, "A_integrated.csv")
    a_faa = os.path.join(tmp_dir, "A.faa")
    _integrated(a_int, "A", [("a1", "T2SS", "Protease")])
    _write_fasta(a_faa, {"a1": "MKAA"})

    called = []
    monkeypatch.setattr(runner_mod, "run_script", lambda *a, **k: called.append(1) or (0, "", ""))

    res = runner_mod.run_cross_genome_orthologs([("A", a_int, a_faa)], output_dir=tmp_dir)
    assert res["n_proteins"] == 1
    assert res["n_groups"] == 0
    assert not called  # never reached BLAST


# --------------------------------------------------------------------------- #
# run_cross_genome_ortholog_figure
# --------------------------------------------------------------------------- #
def test_primary_ss_type():
    assert figmod._primary_ss_type("T2SS") == "T2SS"
    assert figmod._primary_ss_type("T2SS;T6SS") == "T2SS"
    assert figmod._primary_ss_type("T5aSS,T1SS") == "T5aSS"
    assert figmod._primary_ss_type("") == ""
    assert figmod._primary_ss_type("nan") == ""
    assert figmod._primary_ss_type(None) == ""


def test_build_protein_maps_prefers_consensus(tmp_dir):
    p = os.path.join(tmp_dir, "x_integrated.csv")
    _write_csv(
        p,
        ["locus_tag", "sample_id", "nearby_ss_types", "broad_consensus_annotation", "broad_annotation"],
        [
            {
                "locus_tag": "x1",
                "sample_id": "G",
                "nearby_ss_types": "T2SS",
                "broad_consensus_annotation": "Consensus",
                "broad_annotation": "Fallback",
            }
        ],
    )
    ann, ss, genome = figmod.build_protein_maps([p])
    assert ann["x1"] == "Consensus"
    assert ss["x1"] == "T2SS"
    assert genome["x1"] == "G"


def test_render_figure_synthetic(tmp_dir):
    groups = os.path.join(tmp_dir, "cross_genome_ortholog_groups.csv")
    _write_csv(
        groups,
        ["ortholog_group", "n_members", "n_genomes", "genomes", "members", "mean_pident"],
        [
            {
                "ortholog_group": "OG_001",
                "n_members": 3,
                "n_genomes": 3,
                "genomes": "A;B;C",
                "members": "a1;b1;c1",
                "mean_pident": 62.0,
            },
            {
                "ortholog_group": "OG_002",
                "n_members": 2,
                "n_genomes": 2,
                "genomes": "A;B",
                "members": "a2;b2",
                "mean_pident": 71.0,
            },
            {
                "ortholog_group": "OG_003",
                "n_members": 1,
                "n_genomes": 1,
                "genomes": "A",
                "members": "a3",
                "mean_pident": 100.0,
            },
        ],
    )
    integrated = []
    for sid, rows in (
        ("A", [("a1", "T2SS", "Protease"), ("a2", "T5aSS", "Adhesin"), ("a3", "T1SS", "Toxin")]),
        ("B", [("b1", "T2SS", "Protease"), ("b2", "T5aSS", "Adhesin")]),
        ("C", [("c1", "T2SS", "Protease")]),
    ):
        pth = os.path.join(tmp_dir, f"{sid}_integrated.csv")
        _integrated(pth, sid, rows)
        integrated.append(pth)

    out = os.path.join(tmp_dir, "07_cross_genome_orthologs.png")
    assert figmod.render(groups, integrated, out) is True
    assert os.path.exists(out) and os.path.getsize(out) > 5000


def test_render_figure_empty_returns_false(tmp_dir):
    groups = os.path.join(tmp_dir, "empty.csv")
    _write_csv(groups, ["ortholog_group", "n_members", "n_genomes", "genomes", "members", "mean_pident"], [])
    out = os.path.join(tmp_dir, "none.png")
    assert figmod.render(groups, [], out) is False
    assert not os.path.exists(out)
