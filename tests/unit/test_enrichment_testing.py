"""Tests for enrichment_testing.py (circular-shift permutation test).

Covers the pure functions (type labels, positivity rules incl. autotransporter
self-detection, the FFT rotation_counts against brute force, window/self masks,
single_test, per-type aggregation with the DSE-T3SS skip) and the CLI driver
end-to-end against synthetic whole-genome prediction TSVs.
"""

import csv
import os

import numpy as np
import pytest
from _helpers import run_script_main, write_tsv
from enrichment_testing import (
    bh_fdr,
    binom_pvalue,
    broad_type,
    components_by_display_type,
    display_type,
    is_dlp_positive,
    is_dlp_self_positive,
    is_dse_positive,
    load_systems,
    positivity_vectors,
    rotation_counts,
    run_enrichment,
    self_mask,
    single_test,
    window_mask,
)
from enrichment_testing import (
    main as enrichment_main,
)


class TestTypeLabels:
    @pytest.mark.parametrize(
        "raw,expected",
        [("T1SS", "T1SS"), ("pT4SSt", "T4SS"), ("T6SSi", "T6SS"), ("T6SSii", "T6SS"), ("Flagellum", "Flagellum")],
    )
    def test_broad_type(self, raw, expected):
        assert broad_type(raw) == expected

    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("T5aSS", "T5aSS"),  # subtypes kept distinct (autotransporters tested differently)
            ("T5bSS", "T5bSS"),
            ("T5cSS", "T5cSS"),
            ("T6SSi", "T6SS"),  # but T6/T4 subtypes still collapse
            ("pT4SSt", "T4SS"),
            ("T1SS", "T1SS"),
        ],
    )
    def test_display_type(self, raw, expected):
        assert display_type(raw) == expected


class TestPositivity:
    def test_dlp_at_threshold_positive(self):
        assert is_dlp_positive({"dlp_extracellular_prob": "0.8"}, 0.8)

    def test_dlp_legacy_column(self):
        assert is_dlp_positive({"extracellular_prob": "0.9"}, 0.8)

    def test_dlp_malformed_negative(self):
        assert not is_dlp_positive({"dlp_extracellular_prob": "n/a"}, 0.8)

    def test_dse_t3ss_excluded(self):
        assert not is_dse_positive({"dse_ss_type": "T3SS", "dse_max_prob": "0.99"}, 0.8)

    def test_dse_secreted_positive(self):
        assert is_dse_positive({"dse_ss_type": "T2SS", "dse_max_prob": "0.9"}, 0.8)

    def test_self_positive_via_outer_membrane(self):
        # autotransporter: OM-positive even when extracellular is low
        assert is_dlp_self_positive({"outer_membrane_prob": "0.9", "dlp_extracellular_prob": "0.1"}, 0.8)

    def test_self_positive_via_extracellular(self):
        assert is_dlp_self_positive({"outer_membrane_prob": "0.1", "dlp_extracellular_prob": "0.9"}, 0.8)

    def test_self_negative_when_both_low(self):
        assert not is_dlp_self_positive({"outer_membrane_prob": "0.3", "dlp_extracellular_prob": "0.3"}, 0.8)


class TestRotationCounts:
    def test_matches_brute_force(self):
        rng = np.random.default_rng(0)
        for _ in range(200):
            n = int(rng.integers(8, 40))
            pos = (rng.random(n) < 0.3).astype(float)
            mask = (rng.random(n) < 0.3).astype(float)
            c = rotation_counts(pos, mask)
            brute = np.array([int(np.sum(np.roll(pos, k) * mask)) for k in range(n)])
            assert np.array_equal(c, brute)

    def test_offset_zero_is_observed(self):
        pos = np.array([1, 0, 1, 0, 0, 1], dtype=float)
        mask = np.array([1, 1, 0, 0, 1, 0], dtype=float)
        assert int(rotation_counts(pos, mask)[0]) == int(np.sum(pos * mask))


class TestMasks:
    def test_window_excludes_components(self):
        # component at index 5, window 2 -> {3,4,5,6,7}; exclude the component (5)
        m = window_mask([5], n=20, w=2, exclude=[5])
        on = set(np.flatnonzero(m).tolist())
        assert on == {3, 4, 6, 7}

    def test_window_wraps_circularly(self):
        m = window_mask([0], n=10, w=2, exclude=[0])
        assert set(np.flatnonzero(m).tolist()) == {8, 9, 1, 2}

    def test_self_mask_marks_components_only(self):
        m = self_mask([2, 7], n=10)
        assert set(np.flatnonzero(m).tolist()) == {2, 7}


class TestSingleTest:
    def test_planted_enrichment_high_fold(self):
        n = 200
        pos = np.zeros(n)
        comp = [20, 50, 120]  # spaced so flanking positives don't fall on another component
        for p in comp:
            pos[(p + 1) % n] = 1
            pos[(p - 1) % n] = 1
        mask = window_mask(comp, n, 3, exclude=comp)
        r = single_test(pos, mask)
        assert r["observed"] == 6
        assert r["fold"] > 5
        assert r["n_rotations"] == n
        assert r["null"].size == n - 1  # exact all-rotations except identity

    def test_no_positives_fold_zero(self):
        n = 50
        r = single_test(np.zeros(n), self_mask([10], n))
        assert r["observed"] == 0
        assert r["fold"] == 0.0


class TestBhFdr:
    def test_significance_and_monotone(self):
        rows = [{"ss_type": s, "p_perm": p} for s, p in [("a", 0.001), ("b", 0.20), ("c", 0.01), ("d", 0.04)]]
        bh_fdr(rows, pvalue_key="p_perm")
        by = {r["ss_type"]: r for r in rows}
        assert by["a"]["significant"] and by["c"]["significant"]
        assert not by["d"]["significant"] and not by["b"]["significant"]
        qs = [r["qvalue"] for r in sorted(rows, key=lambda r: r["p_perm"])]
        assert qs == sorted(qs)


class TestBinomRetained:
    # binom_pvalue is retained for the validation/comparison scripts, not the test
    def test_obvious_enrichment_low_p(self):
        assert binom_pvalue(8, 10, 0.1) < 1e-4

    def test_degenerate_returns_one(self):
        assert binom_pvalue(0, 0, 0.5) == 1.0


def _build_genome(tmp_dir, n=120):
    """Synthetic whole-genome fixture: one T6SSi system, one T5aSS autotransporter,
    one T3SS system, with planted positives. Returns the four input paths."""
    gene_order = os.path.join(tmp_dir, "gene_order.tsv")
    write_tsv(
        gene_order,
        ["contig", "locus_tag", "position"],
        [{"contig": "c1", "locus_tag": f"g{i:04d}", "position": i} for i in range(n)],
        delimiter="\t",
    )
    ss = os.path.join(tmp_dir, "ss_components.tsv")
    write_tsv(
        ss,
        ["locus_tag", "ss_type", "sys_id", "excluded"],
        [
            {"locus_tag": "g0020", "ss_type": "T6SSi", "sys_id": "sA", "excluded": "False"},
            {"locus_tag": "g0021", "ss_type": "T6SSi", "sys_id": "sA", "excluded": "False"},
            {"locus_tag": "g0060", "ss_type": "T5aSS", "sys_id": "sB", "excluded": "False"},
            {"locus_tag": "g0090", "ss_type": "T3SS", "sys_id": "sC", "excluded": "False"},
            {"locus_tag": "g0091", "ss_type": "T3SS", "sys_id": "sC", "excluded": "False"},
        ],
        delimiter="\t",
    )
    dlp_pos = {"g0019", "g0022"}  # flank the T6SS components
    dse_pos = {"g0019": "T6SS", "g0089": "T3SS"}  # one near T6SS; one is the T3SS comp (DSE excluded)
    dlp_rows, dse_rows = [], []
    for i in range(n):
        lt = f"g{i:04d}"
        dlp_rows.append(
            {
                "locus_tag": lt,
                "dlp_extracellular_prob": "0.95" if lt in dlp_pos else "0.05",
                "outer_membrane_prob": "0.95" if lt == "g0060" else "0.05",  # autotransporter self
            }
        )
        ty = dse_pos.get(lt, "Non-secreted")
        dse_rows.append({"locus_tag": lt, "dse_ss_type": ty, "dse_max_prob": "0.95" if lt in dse_pos else "0.1"})
    dlp = write_tsv(
        os.path.join(tmp_dir, "dlp.tsv"),
        ["locus_tag", "dlp_extracellular_prob", "outer_membrane_prob"],
        dlp_rows,
        delimiter="\t",
    )
    dse = write_tsv(
        os.path.join(tmp_dir, "dse.tsv"), ["locus_tag", "dse_ss_type", "dse_max_prob"], dse_rows, delimiter="\t"
    )
    return {"ss": ss, "gene_order": gene_order, "dlp": dlp, "dse": dse}


class TestRunEnrichmentAggregation:
    def test_per_type_self_and_window(self, tmp_dir):
        fx = _build_genome(tmp_dir)
        dlp = {r["locus_tag"]: r for r in csv.DictReader(open(fx["dlp"]), delimiter="\t")}
        dse = {r["locus_tag"]: r for r in csv.DictReader(open(fx["dse"]), delimiter="\t")}
        systems, ss_type_of_sys = load_systems(fx["ss"])
        # order
        from enrichment_testing import gene_order_flat

        order = gene_order_flat(fx["gene_order"])
        idx, dlp_vec, dse_vec, dlp_self = positivity_vectors(order, dlp, dse, 0.8)
        vecs = {"dlp": dlp_vec, "dse": dse_vec, "dlp_self": dlp_self}
        by_type, all_components = components_by_display_type(systems, ss_type_of_sys)
        all_comp_idx = [idx[lt] for lt in all_components if lt in idx]
        rows = run_enrichment(order, idx, vecs, by_type, all_comp_idx, window=3)

        by = {(r["ss_type"], r["tool"]): r for r in rows}
        # T6SS window: two flanking positives counted (components excluded)
        assert by[("T6SS", "DLP")]["mode"] == "window"
        assert by[("T6SS", "DLP")]["observed"] == 2
        # T5aSS self-detection: the component itself is OM-positive
        assert by[("T5aSS", "DLP")]["mode"] == "self"
        assert by[("T5aSS", "DLP")]["observed"] == 1
        # DSE on T3SS is skipped (DSE can't call T3SS)
        assert by[("T3SS", "DSE")]["skip"] is True
        assert by[("T3SS", "DLP")]["skip"] is False  # DLP T3SS window still tested


class TestCliDriver:
    def test_end_to_end_schema_and_npz(self, monkeypatch, tmp_dir):
        fx = _build_genome(tmp_dir)
        out = os.path.join(tmp_dir, "stats.tsv")
        nulls = os.path.join(tmp_dir, "nulls.npz")
        argv = [
            "enrichment_testing.py",
            "--ss-components", fx["ss"],
            "--gene-order", fx["gene_order"],
            "--dlp", fx["dlp"],
            "--dse", fx["dse"],
            "--window", "3",
            "--conf-threshold", "0.8",
            "--sample", "test",
            "--out", out,
            "--nulls-out", nulls,
        ]  # fmt: skip
        run_script_main(monkeypatch, enrichment_main, argv)
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        cols = set(rows[0].keys())
        assert {"ss_type", "tool", "mode", "observed", "fold", "p_perm", "qvalue", "significant"} <= cols
        assert "M" not in cols and "p_bg" not in cols  # old binomial schema gone
        # T3SS-DSE row present but blank/insignificant
        t3_dse = next(r for r in rows if r["ss_type"] == "T3SS" and r["tool"] == "DSE")
        assert t3_dse["observed"] == "" and t3_dse["significant"] == "False"
        # npz dumped with a key per real (type, tool)
        assert os.path.exists(nulls)
        with np.load(nulls) as npz:
            assert "T6SS__DLP" in npz.files

    def test_no_components_header_only(self, monkeypatch, tmp_dir):
        fx = _build_genome(tmp_dir)
        empty = os.path.join(tmp_dir, "empty_ss.tsv")
        with open(empty, "w") as f:
            f.write("locus_tag\tss_type\tsys_id\texcluded\n")
        out = os.path.join(tmp_dir, "stats.tsv")
        argv = [
            "enrichment_testing.py",
            "--ss-components", empty,
            "--gene-order", fx["gene_order"],
            "--dlp", fx["dlp"],
            "--dse", fx["dse"],
            "--sample", "test",
            "--out", out,
        ]  # fmt: skip
        run_script_main(monkeypatch, enrichment_main, argv)
        assert list(csv.DictReader(open(out), delimiter="\t")) == []
