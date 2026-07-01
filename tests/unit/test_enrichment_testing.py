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
    is_signalp_positive,
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

    @pytest.mark.parametrize(
        "pred,expected",
        [
            ("SP", True),  # Sec/SPI
            ("LIPO", True),  # Sec/SPII
            ("TAT", False),  # Tat/SPI — different export route
            ("TATLIPO", False),
            ("PILIN", False),  # pilus marker, not a secretion signal
            ("OTHER", False),
            ("", False),
        ],
    )
    def test_signalp_sec_class_rule(self, pred, expected):
        assert is_signalp_positive({"signalp_prediction": pred}) is expected

    def test_signalp_missing_row_negative(self):
        assert not is_signalp_positive({})


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
        assert r["n_null"] == n - 1  # null size = the n-1 non-identity rotations
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
    """Synthetic whole-genome fixture: one T6SSi system, one T5aSS autotransporter
    (DLP-self + SignalP positive), one T5cSS autotransporter (SignalP positive but
    DLP-self negative, to pin the combined-uses-SignalP rule), one T3SS system, with
    planted positives. Returns the five input paths (incl. signalp)."""
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
            {"locus_tag": "g0030", "ss_type": "T5cSS", "sys_id": "sD", "excluded": "False"},
            {"locus_tag": "g0060", "ss_type": "T5aSS", "sys_id": "sB", "excluded": "False"},
            {"locus_tag": "g0090", "ss_type": "T3SS", "sys_id": "sC", "excluded": "False"},
            {"locus_tag": "g0091", "ss_type": "T3SS", "sys_id": "sC", "excluded": "False"},
        ],
        delimiter="\t",
    )
    dlp_pos = {"g0019", "g0022"}  # flank the T6SS components
    dse_pos = {"g0019": "T6SS", "g0089": "T3SS"}  # one near T6SS; one is the T3SS comp (DSE excluded)
    # SignalP Sec-signal calls. The two autotransporters deliberately DISAGREE with
    # DLP so the combined DLP-or-SignalP union is testable: T5aSS (g0060) is
    # DLP-self-positive but SignalP-negative (combined must still fire via DLP);
    # T5cSS (g0030) is SignalP-positive but DLP-self-negative (combined fires via
    # SignalP). g0005/g0110 are scattered background Sec signals; g0009 is a TAT
    # (non-Sec) negative control.
    sp_pred = {"g0030": "SP", "g0005": "SP", "g0110": "LIPO", "g0009": "TAT"}
    dlp_rows, dse_rows, sp_rows = [], [], []
    for i in range(n):
        lt = f"g{i:04d}"
        dlp_rows.append(
            {
                "locus_tag": lt,
                "dlp_extracellular_prob": "0.95" if lt in dlp_pos else "0.05",
                "outer_membrane_prob": "0.95" if lt == "g0060" else "0.05",  # only T5aSS is DLP-self positive
            }
        )
        ty = dse_pos.get(lt, "Non-secreted")
        dse_rows.append({"locus_tag": lt, "dse_ss_type": ty, "dse_max_prob": "0.95" if lt in dse_pos else "0.1"})
        sp_rows.append({"locus_tag": lt, "signalp_prediction": sp_pred.get(lt, "OTHER")})
    dlp = write_tsv(
        os.path.join(tmp_dir, "dlp.tsv"),
        ["locus_tag", "dlp_extracellular_prob", "outer_membrane_prob"],
        dlp_rows,
        delimiter="\t",
    )
    dse = write_tsv(
        os.path.join(tmp_dir, "dse.tsv"), ["locus_tag", "dse_ss_type", "dse_max_prob"], dse_rows, delimiter="\t"
    )
    signalp = write_tsv(
        os.path.join(tmp_dir, "signalp.tsv"), ["locus_tag", "signalp_prediction"], sp_rows, delimiter="\t"
    )
    return {"ss": ss, "gene_order": gene_order, "dlp": dlp, "dse": dse, "signalp": signalp}


class TestRunEnrichmentAggregation:
    def test_per_type_self_and_window(self, tmp_dir):
        fx = _build_genome(tmp_dir)
        dlp = {r["locus_tag"]: r for r in csv.DictReader(open(fx["dlp"]), delimiter="\t")}
        dse = {r["locus_tag"]: r for r in csv.DictReader(open(fx["dse"]), delimiter="\t")}
        signalp = {r["locus_tag"]: r for r in csv.DictReader(open(fx["signalp"]), delimiter="\t")}
        systems, ss_type_of_sys = load_systems(fx["ss"])
        # order
        from enrichment_testing import gene_order_flat

        order = gene_order_flat(fx["gene_order"])
        idx, vecs = positivity_vectors(order, dlp, dse, signalp, 0.8)
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

        # SignalP is a full third predictor: emitted for every type incl. T3SS (no
        # skip), self-detection for autotransporters. In this fixture T5aSS is
        # SignalP-negative (only T5cSS carries the signal), so its self observed is 0.
        assert by[("T5aSS", "SignalP")]["mode"] == "self"
        assert by[("T5aSS", "SignalP")]["observed"] == 0
        assert by[("T5cSS", "SignalP")]["observed"] == 1  # autotransporter carries a Sec signal
        assert by[("T3SS", "SignalP")]["skip"] is False  # unlike DSE, SignalP IS tested for T3SS
        assert by[("T6SS", "SignalP")]["mode"] == "window"

        # Combined non-T5 track = DLP-or-DSE: one row per tested type; observed is the
        # OR of the valid predictors, so >= each individual predictor.
        assert by[("T6SS", "COMBINED")]["mode"] == "window"
        assert by[("T6SS", "COMBINED")]["observed"] >= by[("T6SS", "DLP")]["observed"]
        assert by[("T6SS", "COMBINED")]["observed"] >= by[("T6SS", "DSE")]["observed"]
        # T3SS combined uses DLP only (DSE invalid for T3SS).
        assert by[("T3SS", "COMBINED")]["skip"] is False
        assert by[("T3SS", "COMBINED")]["observed"] == by[("T3SS", "DLP")]["observed"]
        # All T5SS combined = DLP-or-SignalP union (NOT DLP-or-DSE, NOT SignalP-alone).
        # The two autotransporters disagree, so the union is provable both ways:
        #  - T5aSS: DLP-self=1, SignalP=0 -> combined follows DLP (1), so SignalP alone
        #    (0) would be wrong.
        assert by[("T5aSS", "DLP")]["observed"] == 1
        assert by[("T5aSS", "COMBINED")]["observed"] == by[("T5aSS", "DLP")]["observed"]
        #  - T5cSS: DLP-self=0, SignalP=1 -> combined follows SignalP (1), so DLP-or-DSE
        #    (0) would be wrong (this is the SignalP rescue).
        assert by[("T5cSS", "DLP")]["observed"] == 0
        assert by[("T5cSS", "COMBINED")]["observed"] == by[("T5cSS", "SignalP")]["observed"] == 1


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
            "--signalp", fx["signalp"],
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
        # SignalP rows emitted (third predictor); it joins the per-tool BH family so its
        # rows carry a qvalue, and it is present even for T3SS (not skipped).
        sp_rows = [r for r in rows if r["tool"] == "SignalP"]
        assert sp_rows and all(r["qvalue"] != "" for r in sp_rows)
        assert any(r["ss_type"] == "T3SS" for r in sp_rows)
        # Combined DLP-or-DSE rows emitted (one per tested type) for the combined figure
        assert any(r["tool"] == "COMBINED" for r in rows)
        # npz dumped with a key per real (type, tool), SignalP included
        assert os.path.exists(nulls)
        with np.load(nulls) as npz:
            assert "T6SS__DLP" in npz.files
            assert "T5aSS__SignalP" in npz.files

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
            "--signalp", fx["signalp"],
            "--sample", "test",
            "--out", out,
        ]  # fmt: skip
        run_script_main(monkeypatch, enrichment_main, argv)
        assert list(csv.DictReader(open(out), delimiter="\t")) == []
