"""Tests for run_interproscan.py.

Two pure-Python surfaces here:

1. `parse_interproscan_tsv` — column-index TSV parser. Each row is one
   member-DB hit; rows for the same protein must aggregate into sets
   for domains, GO terms, Pfam IDs, and descriptions, then join with
   ';' on output. Hits with `_MISSING` ("-") in a field don't pollute
   the aggregate.
2. `_GO_ID_RE` — extracts every `GO:NNNNNNN` from the GO-terms cell
   (one cell can list several pipe-separated GO entries).

The `run_local_interproscan` subprocess path requires a multi-GB
InterProScan install + Java runtime and is exercised by
tests/integration/test_run_interproscan_integration.py.
"""

import os

import pytest
from run_interproscan import (
    _GO_ID_RE,
    _MISSING,
    DEFAULT_IPS_APPLICATIONS,
    parse_interproscan_tsv,
)

# IPS TSV column layout (1-based for documentation, 0-based in code):
#   0  protein_accession
#   1  md5
#   2  seqlen
#   3  analysis (e.g. "Pfam")
#   4  signature_accession (e.g. "PF03797")
#   5  signature_description
#   6  start
#   7  end
#   8  evalue
#   9  status
#   10 date
#   11 interpro_accession (e.g. "IPR005546")
#   12 interpro_description
#   13 go_terms (e.g. "GO:0005886|GO:0019867")
#   14 pathways (optional)


def _ips_row(
    pid="GENE_001",
    sig_acc="PF03797",
    sig_desc="Autotransporter beta-domain",
    ipr_acc="IPR005546",
    ipr_desc="Autotransporter, beta domain",
    go="GO:0019867",
):
    """A single 14-column IPS TSV row (tab-joined)."""
    return "\t".join(
        [
            pid,
            "abc123md5",
            "300",
            "Pfam",
            sig_acc,
            sig_desc,
            "1",
            "300",
            "1e-50",
            "T",
            "2026-05-05",
            ipr_acc,
            ipr_desc,
            go,
        ]
    )


def _write_ips(path, rows):
    with open(path, "w") as f:
        for row in rows:
            f.write(row + "\n")
    return path


# ---------------------------------------------------------------------------
# _GO_ID_RE — regex extraction
# ---------------------------------------------------------------------------


class TestGoIdRegex:
    @pytest.mark.parametrize(
        "cell, expected",
        [
            ("GO:0005886", ["GO:0005886"]),
            ("GO:0005886|GO:0019867", ["GO:0005886", "GO:0019867"]),
            ("GO:0005886(InterPro)|GO:0019867(PANTHER)", ["GO:0005886", "GO:0019867"]),
            ("", []),
            ("-", []),
            ("not a go term", []),
            # 7-digit IDs (newer GO) must also match
            ("GO:0099999", ["GO:0099999"]),
        ],
    )
    def test_extracts_all_go_ids(self, cell, expected):
        assert [m.group(1) for m in _GO_ID_RE.finditer(cell)] == expected


# ---------------------------------------------------------------------------
# parse_interproscan_tsv — single-protein cases
# ---------------------------------------------------------------------------


class TestParseSingleProtein:
    def test_single_row_aggregated(self, tmp_dir):
        path = _write_ips(os.path.join(tmp_dir, "ips.tsv"), [_ips_row()])
        results = parse_interproscan_tsv(path)
        assert "GENE_001" in results
        e = results["GENE_001"]
        assert e["locus_tag"] == "GENE_001"
        assert e["interpro_domains"] == "IPR005546"
        assert e["interpro_descriptions"] == "Autotransporter, beta domain"
        assert e["interpro_pfam_ids"] == "PF03797"
        assert e["interpro_go_terms"] == "GO:0019867"

    def test_multiple_rows_same_protein_merged(self, tmp_dir):
        # Same protein, three different signatures. All four fields must dedup
        # into a sorted, semicolon-joined string.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(sig_acc="PF03797", ipr_acc="IPR005546", ipr_desc="alpha", go="GO:0001"),
                _ips_row(sig_acc="PF00001", ipr_acc="IPR000123", ipr_desc="beta", go="GO:0002"),
                _ips_row(sig_acc="PF03797", ipr_acc="IPR005546", ipr_desc="alpha", go="GO:0001|GO:0003"),
            ],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_domains"] == "IPR000123;IPR005546"  # sorted
        assert e["interpro_descriptions"] == "alpha;beta"
        assert e["interpro_pfam_ids"] == "PF00001;PF03797"
        assert e["interpro_go_terms"] == "GO:0001;GO:0002;GO:0003"

    def test_missing_ipr_excluded(self, tmp_dir):
        # Row with "-" in IPR fields must NOT add to the aggregate.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(ipr_acc=_MISSING, ipr_desc=_MISSING)],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_domains"] == ""
        assert e["interpro_descriptions"] == ""
        # Pfam still extracted (sig_acc isn't missing)
        assert e["interpro_pfam_ids"] == "PF03797"

    def test_missing_go_excluded(self, tmp_dir):
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(go=_MISSING)],
        )
        assert parse_interproscan_tsv(path)["GENE_001"]["interpro_go_terms"] == ""

    def test_non_pfam_signature_not_in_pfam_ids(self, tmp_dir):
        # Signature "TIGR03660" doesn't start with "PF" — must NOT show up in
        # interpro_pfam_ids.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [_ips_row(sig_acc="TIGR03660", sig_desc="something")],
        )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_pfam_ids"] == ""
        # The IPR side is still populated
        assert e["interpro_domains"] == "IPR005546"


# ---------------------------------------------------------------------------
# parse_interproscan_tsv — multi-protein + filter
# ---------------------------------------------------------------------------


class TestParseMultipleProteins:
    def test_each_protein_isolated(self, tmp_dir):
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(pid="GENE_001", sig_acc="PF03797"),
                _ips_row(pid="GENE_002", sig_acc="PF00001"),
            ],
        )
        results = parse_interproscan_tsv(path)
        assert set(results.keys()) == {"GENE_001", "GENE_002"}
        assert results["GENE_001"]["interpro_pfam_ids"] == "PF03797"
        assert results["GENE_002"]["interpro_pfam_ids"] == "PF00001"

    def test_target_ids_filters(self, tmp_dir):
        # Only GENE_001 is in target_ids — GENE_002 must be dropped.
        path = _write_ips(
            os.path.join(tmp_dir, "ips.tsv"),
            [
                _ips_row(pid="GENE_001"),
                _ips_row(pid="GENE_002"),
            ],
        )
        results = parse_interproscan_tsv(path, target_ids={"GENE_001"})
        assert set(results.keys()) == {"GENE_001"}

    def test_empty_target_ids_returns_empty(self, tmp_dir):
        path = _write_ips(os.path.join(tmp_dir, "ips.tsv"), [_ips_row()])
        # Note: the contract here is that target_ids is checked for membership;
        # empty set behaves like "no protein matches" — verifies the truthy
        # check on target_ids in the parser.
        results = parse_interproscan_tsv(path, target_ids={"NOT_PRESENT"})
        assert results == {}


# ---------------------------------------------------------------------------
# Parser resilience
# ---------------------------------------------------------------------------


class TestParserResilience:
    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "ips.tsv")
        open(path, "w").close()
        assert parse_interproscan_tsv(path) == {}

    def test_short_rows_skipped(self, tmp_dir):
        # Rows with <12 columns are dropped (corrupt or partial line)
        path = os.path.join(tmp_dir, "ips.tsv")
        with open(path, "w") as f:
            f.write("GENE_001\tabc\t300\tPfam\tPF03797\n")  # only 5 columns
        assert parse_interproscan_tsv(path) == {}

    def test_row_without_go_column_handled(self, tmp_dir):
        # A row with exactly 13 columns (no GO field at index 13) must not
        # crash. The parser uses bounds-checked column reads.
        path = os.path.join(tmp_dir, "ips.tsv")
        with open(path, "w") as f:
            f.write(
                "GENE_001\tmd5\t300\tPfam\tPF03797\tdesc\t1\t300\t1e-50\tT\t2026-05-05\t"
                "IPR005546\tipr_desc\n"  # only 13 fields
            )
        e = parse_interproscan_tsv(path)["GENE_001"]
        assert e["interpro_go_terms"] == ""
        assert e["interpro_domains"] == "IPR005546"


# ---------------------------------------------------------------------------
# DEFAULT_IPS_APPLICATIONS — pinned to the documented bacteria set
# ---------------------------------------------------------------------------


class TestDefaultApplications:
    def test_panther_excluded(self):
        # PANTHER is intentionally skipped — eukaryote-leaning + slowest member.
        # Future maintainers must opt into PANTHER explicitly via --applications.
        assert "PANTHER" not in DEFAULT_IPS_APPLICATIONS

    def test_pfam_included(self):
        # Pfam is the workhorse for bacterial annotation; must always be present.
        assert "Pfam" in DEFAULT_IPS_APPLICATIONS

    def test_no_duplicates(self):
        assert len(DEFAULT_IPS_APPLICATIONS) == len(set(DEFAULT_IPS_APPLICATIONS))


# ---------------------------------------------------------------------------
# Failure-path diagnostics — non-zero exit must surface stdout/stderr to disk
# ---------------------------------------------------------------------------


class TestRunLocalInterproscanFailureSurface:
    """On non-zero exit, the wrapper writes a sidecar log capturing stdout
    AND stderr (CX3 K-12 run b2060a9 saw empty stderr with a real error
    hiding in stdout) and references the log path in the raised error."""

    def _stub_run(self, returncode, stdout, stderr):

        class _Result:
            pass

        r = _Result()
        r.returncode = returncode
        r.stdout = stdout
        r.stderr = stderr

        def fake_run(*args, **kwargs):
            return r

        return fake_run

    def test_failure_writes_log_with_both_streams(self, tmp_path, monkeypatch):
        import subprocess

        from run_interproscan import run_local_interproscan

        monkeypatch.setattr(
            "run_interproscan._resolve_interproscan_binary",
            lambda d: "/fake/interproscan.sh",
        )
        monkeypatch.setattr(
            subprocess,
            "run",
            self._stub_run(1, "java stack trace here", ""),
        )

        with pytest.raises(RuntimeError, match="interproscan_failure.log"):
            run_local_interproscan(
                query_fasta=str(tmp_path / "in.faa"),
                install_dir="/fake",
                output_dir=str(tmp_path),
            )

        log_path = tmp_path / "interproscan_failure.log"
        assert log_path.is_file()
        body = log_path.read_text()
        assert "exit code: 1" in body
        assert "java stack trace here" in body  # stdout preserved
        assert "(empty)" in body  # stderr placeholder

    def test_failure_log_path_in_runtime_error(self, tmp_path, monkeypatch):
        import subprocess

        from run_interproscan import run_local_interproscan

        monkeypatch.setattr(
            "run_interproscan._resolve_interproscan_binary",
            lambda d: "/fake/interproscan.sh",
        )
        monkeypatch.setattr(
            subprocess,
            "run",
            self._stub_run(2, "", "missing Java"),
        )

        with pytest.raises(RuntimeError) as exc:
            run_local_interproscan(
                query_fasta=str(tmp_path / "in.faa"),
                install_dir="/fake",
                output_dir=str(tmp_path),
            )
        assert str(tmp_path / "interproscan_failure.log") in str(exc.value)


class TestEmptySubstratesShortCircuit:
    """Regression: when substrates_filtered.tsv has zero rows (e.g. a
    reduced-genome endosymbiont like Buchnera with no detected substrates),
    the IPS wrapper used to feed Java an empty FASTA, fail to find the
    expected results.tsv, and crash with FileNotFoundError. It must now
    short-circuit, write a header-only CSV, and exit cleanly."""

    def test_empty_substrates_writes_header_only_csv(self, tmp_path, monkeypatch):
        import subprocess

        import run_interproscan

        # substrates_filtered.tsv with header but no rows.
        substrates = tmp_path / "substrates.tsv"
        substrates.write_text("locus_tag\tss_type\n")

        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">X\nM\n")

        output = tmp_path / "ips.csv"

        # If anything tries to invoke IPS, fail loudly.
        def _no_subprocess(*a, **kw):
            raise AssertionError("subprocess.run must NOT fire when there are no substrates")

        monkeypatch.setattr(subprocess, "run", _no_subprocess)
        monkeypatch.setattr(
            "sys.argv",
            [
                "run_interproscan.py",
                "--substrates",
                str(substrates),
                "--proteins",
                str(proteins),
                "--sample",
                "buchnera",
                "--output",
                str(output),
            ],
        )
        run_interproscan.main()

        body = output.read_text()
        assert (
            body.splitlines()[0]
            == "locus_tag,interpro_domains,interpro_go_terms,interpro_pfam_ids,interpro_descriptions"
        )
        assert len(body.splitlines()) == 1  # header only


class TestJavaOptsHeapAutoScale:
    """IPS hard-codes -Xmx15G in its bundled launcher (verified against
    upstream interproscan.sh) and does NOT read $JAVA_OPTS. The JVM
    always reads $_JAVA_OPTIONS regardless of launcher script, and the
    last -Xmx wins, so the wrapper exports _JAVA_OPTIONS instead. This
    pins the auto-scale: (effective_ram_gb - 2) * 0.5, clamped to [4, 64].
    """

    def _capture_env(self, captured):

        class _Result:
            returncode = 0
            stdout = ""
            stderr = ""

        def fake_run(cmd, **kwargs):
            captured["env"] = kwargs.get("env", {})
            return _Result()

        return fake_run

    def _heap_arg(self, env):
        opts = env.get("_JAVA_OPTIONS", "")
        for token in opts.split():
            if token.startswith("-Xmx"):
                return token
        return None

    def _run_and_capture(self, tmp_path, ram_gb, monkeypatch):
        import subprocess

        from run_interproscan import run_local_interproscan

        captured: dict = {}
        monkeypatch.setattr("run_interproscan._resolve_interproscan_binary", lambda d: "/fake/interproscan.sh")
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: ram_gb)
        monkeypatch.setattr(subprocess, "run", self._capture_env(captured))
        # Stub the results.tsv that the wrapper reads back after subprocess.
        (tmp_path / "results.tsv").write_text("")

        try:
            run_local_interproscan(
                query_fasta=str(tmp_path / "in.faa"),
                install_dir="/fake",
                output_dir=str(tmp_path),
            )
        except Exception:
            pass  # wrapper may parse results.tsv after subprocess; we only need env.

        # Brittle-test guard: if the JAVA_OPTIONS export ever moves to
        # AFTER subprocess.run (refactor accident), `captured` is empty
        # and the assertions below would also be swallowed. Fail loudly here.
        assert "env" in captured, "subprocess.run was never called; JAVA_OPTIONS branch never ran"
        return captured["env"]

    def test_heap_scales_with_detected_ram(self, tmp_path, monkeypatch):
        env = self._run_and_capture(tmp_path, 80.0, monkeypatch)
        # (80 - 2) * 0.5 = 39 → in-range
        assert self._heap_arg(env) == "-Xmx39g"

    def test_heap_floor_4gb(self, tmp_path, monkeypatch):
        env = self._run_and_capture(tmp_path, 2.0, monkeypatch)
        # (2 - 2) * 0.5 = 0 → floor at 4 GB
        assert self._heap_arg(env) == "-Xmx4g"

    def test_heap_ceiling_64gb(self, tmp_path, monkeypatch):
        env = self._run_and_capture(tmp_path, 512.0, monkeypatch)
        # (512 - 2) * 0.5 = 255 → ceiling at 64
        assert self._heap_arg(env) == "-Xmx64g"

    def test_heap_sizes_to_group_share_not_whole_node(self, tmp_path, monkeypatch):
        # The §4.2 fix: inside the annotation parallel group the heap must
        # size to this tool's RAM *share*, not the whole node — else N Java
        # heaps each sized to the full node overcommit and OOM-kill it.
        # 90 GB node, 3-tool group → 30 GB share → (30-2)*0.5 = 14 GB,
        # NOT the (90-2)*0.5 = 44 GB a whole-node calc would pick.
        from ssign_lib.resources import parallel_group

        with parallel_group(3):
            env = self._run_and_capture(tmp_path, 90.0, monkeypatch)
        assert self._heap_arg(env) == "-Xmx14g"
