"""Unit tests for run_eggnog.py — parser for emapper .annotations output.

Does not exercise the `emapper.py` subprocess itself (that requires a
~50 GB EggNOG database and a working eggnog-mapper install). Tests target
the pure-Python parser that reads the tool's TSV output.
"""

import os

import pytest
from run_eggnog import (
    _autodetect_block_size,
    _autodetect_dbmem,
    _autodetect_index_chunks,
    _build_emapper_cmd,
    _split_rich_field,
    _stage_eggnog_db_to_local,
    load_substrate_ids,
    parse_eggnog_annotations,
    write_substrates_only_fasta,
)


class TestBuildEmapperCmd:
    """Verify the emapper.py argv assembly.

    `--dbmem` loads ~44 GB resident; gated on host RAM by default so 32 GB
    nodes (e.g. Imperial CX3 GPU partition) don't OOM-kill emapper mid-load
    with no clean error. Explicit True/False overrides the auto-detect.
    """

    def _common_args(self):
        return dict(
            proteins_fasta="/tmp/proteins.faa",
            db_path="/db/eggnog",
            sample_id="sample",
            output_dir="/tmp/out",
        )

    def test_dbmem_forced_on(self):
        cmd = _build_emapper_cmd(**self._common_args(), dbmem=True)
        assert cmd.count("--dbmem") == 1

    def test_dbmem_forced_off(self):
        cmd = _build_emapper_cmd(**self._common_args(), dbmem=False)
        assert "--dbmem" not in cmd

    def test_dbmem_autodetect_low_ram(self, monkeypatch):
        monkeypatch.setattr("run_eggnog._autodetect_dbmem", lambda: False)
        cmd = _build_emapper_cmd(**self._common_args())  # dbmem=None
        assert "--dbmem" not in cmd

    def test_dbmem_autodetect_high_ram(self, monkeypatch):
        monkeypatch.setattr("run_eggnog._autodetect_dbmem", lambda: True)
        cmd = _build_emapper_cmd(**self._common_args())
        assert cmd.count("--dbmem") == 1

    def test_sensmode_default_is_sensitive(self):
        cmd = _build_emapper_cmd(**self._common_args())
        i = cmd.index("--sensmode")
        assert cmd[i + 1] == "sensitive"

    def test_required_flags_present(self):
        """Sanity: the existing fixed flags survive the refactor."""
        cmd = _build_emapper_cmd(**self._common_args())
        for flag in (
            "-i",
            "--itype",
            "--output_dir",
            "-o",
            "--data_dir",
            "--cpu",
            "--tax_scope",
            "--sensmode",
            "--override",
            "--block_size",
            "--index_chunks",
        ):
            assert flag in cmd, f"expected {flag} in cmd, got {cmd}"

    def test_block_size_forced(self):
        cmd = _build_emapper_cmd(**self._common_args(), block_size=6)
        i = cmd.index("--block_size")
        assert cmd[i + 1] == "6"

    def test_index_chunks_forced(self):
        cmd = _build_emapper_cmd(**self._common_args(), index_chunks=1)
        i = cmd.index("--index_chunks")
        assert cmd[i + 1] == "1"

    def test_autodetect_emits_default_when_low_ram(self, monkeypatch):
        """Below the index_chunks threshold and the block_size headroom,
        autodetect should emit DIAMOND's own defaults (2 and 4)."""
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 16)
        cmd = _build_emapper_cmd(**self._common_args(), dbmem=False)
        assert cmd[cmd.index("--block_size") + 1] == "2"
        assert cmd[cmd.index("--index_chunks") + 1] == "4"


class TestAutodetectBlockSize:
    """Gating logic for the DIAMOND --block_size passthrough."""

    def test_low_ram_keeps_default(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 16)
        # 16 GB - 16 reserved = 0 GB free → keep default of 2
        assert _autodetect_block_size(dbmem=False) == 2

    def test_dbmem_swallows_most_ram(self, monkeypatch):
        # 50 GB - 44 reserved (dbmem) = 6 GB free → keep default of 2
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 50)
        assert _autodetect_block_size(dbmem=True) == 2

    def test_high_ram_no_dbmem_bumps(self, monkeypatch):
        # 64 GB - 16 reserved = 48 GB free → 48/6 = 8, capped at 8.
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 64)
        assert _autodetect_block_size(dbmem=False) == 8

    def test_capped_at_8(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 256)
        assert _autodetect_block_size(dbmem=False) == 8

    def test_dbmem_on_big_node(self, monkeypatch):
        # 80 GB - 44 reserved = 36 GB free → 36/6 = 6.
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 80)
        assert _autodetect_block_size(dbmem=True) == 6


class TestAutodetectIndexChunks:
    def test_low_ram_keeps_default_4(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 32)
        assert _autodetect_index_chunks() == 4

    def test_high_ram_drops_to_1(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 64)
        assert _autodetect_index_chunks() == 1

    def test_threshold_inclusive(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 60)
        assert _autodetect_index_chunks() == 1


class TestAutodetectParallelGroupAware:
    """EggNOG sizes memory off its RAM *share*, not the whole node. In the
    full-tier 6-tool annotation group a node-sized --dbmem/--block_size
    overcommits and OOM-kills DIAMOND; dividing by the group size fixes it.
    Same 120 GB node, same effective_ram_gb — only the group env differs."""

    def test_dbmem_off_when_group_shrinks_share_below_threshold(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 120)
        # Alone on the node: 120 GB >= 50 → --dbmem on.
        monkeypatch.delenv("SSIGN_PARALLEL_GROUP_SIZE", raising=False)
        assert _autodetect_dbmem() is True
        # In the 6-tool group: 120/6 = 20 GB < 50 → --dbmem off (no overcommit).
        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "6")
        assert _autodetect_dbmem() is False

    def test_block_size_and_chunks_shrink_in_group(self, monkeypatch):
        monkeypatch.setattr("ssign_lib.resources.effective_ram_gb", lambda: 120)
        monkeypatch.setenv("SSIGN_PARALLEL_GROUP_SIZE", "6")
        # share = 20 GB: no --dbmem, block_size falls to the default 2
        # (20-16=4 GB free < 2*6), index_chunks back to 4 (20 < 60).
        assert _autodetect_block_size(dbmem=False) == 2
        assert _autodetect_index_chunks() == 4


# Representative emapper 2.1.x .annotations fixture. Lines starting with
# '##' are comments emitted by emapper before and after the data block;
# the parser must skip them. The column header line begins with '#query'.
_EGGNOG_FIXTURE = (
    "## emapper-2.1.12\n"
    "## /path/to/emapper.py -i proteins.faa -o sample --data_dir /db\n"
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
    "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\t"
    "KEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\t"
    "KEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
    "GENE_00001\t83333.eco0001\t1e-180\t500.0\tCOG0012@1|root\t2|Bacteria\t"
    "M\tDNA-directed RNA polymerase subunit beta\trpoB\tGO:0003677\t"
    "2.7.7.6\tko:K03040\t-\t-\t-\t-\t-\t-\t-\t-\tRNA_pol_Rpb1\n"
    "GENE_00003\t83333.eco0003\t2e-50\t200.0\tCOG0084@1|root\t2|Bacteria\t"
    "L\tUncharacterized conserved protein\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"
    "-\t-\t-\t-\n"
    "## 2 queries annotated\n"
)


class TestStageEggnogDb:
    """Verify the local-DB staging fix for shared-FS hang on emapper mmap."""

    def _make_db(self, src, files=("eggnog.db", "eggnog_proteins.dmnd"), nbytes=128):
        for name in files:
            with open(os.path.join(src, name), "wb") as f:
                f.write(b"\0" * nbytes)

    def test_copies_required_files(self, tmp_dir):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        os.makedirs(src)
        os.makedirs(cache)
        self._make_db(src)
        out = _stage_eggnog_db_to_local(src, cache)
        assert out == os.path.join(cache, "eggnog")
        assert os.path.getsize(os.path.join(out, "eggnog.db")) == 128
        assert os.path.getsize(os.path.join(out, "eggnog_proteins.dmnd")) == 128

    def test_idempotent_on_size_match(self, tmp_dir):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        os.makedirs(src)
        os.makedirs(cache)
        self._make_db(src)
        _stage_eggnog_db_to_local(src, cache)
        mtime = os.path.getmtime(os.path.join(cache, "eggnog", "eggnog.db"))
        _stage_eggnog_db_to_local(src, cache)  # second call must not re-copy
        assert os.path.getmtime(os.path.join(cache, "eggnog", "eggnog.db")) == mtime

    def test_raises_if_required_file_missing(self, tmp_dir):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        os.makedirs(src)
        os.makedirs(cache)
        self._make_db(src, files=("eggnog.db",))  # missing eggnog_proteins.dmnd
        with pytest.raises(FileNotFoundError, match="eggnog_proteins.dmnd"):
            _stage_eggnog_db_to_local(src, cache)

    def test_skips_missing_optional_files(self, tmp_dir):
        src, cache = os.path.join(tmp_dir, "src"), os.path.join(tmp_dir, "cache")
        os.makedirs(src)
        os.makedirs(cache)
        self._make_db(src)  # only the two required files
        out = _stage_eggnog_db_to_local(src, cache)
        # optionals shouldn't be there, but the call must succeed
        assert not os.path.exists(os.path.join(out, "eggnog.taxa.db"))


class TestParseEggnogAnnotations:
    def test_skips_comment_lines(self, tmp_dir):
        """Lines starting with '##' before and after the data block are skipped."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        assert len(entries) == 2  # only the two data rows, no ## comments

    def test_extracts_expected_fields(self, tmp_dir):
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e1 = next(e for e in entries if e["locus_tag"] == "GENE_00001")
        assert e1["seed_ortholog"] == "83333.eco0001"
        assert e1["evalue"] == "1e-180"
        assert e1["eggnog_description"] == "DNA-directed RNA polymerase subunit beta"
        assert e1["preferred_name"] == "rpoB"

    def test_dash_description_preserved_as_dash(self, tmp_dir):
        """A literal '-' in the Description column is a real emapper value
        (no description found) and should pass through unchanged."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e3 = next(e for e in entries if e["locus_tag"] == "GENE_00003")
        # GENE_00003 has a real description in the fixture
        assert e3["eggnog_description"] == "Uncharacterized conserved protein"
        # But its Preferred_name is "-" — should be kept as empty string or "-"
        assert e3["preferred_name"] == "-"

    def test_empty_description_defaults_to_dash(self, tmp_dir):
        """An empty Description field falls back to '-' for consistency."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "GENE_X\t83333.x\t1e-10\t50\tCOG@1\t2|B\tS\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n"
            )

        entries = parse_eggnog_annotations(path)
        assert entries[0]["eggnog_description"] == "-"

    def test_empty_annotations_returns_empty_list(self, tmp_dir):
        """A valid emapper output with zero annotated queries returns []."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "## emapper-2.1.12\n"
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "## 0 queries annotated\n"
            )

        entries = parse_eggnog_annotations(path)
        assert entries == []

    def test_output_shape_is_stable(self, tmp_dir):
        """Every entry must contain the expected fields for downstream consumers."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        required = {
            "locus_tag",
            "seed_ortholog",
            "evalue",
            "eggnog_description",
            "preferred_name",
            "cog_category",
            "ec_numbers",
            "kegg_ko",
            "go_terms",
            "pfam_ids",
        }
        for e in entries:
            assert required <= set(e.keys())


class TestRichFieldsSurfaced:
    def test_annotated_row_populates_rich_fields(self, tmp_dir):
        """GENE_00001 in the fixture has EC, KEGG_ko, GOs, PFAMs, COG_category."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e1 = next(e for e in entries if e["locus_tag"] == "GENE_00001")
        assert e1["cog_category"] == "M"
        assert e1["ec_numbers"] == ["2.7.7.6"]
        assert e1["kegg_ko"] == ["ko:K03040"]
        assert e1["go_terms"] == ["GO:0003677"]
        assert e1["pfam_ids"] == ["RNA_pol_Rpb1"]

    def test_dash_sentinels_become_empty(self, tmp_dir):
        """GENE_00003 has COG_category but '-' for every multi-value field."""
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(_EGGNOG_FIXTURE)

        entries = parse_eggnog_annotations(path)
        e3 = next(e for e in entries if e["locus_tag"] == "GENE_00003")
        assert e3["cog_category"] == "L"
        assert e3["ec_numbers"] == []
        assert e3["kegg_ko"] == []
        assert e3["go_terms"] == []
        assert e3["pfam_ids"] == []

    def test_multi_value_row_splits_on_comma(self, tmp_dir):
        path = os.path.join(tmp_dir, "sample.emapper.annotations")
        with open(path, "w") as f:
            f.write(
                "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\t"
                "max_annot_lvl\tCOG_category\tDescription\tPreferred_name\t"
                "GOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\t"
                "KEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n"
                "GENE_X\t83333.x\t1e-10\t50\tCOG@1\t2|B\tM\tmulti\tname\t"
                "GO:0001,GO:0002\t1.1.1.1,2.2.2.2\tko:K1,ko:K2\t-\t-\t-\t"
                "-\t-\t-\t-\t-\tPF1234,PF5678\n"
            )

        entries = parse_eggnog_annotations(path)
        e = entries[0]
        assert e["go_terms"] == ["GO:0001", "GO:0002"]
        assert e["ec_numbers"] == ["1.1.1.1", "2.2.2.2"]
        assert e["kegg_ko"] == ["ko:K1", "ko:K2"]
        assert e["pfam_ids"] == ["PF1234", "PF5678"]


class TestSplitRichField:
    def test_none_returns_empty_list(self):
        assert _split_rich_field(None) == []

    def test_empty_string_returns_empty_list(self):
        assert _split_rich_field("") == []
        assert _split_rich_field("   ") == []

    def test_dash_sentinel_returns_empty_list(self):
        assert _split_rich_field("-") == []

    def test_single_value_returned_as_one_element_list(self):
        assert _split_rich_field("1.1.1.1") == ["1.1.1.1"]

    def test_comma_separated_values_split(self):
        assert _split_rich_field("a,b,c") == ["a", "b", "c"]

    def test_whitespace_around_entries_stripped(self):
        assert _split_rich_field(" a , b , c ") == ["a", "b", "c"]

    def test_empty_fragments_dropped(self):
        assert _split_rich_field("a,,b") == ["a", "b"]


class TestLoadSubstrateIds:
    """Substrate-TSV locus_tag extraction for filtering the protein FASTA."""

    def test_reads_locus_tags_from_tsv(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tother_col\n")
            f.write("GENE_00001\tx\n")
            f.write("GENE_00002\ty\n")
            f.write("GENE_00003\tz\n")
        assert load_substrate_ids(path) == {"GENE_00001", "GENE_00002", "GENE_00003"}

    def test_empty_locus_tag_rows_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tproduct\n")
            f.write("G1\thit\n")
            f.write("\tghost\n")  # blank locus_tag — skip
            f.write("G2\thit\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_empty_file_returns_empty_set(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\n")  # header only
        assert load_substrate_ids(path) == set()


class TestWriteSubstratesOnlyFasta:
    def _write_proteins(self, path, seqs):
        with open(path, "w") as f:
            for k, v in seqs.items():
                f.write(f">{k}\n{v}\n")

    def test_writes_only_requested_ids(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(
            src,
            {"G1": "MKTLL", "G2": "MFVFL", "G3": "MQKAL", "G4": "MSRLK"},
        )
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G2", "G4"}, out)
        assert n == 2
        with open(out) as f:
            body = f.read()
        assert ">G2" in body
        assert ">G4" in body
        assert ">G1" not in body
        assert ">G3" not in body

    def test_missing_ids_are_silently_dropped(self, tmp_dir):
        """If the substrates TSV references a locus_tag that doesn't exist
        in the protein FASTA (possible if CDS sets drifted), skip rather
        than crash — the output just has fewer rows."""
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(src, {"G1": "MKT"})
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "MISSING"}, out)
        assert n == 1

    def test_empty_substrate_set_writes_empty_fasta(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        self._write_proteins(src, {"G1": "MKT"})
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, set(), out)
        assert n == 0
        assert os.path.exists(out)
        assert os.path.getsize(out) == 0


class TestRunEmapperThreadDefault:
    """`run_emapper(..., threads=None)` must auto-detect from
    `effective_cpu_count()`. Argparse default already did this; the
    Python-API signature pinned 4 and starved direct callers on
    multi-core machines."""

    def _stub_subprocess_run(self, captured):
        class _Done:
            returncode = 0
            stdout = ""
            stderr = ""

        def _fake(cmd, **kwargs):
            captured.append(list(cmd))
            return _Done()

        return _fake

    def _touch_annotations(self, tmp_path, sample_id):
        # run_emapper reads back the .emapper.annotations file after
        # subprocess returns; stub it so the wrapper doesn't FileNotFoundError.
        open(os.path.join(tmp_path, f"{sample_id}.emapper.annotations"), "w").close()

    def test_threads_none_resolves_to_effective_cpu_count(self, monkeypatch, tmp_path):
        import subprocess

        from run_eggnog import run_emapper
        from ssign_lib.resources import effective_cpu_count as _cpu

        captured = []
        monkeypatch.setattr(subprocess, "run", self._stub_subprocess_run(captured))
        self._touch_annotations(tmp_path, "s")
        # Force dbmem=False so the cmd path doesn't depend on host RAM.
        run_emapper(str(tmp_path / "p.faa"), "/db", "s", str(tmp_path), threads=None, dbmem=False)

        assert len(captured) == 1
        cmd = captured[0]
        i = cmd.index("--cpu")
        assert cmd[i + 1] == str(_cpu())

    def test_explicit_threads_passes_through(self, monkeypatch, tmp_path):
        import subprocess

        from run_eggnog import run_emapper

        captured = []
        monkeypatch.setattr(subprocess, "run", self._stub_subprocess_run(captured))
        self._touch_annotations(tmp_path, "s")
        run_emapper(str(tmp_path / "p.faa"), "/db", "s", str(tmp_path), threads=7, dbmem=False)

        cmd = captured[0]
        i = cmd.index("--cpu")
        assert cmd[i + 1] == "7"
