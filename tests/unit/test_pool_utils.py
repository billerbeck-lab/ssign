"""Tests for ssign_app.core._pool_utils.

The load-bearing invariant: for any sample_ids that don't contain the
SEPARATOR and any input records, ``split(pool(records)) == records``.
Tests cover the four pool/split paths (FASTA + TSV) plus the
helper edge cases (collisions across genomes, tags containing the
separator, blank rows, missing id columns).
"""

import csv

import pytest
from hypothesis import given
from hypothesis import strategies as st

from ssign_app.core._pool_utils import (
    SEPARATOR,
    make_prefixed_id,
    pool_fastas,
    pool_tsvs,
    split_fasta_by_source,
    split_prefixed_id,
    split_tsv_by_source,
    validate_sample_id,
)

# ---------------------------------------------------------------------------
# Hypothesis strategies
# ---------------------------------------------------------------------------

# sample_id: realistic ASCII, no SEPARATOR, no whitespace. Single underscore
# is fine (e.g. "ecoli_k12") — the prohibited token is the DOUBLE underscore.
sample_ids = st.text(
    alphabet=st.characters(
        whitelist_categories=("Lu", "Ll", "Nd"),
        whitelist_characters="-.",
    ),
    min_size=1,
    max_size=16,
).filter(lambda s: SEPARATOR not in s)

# Locus tags MAY contain SEPARATOR (e.g. legacy GenBank locus tags). That's
# the whole point of the prefix scheme: only the leading <sample_id>__ is
# stripped at split time.
locus_tags = st.text(
    alphabet=st.characters(
        whitelist_categories=("Lu", "Ll", "Nd"),
        whitelist_characters="-_.",
    ),
    min_size=1,
    max_size=20,
)

sequences = st.text(
    alphabet=st.sampled_from("ACDEFGHIKLMNPQRSTVWY*"),
    min_size=1,
    max_size=80,
)

# TSV cell values: printable text only. Control characters (category Cc,
# which covers tab / CR / LF and NUL) don't round-trip through a plain csv
# TSV, and surrogates (Cs) aren't UTF-8 encodable; real annotation data has
# neither, and this layer tests pool/split, not csv's control-char handling.
tsv_cell_values = st.text(
    alphabet=st.characters(
        blacklist_categories=("Cs", "Cc"),
    ),
    min_size=0,
    max_size=12,
)

# ---------------------------------------------------------------------------
# make_prefixed_id / split_prefixed_id
# ---------------------------------------------------------------------------


class TestMakeAndSplitPrefixedId:
    def test_round_trip_simple(self):
        prefixed = make_prefixed_id("ecoli", "BAKTA_001")
        assert split_prefixed_id(prefixed) == ("ecoli", "BAKTA_001")

    def test_locus_tag_can_contain_separator(self):
        prefixed = make_prefixed_id("ecoli", f"FOO{SEPARATOR}BAR")
        sample_id, tag = split_prefixed_id(prefixed)
        assert sample_id == "ecoli"
        assert tag == f"FOO{SEPARATOR}BAR"

    def test_missing_separator_raises(self):
        with pytest.raises(ValueError, match="prefixed"):
            split_prefixed_id("no-separator-here")

    @given(sample_ids, locus_tags)
    def test_round_trip_property(self, sample_id, locus_tag):
        prefixed = make_prefixed_id(sample_id, locus_tag)
        assert split_prefixed_id(prefixed) == (sample_id, locus_tag)


class TestValidateSampleId:
    def test_empty_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            validate_sample_id("")

    def test_separator_in_sample_id_raises(self):
        with pytest.raises(ValueError, match="separator"):
            validate_sample_id(f"foo{SEPARATOR}bar")

    def test_normal_sample_id_passes(self):
        validate_sample_id("ecoli_k12")  # single underscore is fine
        validate_sample_id("pao1")
        validate_sample_id("genome.1")


# ---------------------------------------------------------------------------
# pool_fastas + split_fasta_by_source
# ---------------------------------------------------------------------------


def _parse_fasta_to_dict(path):
    """Tiny parser local to the test: avoids depending on read_fasta's
    exact tokenisation, so a round-trip failure can't be masked by it."""
    result = {}
    current = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    result[current] = "".join(buf)
                current = line[1:].split()[0]
                buf = []
            elif current is not None:
                buf.append(line.strip())
    if current is not None:
        result[current] = "".join(buf)
    return result


class TestPoolAndSplitFasta:
    def test_pool_then_split_recovers_per_genome(self, tmp_path):
        a = tmp_path / "a.faa"
        a.write_text(">A1\nMKT\n>A2\nMFV\n")
        b = tmp_path / "b.faa"
        b.write_text(">B1\nMGA\n")

        pooled = tmp_path / "pool.faa"
        n = pool_fastas([("ecoli", a), ("pao1", b)], pooled)
        assert n == 3

        pooled_text = pooled.read_text()
        assert f">ecoli{SEPARATOR}A1" in pooled_text
        assert f">ecoli{SEPARATOR}A2" in pooled_text
        assert f">pao1{SEPARATOR}B1" in pooled_text

        split = split_fasta_by_source(pooled, tmp_path / "out")
        assert set(split.keys()) == {"ecoli", "pao1"}

        ecoli = _parse_fasta_to_dict(split["ecoli"])
        pao1 = _parse_fasta_to_dict(split["pao1"])
        assert ecoli == {"A1": "MKT", "A2": "MFV"}
        assert pao1 == {"B1": "MGA"}

    def test_pool_rejects_separator_in_sample_id(self, tmp_path):
        a = tmp_path / "a.faa"
        a.write_text(">A1\nMKT\n")
        with pytest.raises(ValueError, match="separator"):
            pool_fastas([(f"bad{SEPARATOR}id", a)], tmp_path / "pool.faa")

    def test_locus_tags_with_separator_survive_roundtrip(self, tmp_path):
        a = tmp_path / "a.faa"
        a.write_text(f">FOO{SEPARATOR}BAR\nMKT\n")

        pooled = tmp_path / "pool.faa"
        pool_fastas([("g1", a)], pooled)

        split = split_fasta_by_source(pooled, tmp_path / "out")
        recovered = _parse_fasta_to_dict(split["g1"])
        assert recovered == {f"FOO{SEPARATOR}BAR": "MKT"}

    def test_collision_across_genomes(self, tmp_path):
        a = tmp_path / "a.faa"
        a.write_text(">XYZ_001\nMKT\n")
        b = tmp_path / "b.faa"
        b.write_text(">XYZ_001\nMFV\n")

        pooled = tmp_path / "pool.faa"
        pool_fastas([("ga", a), ("gb", b)], pooled)

        split = split_fasta_by_source(pooled, tmp_path / "out")
        assert _parse_fasta_to_dict(split["ga"]) == {"XYZ_001": "MKT"}
        assert _parse_fasta_to_dict(split["gb"]) == {"XYZ_001": "MFV"}

    def test_empty_sources_writes_empty_file(self, tmp_path):
        pooled = tmp_path / "pool.faa"
        n = pool_fastas([], pooled)
        assert n == 0
        assert pooled.exists() and pooled.read_text() == ""

    def test_fasta_header_metadata_preserved(self, tmp_path):
        # FASTA headers often have a description after the ID. The
        # round-trip must preserve everything after the first whitespace.
        a = tmp_path / "a.faa"
        a.write_text(">A1 hypothetical protein\nMKT\n")

        pooled = tmp_path / "pool.faa"
        pool_fastas([("g1", a)], pooled)
        assert "hypothetical protein" in pooled.read_text()

        split = split_fasta_by_source(pooled, tmp_path / "out")
        assert "hypothetical protein" in split["g1"].read_text()

    def test_split_warns_on_unprefixed_records(self, tmp_path, caplog):
        pooled = tmp_path / "pool.faa"
        pooled.write_text(">NOPREFIX\nMKT\n")
        with caplog.at_level("WARNING"):
            paths = split_fasta_by_source(pooled, tmp_path / "out")
        assert paths == {}
        assert "missing prefix" in caplog.text


@given(
    genomes=st.lists(
        st.tuples(
            sample_ids,
            st.dictionaries(locus_tags, sequences, min_size=1, max_size=6),
        ),
        min_size=1,
        max_size=4,
        unique_by=lambda t: t[0],
    )
)
def test_fasta_round_trip_property(genomes, tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("rtfasta")
    sources = []
    for sample_id, id_to_seq in genomes:
        src = tmp_path / f"{sample_id}.faa"
        with open(src, "w") as f:
            for tag, seq in id_to_seq.items():
                f.write(f">{tag}\n{seq}\n")
        sources.append((sample_id, src))

    pooled = tmp_path / "pool.faa"
    pool_fastas(sources, pooled)
    split = split_fasta_by_source(pooled, tmp_path / "out")

    for sample_id, id_to_seq in genomes:
        recovered = _parse_fasta_to_dict(split[sample_id])
        assert recovered == id_to_seq, f"Mismatch for {sample_id}"


# ---------------------------------------------------------------------------
# pool_tsvs + split_tsv_by_source
# ---------------------------------------------------------------------------


def _write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)  # default comma
        w.writeheader()
        w.writerows(rows)


def _read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))  # default comma


class TestPoolAndSplitTsv:
    def test_pool_then_split_recovers_per_genome(self, tmp_path):
        a = tmp_path / "a.tsv"
        _write_tsv(a, ["locus_tag", "x"], [{"locus_tag": "A1", "x": "11"}, {"locus_tag": "A2", "x": "12"}])
        b = tmp_path / "b.tsv"
        _write_tsv(b, ["locus_tag", "x"], [{"locus_tag": "B1", "x": "21"}])

        pooled = tmp_path / "pool.tsv"
        n = pool_tsvs([("g1", a), ("g2", b)], pooled)
        assert n == 3

        split = split_tsv_by_source(pooled, tmp_path / "out")
        a_back = _read_tsv(split["g1"])
        b_back = _read_tsv(split["g2"])
        assert a_back == [{"locus_tag": "A1", "x": "11"}, {"locus_tag": "A2", "x": "12"}]
        assert b_back == [{"locus_tag": "B1", "x": "21"}]

    def test_pool_unions_columns(self, tmp_path):
        # A has [locus_tag, x]; B has [locus_tag, x, y]. Pooled has all three,
        # and A's rows get an empty y.
        a = tmp_path / "a.tsv"
        _write_tsv(a, ["locus_tag", "x"], [{"locus_tag": "A1", "x": "11"}])
        b = tmp_path / "b.tsv"
        _write_tsv(b, ["locus_tag", "x", "y"], [{"locus_tag": "B1", "x": "21", "y": "Q"}])

        pooled = tmp_path / "pool.tsv"
        pool_tsvs([("g1", a), ("g2", b)], pooled)

        rows = _read_tsv(pooled)
        assert {r["locus_tag"] for r in rows} == {f"g1{SEPARATOR}A1", f"g2{SEPARATOR}B1"}
        a_row = next(r for r in rows if r["locus_tag"].endswith("A1"))
        b_row = next(r for r in rows if r["locus_tag"].endswith("B1"))
        assert a_row["y"] == ""
        assert b_row["y"] == "Q"

    def test_pool_skips_blank_locus_tag_rows(self, tmp_path):
        # Mirrors the load_substrate_ids tolerance for blank-tag artefacts.
        a = tmp_path / "a.tsv"
        a.write_text("locus_tag\tx\nA1\t11\n\t99\nA2\t12\n")

        pooled = tmp_path / "pool.tsv"
        n = pool_tsvs([("g1", a)], pooled)
        assert n == 2

    def test_pool_missing_id_column_raises(self, tmp_path):
        a = tmp_path / "a.tsv"
        _write_tsv(a, ["other_col", "x"], [{"other_col": "A1", "x": "11"}])
        with pytest.raises(ValueError, match="locus_tag"):
            pool_tsvs([("g1", a)], tmp_path / "pool.tsv")

    def test_pool_rejects_separator_in_sample_id(self, tmp_path):
        a = tmp_path / "a.tsv"
        _write_tsv(a, ["locus_tag", "x"], [{"locus_tag": "A1", "x": "11"}])
        with pytest.raises(ValueError, match="separator"):
            pool_tsvs([(f"bad{SEPARATOR}id", a)], tmp_path / "pool.tsv")

    def test_split_warns_on_unprefixed_rows(self, tmp_path, caplog):
        pooled = tmp_path / "pool.tsv"
        _write_tsv(pooled, ["locus_tag", "x"], [{"locus_tag": "UNPREFIXED", "x": "11"}])
        with caplog.at_level("WARNING"):
            split = split_tsv_by_source(pooled, tmp_path / "out")
        assert split == {}
        assert "prefix" in caplog.text

    def test_split_custom_id_column(self, tmp_path):
        # Some tools emit seq_id rather than locus_tag.
        pooled = tmp_path / "pool.tsv"
        _write_tsv(
            pooled,
            ["seq_id", "score"],
            [
                {"seq_id": f"g1{SEPARATOR}P1", "score": "0.9"},
                {"seq_id": f"g2{SEPARATOR}P1", "score": "0.4"},
            ],
        )
        split = split_tsv_by_source(pooled, tmp_path / "out", id_column="seq_id")
        assert set(split.keys()) == {"g1", "g2"}
        g1 = _read_tsv(split["g1"])
        assert g1 == [{"seq_id": "P1", "score": "0.9"}]

    def test_collision_across_genomes_tsv(self, tmp_path):
        a = tmp_path / "a.tsv"
        _write_tsv(a, ["locus_tag", "x"], [{"locus_tag": "XYZ_001", "x": "a"}])
        b = tmp_path / "b.tsv"
        _write_tsv(b, ["locus_tag", "x"], [{"locus_tag": "XYZ_001", "x": "b"}])

        pooled = tmp_path / "pool.tsv"
        pool_tsvs([("ga", a), ("gb", b)], pooled)
        split = split_tsv_by_source(pooled, tmp_path / "out")
        assert _read_tsv(split["ga"]) == [{"locus_tag": "XYZ_001", "x": "a"}]
        assert _read_tsv(split["gb"]) == [{"locus_tag": "XYZ_001", "x": "b"}]

    # ------------------------------------------------------------------
    # CSV-format handling — regression for the 2026-06-05 multi-genome
    # crash where _pool_interproscan.csv was read with a tab delimiter
    # and the entire comma-separated header parsed as one field.
    # ------------------------------------------------------------------

    def test_split_csv_input_preserves_csv_output(self, tmp_path):
        """A `.csv` pooled file is read comma-delimited and split files keep `.csv`."""
        pooled = tmp_path / "_pool_interproscan.csv"
        _write_csv(
            pooled,
            ["locus_tag", "interpro_domains", "interpro_pfam_ids"],
            [
                {"locus_tag": f"g1{SEPARATOR}P1", "interpro_domains": "Sec61", "interpro_pfam_ids": "PF00344"},
                {"locus_tag": f"g2{SEPARATOR}P2", "interpro_domains": "TonB", "interpro_pfam_ids": "PF03544"},
            ],
        )
        split = split_tsv_by_source(pooled, tmp_path / "out")
        assert set(split.keys()) == {"g1", "g2"}
        for sid, path in split.items():
            assert path.suffix == ".csv", f"{sid} output should keep .csv extension"

        g1 = _read_csv(split["g1"])
        assert g1 == [{"locus_tag": "P1", "interpro_domains": "Sec61", "interpro_pfam_ids": "PF00344"}]

    def test_pool_tsvs_csv_round_trip(self, tmp_path):
        """Two CSV inputs round-trip through pool + split with comma delimiters end-to-end."""
        a = tmp_path / "a.csv"
        _write_csv(a, ["locus_tag", "x"], [{"locus_tag": "A1", "x": "11"}])
        b = tmp_path / "b.csv"
        _write_csv(b, ["locus_tag", "x"], [{"locus_tag": "B1", "x": "21"}])

        pooled = tmp_path / "pool.csv"
        n = pool_tsvs([("g1", a), ("g2", b)], pooled)
        assert n == 2
        # Pooled file must be readable as comma-delimited.
        pooled_rows = _read_csv(pooled)
        assert {r["locus_tag"] for r in pooled_rows} == {f"g1{SEPARATOR}A1", f"g2{SEPARATOR}B1"}

        split = split_tsv_by_source(pooled, tmp_path / "out")
        assert _read_csv(split["g1"]) == [{"locus_tag": "A1", "x": "11"}]
        assert _read_csv(split["g2"]) == [{"locus_tag": "B1", "x": "21"}]

    def test_pool_handles_mixed_csv_and_tsv_sources(self, tmp_path):
        """If one source is CSV and another is TSV, each is read with its own delimiter.

        Not a realistic ssign workflow (per-tool outputs are uniform) but
        the helper should not assume the sources share a format.
        """
        a = tmp_path / "a.csv"
        _write_csv(a, ["locus_tag", "x"], [{"locus_tag": "A1", "x": "11"}])
        b = tmp_path / "b.tsv"
        _write_tsv(b, ["locus_tag", "x"], [{"locus_tag": "B1", "x": "21"}])

        pooled = tmp_path / "pool.tsv"
        n = pool_tsvs([("g1", a), ("g2", b)], pooled)
        assert n == 2
        pooled_rows = _read_tsv(pooled)
        assert {r["locus_tag"] for r in pooled_rows} == {f"g1{SEPARATOR}A1", f"g2{SEPARATOR}B1"}


@given(
    genomes=st.lists(
        st.tuples(
            sample_ids,
            st.dictionaries(
                locus_tags,
                tsv_cell_values,
                min_size=1,
                max_size=6,
            ),
        ),
        min_size=1,
        max_size=4,
        unique_by=lambda t: t[0],
    )
)
def test_tsv_round_trip_property(genomes, tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("rttsv")
    sources = []
    for sample_id, id_to_val in genomes:
        src = tmp_path / f"{sample_id}.tsv"
        _write_tsv(
            src,
            ["locus_tag", "val"],
            [{"locus_tag": tag, "val": val} for tag, val in id_to_val.items()],
        )
        sources.append((sample_id, src))

    pooled = tmp_path / "pool.tsv"
    pool_tsvs(sources, pooled)
    split = split_tsv_by_source(pooled, tmp_path / "out")

    for sample_id, id_to_val in genomes:
        rows = _read_tsv(split[sample_id])
        recovered = {r["locus_tag"]: r["val"] for r in rows}
        assert recovered == id_to_val, f"Mismatch for {sample_id}"
