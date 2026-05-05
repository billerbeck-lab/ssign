"""Tests for extract_proteins.py.

The biggest module in Group B — entry point for parsing genome
annotations into ssign's canonical gene_info schema. Four format paths
(GenBank / GFF3+FASTA / nucleotide FASTA via Pyrodigal / protein FASTA)
plus organism-name inference and locus_tag deduplication.

Pyrodigal-driven contig prediction is integration-test territory
(needs the gene-finder model trained on real DNA) and is exercised by
tests/integration/test_pipeline_fixture.py. Here we cover the
GenBank / GFF3 / .faa paths, the dispatch logic, the organism-inference
chain, and the dedup/FASTA-write steps.
"""

import json
import os
import sys

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

from _helpers import read_tsv_rows, run_script_main  # noqa: E402
from extract_proteins import (  # noqa: E402
    extract_from_genbank,
    extract_from_gff3,
)
from extract_proteins import (
    main as extract_main,
)

# ---------------------------------------------------------------------------
# GenBank fixture builder
# ---------------------------------------------------------------------------


def _make_cds(start, end, strand, qualifiers):
    """Build a CDS SeqFeature with the given location + qualifier dict."""
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers=qualifiers,
    )


def _write_genbank(path, contig_id, cds_list, organism="", source_organism=""):
    """Write a small GenBank file with the given CDS features."""
    seq_len = max((cds.location.end for cds in cds_list), default=1000)
    record = SeqRecord(
        Seq("A" * seq_len),
        id=contig_id,
        name=contig_id,
        annotations={"molecule_type": "DNA"},
    )
    if organism:
        record.annotations["organism"] = organism
    if source_organism:
        record.features.append(
            SeqFeature(
                FeatureLocation(0, seq_len, strand=1),
                type="source",
                qualifiers={"organism": [source_organism]},
            )
        )
    record.features.extend(cds_list)
    with open(path, "w") as f:
        SeqIO.write([record], f, "genbank")
    return path


# ---------------------------------------------------------------------------
# extract_from_genbank
# ---------------------------------------------------------------------------


class TestExtractFromGenbank:
    def test_basic_cds_extracted(self, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKTLLLT"],
                        "product": ["test protein"],
                    },
                )
            ],
        )
        entries = list(extract_from_genbank(gbff, "sample"))
        assert len(entries) == 1
        e = entries[0]
        assert e["locus_tag"] == "GENE_001"
        assert e["product"] == "test protein"
        assert e["sequence"] == "MKTLLLT"
        assert e["contig"] == "contig_1"
        assert e["strand"] == "+"

    @pytest.mark.parametrize(
        "qualifiers",
        [
            # Pseudogene markers (both qualifier names trigger the skip)
            {"locus_tag": ["X"], "translation": ["MKT"], "pseudo": [""]},
            {"locus_tag": ["X"], "translation": ["MKT"], "pseudogene": ["unknown"]},
            # Missing translation
            {"locus_tag": ["X"]},
            # Missing both locus_tag and protein_id (no usable identifier)
            {"translation": ["MKT"]},
        ],
        ids=["pseudo", "pseudogene", "no-translation", "no-id"],
    )
    def test_unsuitable_cds_skipped(self, tmp_dir, qualifiers):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [_make_cds(0, 30, 1, qualifiers)],
        )
        assert list(extract_from_genbank(gbff, "sample")) == []

    def test_locus_tag_fallback_to_protein_id(self, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "protein_id": ["WP_000000001.1"],
                        "translation": ["MKT"],
                    },
                )
            ],
        )
        entries = list(extract_from_genbank(gbff, "sample"))
        assert len(entries) == 1
        assert entries[0]["locus_tag"] == "WP_000000001.1"
        assert entries[0]["protein_id"] == "WP_000000001.1"

    def test_minus_strand_decoded(self, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    -1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                    },
                )
            ],
        )
        entries = list(extract_from_genbank(gbff, "sample"))
        assert entries[0]["strand"] == "-"

    def test_default_product_when_missing(self, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                    },
                )
            ],
        )
        entries = list(extract_from_genbank(gbff, "sample"))
        assert entries[0]["product"] == "hypothetical protein"

    def test_non_cds_features_skipped(self, tmp_dir):
        gbff_path = os.path.join(tmp_dir, "test.gbff")
        record = SeqRecord(
            Seq("A" * 1000),
            id="contig_1",
            name="contig_1",
            annotations={"molecule_type": "DNA"},
        )
        # tRNA feature — must not be extracted as a CDS
        record.features.append(
            SeqFeature(
                FeatureLocation(0, 30, strand=1),
                type="tRNA",
                qualifiers={"locus_tag": ["TRNA_001"]},
            )
        )
        record.features.append(
            _make_cds(
                100,
                130,
                1,
                {
                    "locus_tag": ["GENE_001"],
                    "translation": ["MKT"],
                },
            )
        )
        with open(gbff_path, "w") as f:
            SeqIO.write([record], f, "genbank")
        entries = list(extract_from_genbank(gbff_path, "sample"))
        assert {e["locus_tag"] for e in entries} == {"GENE_001"}


# ---------------------------------------------------------------------------
# extract_from_gff3
# ---------------------------------------------------------------------------


def _write_gff3_pair(tmp_dir, contig_seq, cds_records):
    """cds_records: list of (start_1based, end, strand, attrs_dict)."""
    fasta_path = os.path.join(tmp_dir, "genome.fna")
    with open(fasta_path, "w") as f:
        f.write(">contig_1\n")
        f.write(contig_seq + "\n")
    gff_path = os.path.join(tmp_dir, "genome.gff3")
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
        for start, end, strand, attrs in cds_records:
            attr_str = ";".join(f"{k}={v}" for k, v in attrs.items())
            f.write(f"contig_1\t.\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attr_str}\n")
    return gff_path, fasta_path


class TestExtractFromGff3:
    def test_one_based_to_zero_based_start_conversion(self, tmp_dir):
        # GFF3 start=1 → entry start=0
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [(1, 30, "+", {"locus_tag": "GENE_001"})],
        )
        entries = list(extract_from_gff3(gff, fasta, "sample"))
        assert entries[0]["start"] == 0
        assert entries[0]["end"] == 30

    def test_url_decoded_product(self, tmp_dir):
        # GFF3 percent-encodes spaces and commas in attribute values
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [
                (
                    1,
                    30,
                    "+",
                    {
                        "locus_tag": "GENE_001",
                        "product": "type%20IV%20secretion%2C%20component",
                    },
                )
            ],
        )
        entries = list(extract_from_gff3(gff, fasta, "sample"))
        assert entries[0]["product"] == "type IV secretion, component"

    def test_translation_attribute_used_when_present(self, tmp_dir):
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [
                (
                    1,
                    30,
                    "+",
                    {
                        "locus_tag": "GENE_001",
                        "translation": "MPROTEIN",
                    },
                )
            ],
        )
        entries = list(extract_from_gff3(gff, fasta, "sample"))
        assert entries[0]["sequence"] == "MPROTEIN"

    def test_minus_strand_uses_revcomp(self, tmp_dir):
        # ATGAAATAA on minus strand → revcomp = TTATTTCAT → translate "LFH"
        # We don't pin the exact translation; we pin the property that the
        # resulting translation differs from the plus-strand translation.
        seq = "ATGAAATAA" + "AAA" * 10
        gff_plus, fasta = _write_gff3_pair(
            tmp_dir,
            seq,
            [(1, 9, "+", {"locus_tag": "GENE_PLUS"})],
        )
        plus_entries = list(extract_from_gff3(gff_plus, fasta, "sample"))

        gff_minus, _ = _write_gff3_pair(
            tmp_dir,
            seq,
            [(1, 9, "-", {"locus_tag": "GENE_MINUS"})],
        )
        minus_entries = list(extract_from_gff3(gff_minus, fasta, "sample"))

        # Plus-strand translation starts with M (start codon ATG); minus does not
        assert plus_entries[0]["sequence"].startswith("M")
        assert plus_entries[0]["sequence"] != minus_entries[0]["sequence"]

    def test_locus_tag_fallback_to_id(self, tmp_dir):
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [
                (
                    1,
                    30,
                    "+",
                    {
                        "ID": "cds-XYZ",
                        "translation": "MPROTEIN",
                    },
                )
            ],
        )
        entries = list(extract_from_gff3(gff, fasta, "sample"))
        assert entries[0]["locus_tag"] == "cds-XYZ"

    def test_non_cds_lines_skipped(self, tmp_dir):
        fasta_path = os.path.join(tmp_dir, "genome.fna")
        with open(fasta_path, "w") as f:
            f.write(">contig_1\nATG" + "AAA" * 10 + "\n")
        gff_path = os.path.join(tmp_dir, "genome.gff3")
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig_1\t.\tgene\t1\t30\t.\t+\t.\tID=gene_001\n")
            f.write("contig_1\t.\tCDS\t1\t30\t.\t+\t0\tlocus_tag=GENE_001;translation=MKT\n")
        entries = list(extract_from_gff3(gff_path, fasta_path, "sample"))
        assert {e["locus_tag"] for e in entries} == {"GENE_001"}

    def test_comment_lines_skipped(self, tmp_dir):
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [(1, 30, "+", {"locus_tag": "GENE_001", "translation": "MKT"})],
        )
        # _write_gff3_pair already prepends "##gff-version 3"; verify it's tolerated
        entries = list(extract_from_gff3(gff, fasta, "sample"))
        assert len(entries) == 1


# ---------------------------------------------------------------------------
# main() integration — dispatch, organism inference, dedup, output formats
# ---------------------------------------------------------------------------


def _run_main(monkeypatch, tmp_dir, input_path, *, fasta="", original_filename="", emit_metadata=False):
    out_proteins = os.path.join(tmp_dir, "proteins.faa")
    out_gene_info = os.path.join(tmp_dir, "gene_info.tsv")
    argv = [
        "extract_proteins",
        "--input",
        input_path,
        "--sample",
        "test_sample",
        "--out-proteins",
        out_proteins,
        "--out-gene-info",
        out_gene_info,
    ]
    if fasta:
        argv += ["--fasta", fasta]
    if original_filename:
        argv += ["--original-filename", original_filename]
    out_metadata = None
    if emit_metadata:
        out_metadata = os.path.join(tmp_dir, "metadata.json")
        argv += ["--out-metadata", out_metadata]
    run_script_main(monkeypatch, extract_main, argv)
    return out_proteins, out_gene_info, out_metadata


class TestMainGenbankDispatch:
    def test_genbank_input_extracted(self, monkeypatch, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                    },
                )
            ],
            organism="Xanthomonas campestris",
        )
        _, gene_info, _ = _run_main(monkeypatch, tmp_dir, gbff)
        rows = read_tsv_rows(gene_info)
        assert {r["locus_tag"] for r in rows} == {"GENE_001"}

    def test_organism_extracted_from_genbank_record(self, monkeypatch, tmp_dir):
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                    },
                )
            ],
            organism="Xanthomonas campestris",
        )
        _, _, metadata = _run_main(monkeypatch, tmp_dir, gbff, emit_metadata=True)
        with open(metadata) as f:
            meta = json.load(f)
        assert meta["organism"] == "Xanthomonas campestris"

    def test_genus_only_record_enriched_from_filename(self, monkeypatch, tmp_dir):
        # Record has only "Xanthomonas" (genus); filename has full binomial.
        gbff = _write_genbank(
            os.path.join(tmp_dir, "Xanthomonas_campestris_genomic.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                    },
                )
            ],
            organism="Xanthomonas",
        )
        _, _, metadata = _run_main(monkeypatch, tmp_dir, gbff, emit_metadata=True)
        with open(metadata) as f:
            meta = json.load(f)
        assert meta["organism"] == "Xanthomonas campestris"

    def test_dedup_by_locus_tag_keeps_first(self, monkeypatch, tmp_dir):
        # Two CDS with the same locus_tag — only first is kept
        gbff = _write_genbank(
            os.path.join(tmp_dir, "test.gbff"),
            "contig_1",
            [
                _make_cds(
                    0,
                    30,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MKT"],
                        "product": ["first"],
                    },
                ),
                _make_cds(
                    100,
                    130,
                    1,
                    {
                        "locus_tag": ["GENE_001"],
                        "translation": ["MFV"],
                        "product": ["duplicate"],
                    },
                ),
            ],
        )
        _, gene_info, _ = _run_main(monkeypatch, tmp_dir, gbff)
        rows = read_tsv_rows(gene_info)
        assert len(rows) == 1
        assert rows[0]["product"] == "first"


class TestMainGff3Dispatch:
    def test_gff3_with_fasta_pair(self, monkeypatch, tmp_dir):
        gff, fasta = _write_gff3_pair(
            tmp_dir,
            "ATG" + "AAA" * 10,
            [(1, 30, "+", {"locus_tag": "GENE_001", "translation": "MKT"})],
        )
        _, gene_info, _ = _run_main(monkeypatch, tmp_dir, gff, fasta=fasta)
        rows = read_tsv_rows(gene_info)
        assert {r["locus_tag"] for r in rows} == {"GENE_001"}

    def test_gff3_without_fasta_exits_with_error(self, monkeypatch, tmp_dir):
        gff_path = os.path.join(tmp_dir, "genome.gff3")
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
        with pytest.raises(SystemExit) as exc_info:
            _run_main(monkeypatch, tmp_dir, gff_path)
        assert exc_info.value.code == 1


class TestMainProteinFastaDispatch:
    def test_faa_input_read_directly(self, monkeypatch, tmp_dir):
        faa = os.path.join(tmp_dir, "proteins.faa")
        with open(faa, "w") as f:
            f.write(">P1 some product\nMKTLLLT\n")
            f.write(">P2\nMFVFLVL\n")
        _, gene_info, _ = _run_main(monkeypatch, tmp_dir, faa)
        rows = read_tsv_rows(gene_info)
        by_locus = {r["locus_tag"]: r for r in rows}
        assert set(by_locus.keys()) == {"P1", "P2"}
        # Description is preserved as the product when present
        assert "some product" in by_locus["P1"]["product"]


class TestMainErrorPaths:
    def test_unsupported_extension_exits(self, monkeypatch, tmp_dir):
        path = os.path.join(tmp_dir, "thing.unknown")
        with open(path, "w") as f:
            f.write("nope\n")
        with pytest.raises(SystemExit) as exc_info:
            _run_main(monkeypatch, tmp_dir, path)
        assert exc_info.value.code == 1
