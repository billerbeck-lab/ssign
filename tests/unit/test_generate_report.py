"""Text summary report: the distilled enrichment block + coverage header.

Covers generate_report._summarize_enrichment (the COMBINED-track distillation of
the per-type permutation test) and the generate_text_report wiring, without
needing a real pipeline run.
"""

import csv
import os

import generate_report as gr


def _write(path, header, rows, sep="\t"):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter=sep)
        w.writerow(header)
        w.writerows(rows)


_ENR_HEADER = ["ss_type", "tool", "mode", "fold", "qvalue", "significant"]


def _enr_tsv(path):
    # COMBINED rows (the headline track) + noise rows for other tools that must be ignored.
    _write(
        path,
        _ENR_HEADER,
        [
            ["T5aSS", "COMBINED", "self", "6.8", "0.003", "True"],
            ["T5aSS", "COMBINED", "window", "1.5", "0.47", "False"],
            ["T1SS", "COMBINED", "window", "4.9", "0.16", "False"],
            ["T5aSS", "DLP", "self", "17.9", "0.012", "True"],  # non-COMBINED: ignored
            ["T1SS", "SignalP", "window", "0.3", "0.95", "False"],  # ignored
        ],
    )


def test_summarize_enrichment_labels_sort_and_significance(tmp_dir):
    p = os.path.join(tmp_dir, "e.tsv")
    _enr_tsv(p)
    out = gr._summarize_enrichment(p, n_genomes=1)
    body = "\n".join(out)

    # Only the 3 COMBINED rows, DLP/SignalP dropped.
    assert sum(1 for line in out if "x   q=" in line) == 3
    # T5aSS has a self row, so its window row is the hitchhiker; T1SS stays bare.
    assert "T5aSS (self)" in body
    assert "T5aSS (hitchhiker)" in body
    assert "T1SS " in body and "T1SS (" not in body
    # Sorted by q ascending: the significant self row is first.
    first_data = next(line for line in out if "x   q=" in line)
    assert "T5aSS (self)" in first_data and first_data.rstrip().endswith("*")
    # Non-significant rows carry no star.
    assert not any(line.rstrip().endswith("*") for line in out if "T1SS " in line)
    assert any("Benjamini-Hochberg q < 0.05" in line for line in out)


def test_summarize_enrichment_power_note_when_nothing_significant(tmp_dir):
    p = os.path.join(tmp_dir, "e.tsv")
    _write(p, _ENR_HEADER, [["T1SS", "COMBINED", "window", "2.0", "0.30", "False"]])
    single = gr._summarize_enrichment(p, n_genomes=1)
    pooled = gr._summarize_enrichment(p, n_genomes=5)
    assert any("single-genome power is low" in line for line in single)
    assert any("nothing reached q < 0.05" in line for line in pooled)
    assert not any("single-genome" in line for line in pooled)


def test_summarize_enrichment_degenerate_inputs(tmp_dir):
    empty = os.path.join(tmp_dir, "empty.tsv")
    open(empty, "w").close()
    assert gr._summarize_enrichment(empty) == []
    # header only, no COMBINED track
    nocomb = os.path.join(tmp_dir, "nc.tsv")
    _write(nocomb, _ENR_HEADER, [["T1SS", "DLP", "window", "3.0", "0.2", "False"]])
    assert gr._summarize_enrichment(nocomb) == []
    assert gr._summarize_enrichment(os.path.join(tmp_dir, "nope.tsv")) == []


def test_generate_text_report_coverage_header_and_enrichment(tmp_dir):
    master = os.path.join(tmp_dir, "m.csv")
    _write(
        master,
        ["locus_tag", "sample_id", "nearby_ss_types", "signalp_prediction", "eggnog_desc"],
        [
            ["a1", "G", "T5aSS", "SP", "adhesin"],
            ["a2", "G", "T1SS", "SP", ""],
        ],
        sep=",",
    )
    enr = os.path.join(tmp_dir, "e.tsv")
    _enr_tsv(enr)
    out = os.path.join(tmp_dir, "summary.txt")
    gr.generate_text_report([master], enr, out, tier="full")
    text = open(out).read()

    assert "Annotation coverage:" in text
    assert "(of 2 secreted proteins)" not in text  # dropped
    assert "SignalP" in text and "2/2 (100%)" in text
    assert "Secretion-system enrichment (per-type permutation test, COMBINED predictor):" in text
    assert "Enrichment results:" not in text  # old raw-dump header is gone
