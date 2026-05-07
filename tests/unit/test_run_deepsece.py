"""Regression + comprehensive tests for run_deepsece.py.

Three layers:

1. parse_deepsece_output: Critical Bug Fix #4 (T3SS preservation),
   column rename, delimiter auto-detection (comma + tab).
2. Model-label contract: SS_MAP / PREDICTED_LABELS pinning.
3. Checkpoint management: _validate_checkpoint, _download_with_retries,
   _ensure_checkpoint, env-var override.

The full run_deepsece() inference path requires torch + ESM + a 2.5 GB
checkpoint and is exercised by tests/integration/test_run_deepsece_integration.py.
"""

import os
import sys

import pytest

# Production module

import run_deepsece as rd  # noqa: E402
from _helpers import DSE_RAW_FIELDS, make_dse_raw_row, write_tsv  # noqa: E402
from run_deepsece import (  # noqa: E402
    MIN_CHECKPOINT_BYTES,
    PREDICTED_LABELS,
    SS_MAP,
    _download_with_retries,
    _ensure_checkpoint,
    _validate_checkpoint,
    parse_deepsece_output,
)


def _write_dse_csv(tmp_dir, rows):
    return write_tsv(os.path.join(tmp_dir, "dse.csv"), DSE_RAW_FIELDS, rows, delimiter=",")


class TestT3ssRowsPreserved:
    def test_t3ss_row_present_after_parse(self, tmp_dir):
        path = _write_dse_csv(
            tmp_dir,
            [
                make_dse_raw_row("GENE_0001", ss_type="T1SS", max_prob=0.95, T1_prob=0.95),
                make_dse_raw_row("GENE_0002", ss_type="T3SS", max_prob=0.91, T3_prob=0.91),
                make_dse_raw_row("GENE_0003", ss_type="Non-secreted", nonsec_prob=0.99),
            ],
        )
        by_locus = {e["locus_tag"]: e for e in parse_deepsece_output(path)}
        assert by_locus["GENE_0002"]["dse_ss_type"] == "T3SS"
        assert by_locus["GENE_0002"]["dse_max_prob"] == "0.91"

    def test_all_ss_type_calls_pass_through(self, tmp_dir):
        rows = [
            make_dse_raw_row(f"GENE_{i:04d}", ss_type=ss_type, max_prob=0.9)
            for i, ss_type in enumerate(["Non-secreted", "T1SS", "T2SS", "T3SS", "T4SS", "T6SS"])
        ]
        entries = parse_deepsece_output(_write_dse_csv(tmp_dir, rows))
        assert len(entries) == 6
        assert {e["dse_ss_type"] for e in entries} == {
            "Non-secreted",
            "T1SS",
            "T2SS",
            "T3SS",
            "T4SS",
            "T6SS",
        }


class TestColumnMapping:
    def test_protein_id_renamed_to_locus_tag(self, tmp_dir):
        path = _write_dse_csv(tmp_dir, [make_dse_raw_row("BIMENO_04457", ss_type="T1SS")])
        entry = parse_deepsece_output(path)[0]
        assert entry["locus_tag"] == "BIMENO_04457"
        assert "protein_id" not in entry

    def test_t_prob_columns_renamed_with_dse_prefix(self, tmp_dir):
        path = _write_dse_csv(tmp_dir, [make_dse_raw_row("GENE_0001", T3_prob=0.7, T6_prob=0.3)])
        entry = parse_deepsece_output(path)[0]
        assert entry["dse_T3_prob"] == "0.7"
        assert entry["dse_T6_prob"] == "0.3"


class TestModelLabelContract:
    """If the upstream model retrains and label order shifts, every downstream
    T3SS guard breaks silently. Pin the contract."""

    def test_predicted_labels_position_3_is_t3ss(self):
        assert PREDICTED_LABELS[3] == "III"
        assert SS_MAP["III"] == "T3SS"

    def test_ss_map_covers_every_predicted_label(self):
        for label in PREDICTED_LABELS:
            assert label in SS_MAP, f"PREDICTED_LABEL {label!r} missing from SS_MAP"

    def test_ss_map_t_indices_match_t_prob_columns(self):
        # dse_T3_prob means "T3SS probability" downstream — pin the alignment.
        assert SS_MAP["I"] == "T1SS"
        assert SS_MAP["II"] == "T2SS"
        assert SS_MAP["III"] == "T3SS"
        assert SS_MAP["IV"] == "T4SS"
        assert SS_MAP["VI"] == "T6SS"


# ---------------------------------------------------------------------------
# Delimiter auto-detection (the bug 4.c surfaced; fixed in this commit)
# ---------------------------------------------------------------------------


class TestDelimiterDetection:
    """parse_deepsece_output tries comma first, then tab. Without the
    locus_tag-non-empty guard, feeding a tab file silently emitted rows
    with empty fields instead of falling through to the tab pass."""

    def test_comma_separated_input(self, tmp_dir):
        path = _write_dse_csv(tmp_dir, [make_dse_raw_row("GENE_0001", ss_type="T1SS")])
        entry = parse_deepsece_output(path)[0]
        assert entry["locus_tag"] == "GENE_0001"

    def test_tab_separated_input(self, tmp_dir):
        path = write_tsv(
            os.path.join(tmp_dir, "dse.tsv"),
            DSE_RAW_FIELDS,
            [make_dse_raw_row("GENE_0001", ss_type="T1SS")],
            delimiter="\t",
        )
        entry = parse_deepsece_output(path)[0]
        assert entry["locus_tag"] == "GENE_0001"

    def test_empty_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "empty.csv")
        with open(path, "w") as f:
            f.write("")
        assert parse_deepsece_output(path) == []

    def test_header_only_file_returns_empty(self, tmp_dir):
        path = os.path.join(tmp_dir, "header_only.csv")
        with open(path, "w") as f:
            f.write(",".join(DSE_RAW_FIELDS) + "\n")
        assert parse_deepsece_output(path) == []


# ---------------------------------------------------------------------------
# _validate_checkpoint — file size sanity
# ---------------------------------------------------------------------------


class TestValidateCheckpoint:
    def test_missing_file_returns_false(self, tmp_dir):
        assert _validate_checkpoint(os.path.join(tmp_dir, "nope.pt")) is False

    def test_truncated_file_removed_and_false(self, tmp_dir):
        path = os.path.join(tmp_dir, "truncated.pt")
        with open(path, "wb") as f:
            f.write(b"\x00" * 100)  # well below MIN_CHECKPOINT_BYTES
        assert _validate_checkpoint(path) is False
        assert not os.path.exists(path)  # truncated file deleted

    def test_full_size_file_returns_true(self, tmp_dir):
        path = os.path.join(tmp_dir, "good.pt")
        with open(path, "wb") as f:
            # Write just over the floor — sparse file is fine for size check
            f.seek(MIN_CHECKPOINT_BYTES + 1024)
            f.write(b"\x00")
        assert _validate_checkpoint(path) is True


# ---------------------------------------------------------------------------
# _download_with_retries — urllib mocking
# ---------------------------------------------------------------------------


class TestDownloadWithRetries:
    def _fake_urlretrieve_succeed(self, dest_size):
        """Builds a urlretrieve fake that writes a file of `dest_size` bytes."""

        def _impl(url, dest, reporthook=None):
            with open(dest, "wb") as f:
                f.seek(dest_size - 1)
                f.write(b"\x00")

        return _impl

    def test_first_attempt_succeeds(self, tmp_dir, monkeypatch):
        dest = os.path.join(tmp_dir, "checkpoint.pt")
        monkeypatch.setattr(
            rd.urllib.request,
            "urlretrieve",
            self._fake_urlretrieve_succeed(MIN_CHECKPOINT_BYTES + 1024),
        )
        assert _download_with_retries("http://example.com/x.pt", dest) is True
        assert os.path.exists(dest)

    def test_retry_after_truncated_download(self, tmp_dir, monkeypatch):
        # First attempt produces a too-small file; second produces full size.
        attempts = []
        full = self._fake_urlretrieve_succeed(MIN_CHECKPOINT_BYTES + 1024)

        def _impl(url, dest, reporthook=None):
            attempts.append(dest)
            if len(attempts) == 1:
                with open(dest, "wb") as f:
                    f.write(b"\x00" * 100)  # truncated
            else:
                full(url, dest, reporthook)

        monkeypatch.setattr(rd.urllib.request, "urlretrieve", _impl)
        # Skip the sleep in retry loop
        monkeypatch.setattr(rd.time, "sleep", lambda *_a, **_k: None)
        dest = os.path.join(tmp_dir, "checkpoint.pt")
        assert _download_with_retries("http://example.com/x.pt", dest) is True
        assert len(attempts) == 2

    def test_all_attempts_fail_returns_false(self, tmp_dir, monkeypatch):
        def _impl(url, dest, reporthook=None):
            raise OSError("network down")

        monkeypatch.setattr(rd.urllib.request, "urlretrieve", _impl)
        monkeypatch.setattr(rd.time, "sleep", lambda *_a, **_k: None)
        assert (
            _download_with_retries(
                "http://example.com/x.pt",
                os.path.join(tmp_dir, "checkpoint.pt"),
            )
            is False
        )


# ---------------------------------------------------------------------------
# _ensure_checkpoint — caching + URL fallback + error path
# ---------------------------------------------------------------------------


class TestEnsureCheckpoint:
    def test_user_provided_path_validated_and_returned(self, tmp_dir, monkeypatch):
        # Pre-stage a valid checkpoint at user path → no download triggered
        path = os.path.join(tmp_dir, "user.pt")
        with open(path, "wb") as f:
            f.seek(MIN_CHECKPOINT_BYTES + 1024)
            f.write(b"\x00")

        called = {"download": 0}

        def fake_download(url, dest):
            called["download"] += 1
            return True

        monkeypatch.setattr(rd, "_download_with_retries", fake_download)
        result = _ensure_checkpoint(path)
        assert result == path
        assert called["download"] == 0

    def test_default_path_cached_checkpoint_used(self, tmp_dir, monkeypatch):
        # Stage a valid checkpoint at the DEFAULT path; ensure no download.
        default = os.path.join(tmp_dir, "default.pt")
        with open(default, "wb") as f:
            f.seek(MIN_CHECKPOINT_BYTES + 1024)
            f.write(b"\x00")
        monkeypatch.setattr(rd, "DEFAULT_CHECKPOINT", default)

        called = {"download": 0}

        def fake_download(url, dest):
            called["download"] += 1
            return True

        monkeypatch.setattr(rd, "_download_with_retries", fake_download)
        result = _ensure_checkpoint()
        assert result == default
        assert called["download"] == 0

    def test_first_url_fails_second_succeeds(self, tmp_dir, monkeypatch):
        default = os.path.join(tmp_dir, "default.pt")
        monkeypatch.setattr(rd, "DEFAULT_CHECKPOINT", default)
        monkeypatch.setattr(rd, "CHECKPOINT_URLS", ["http://first/x.pt", "http://second/x.pt"])

        attempts = []

        def fake_download(url, dest):
            attempts.append(url)
            if "second" in url:
                # Materialise a "valid" checkpoint
                with open(dest, "wb") as f:
                    f.seek(MIN_CHECKPOINT_BYTES + 1024)
                    f.write(b"\x00")
                return True
            return False

        monkeypatch.setattr(rd, "_download_with_retries", fake_download)
        result = _ensure_checkpoint()
        assert result == default
        assert attempts == ["http://first/x.pt", "http://second/x.pt"]

    def test_all_urls_fail_raises_runtime_error(self, tmp_dir, monkeypatch):
        default = os.path.join(tmp_dir, "default.pt")
        monkeypatch.setattr(rd, "DEFAULT_CHECKPOINT", default)
        monkeypatch.setattr(rd, "CHECKPOINT_URLS", ["http://a/x.pt", "http://b/x.pt"])
        monkeypatch.setattr(rd, "_download_with_retries", lambda url, dest: False)

        with pytest.raises(RuntimeError, match="Could not download DeepSecE checkpoint"):
            _ensure_checkpoint()
