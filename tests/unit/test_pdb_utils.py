"""Tests for ssign_lib.pdb_utils — PDB parsing, pLDDT normalization, validation.

These tests cover the three most critical bug fixes in the pipeline:
1. pLDDT scale detection (0-1 vs 0-100) — auto-normalization
2. PDB line length guard (>= 66, NOT > 66)
3. B-factor rewriting preserves PDB column widths
"""

import pytest
from ssign_lib.pdb_utils import (
    extract_mean_plddt,
    normalize_pdb_bfactors,
    validate_pdb_structure,
)


class TestExtractMeanPlddt:
    def test_01_scale_detected_and_normalized(self, sample_pdb_01_scale):
        """Values on 0-1 scale should be detected and multiplied by 100."""
        mean = extract_mean_plddt(sample_pdb_01_scale, normalize=True)
        # (0.85 + 0.72 + 0.90 + 0.68) / 4 * 100 = 78.75
        assert abs(mean - 78.75) < 0.01

    def test_100_scale_passthrough(self, sample_pdb_100_scale):
        """Values on 0-100 scale should NOT be re-scaled."""
        mean = extract_mean_plddt(sample_pdb_100_scale, normalize=True)
        # (85 + 72 + 90 + 68) / 4 = 78.75
        assert abs(mean - 78.75) < 0.01

    def test_normalize_false_returns_raw(self, sample_pdb_01_scale):
        """With normalize=False, 0-1 values are returned as-is."""
        mean = extract_mean_plddt(sample_pdb_01_scale, normalize=False)
        assert abs(mean - 0.7875) < 0.001

    def test_no_atom_lines_returns_zero(self):
        pdb = "REMARK  test file\nEND\n"
        assert extract_mean_plddt(pdb) == 0.0

    def test_exact_66_char_line(self, pdb_exact_66_chars):
        """Line of exactly 66 chars must be parsed (>= 66 guard, not > 66)."""
        mean = extract_mean_plddt(pdb_exact_66_chars)
        assert abs(mean - 75.0) < 0.01

    def test_short_line_ignored(self):
        """Lines shorter than 66 chars should be skipped, not crash."""
        pdb = "ATOM      1  N   ALA A   1       1.000   2.000   3.000\nEND\n"
        mean = extract_mean_plddt(pdb)
        assert mean == 0.0


class TestNormalizePdbBfactors:
    def test_01_scale_normalized(self, sample_pdb_01_scale):
        normalized, was_changed = normalize_pdb_bfactors(sample_pdb_01_scale)
        assert was_changed is True
        # After normalization, extracting should give 0-100 values
        mean = extract_mean_plddt(normalized, normalize=False)
        assert abs(mean - 78.75) < 0.01

    def test_100_scale_unchanged(self, sample_pdb_100_scale):
        normalized, was_changed = normalize_pdb_bfactors(sample_pdb_100_scale)
        assert was_changed is False
        assert normalized == sample_pdb_100_scale

    def test_preserves_column_widths(self, sample_pdb_01_scale):
        """B-factor must stay in columns 60:66 (6 chars, right-justified)."""
        normalized, _ = normalize_pdb_bfactors(sample_pdb_01_scale)
        for line in normalized.split("\n"):
            if line.startswith("ATOM") and len(line) >= 66:
                bfactor_str = line[60:66]
                assert len(bfactor_str) == 6
                # Should be parseable as float
                float(bfactor_str)

    def test_no_atoms_returns_unchanged(self):
        pdb = "REMARK  test\nEND\n"
        result, changed = normalize_pdb_bfactors(pdb)
        assert changed is False
        assert result == pdb


class TestValidatePdbStructure:
    def test_valid_structure(self, sample_pdb_100_scale):
        result = validate_pdb_structure(sample_pdb_100_scale)
        assert result["is_valid"] is True
        assert result["atom_count"] == 4
        assert result["residue_count"] == 1  # All same residue (ALA A 1)
        assert "A" in result["chain_ids"]

    def test_empty_pdb(self):
        result = validate_pdb_structure("REMARK test\nEND\n")
        assert result["is_valid"] is False
        assert result["atom_count"] == 0
