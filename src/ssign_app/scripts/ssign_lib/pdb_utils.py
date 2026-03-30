"""PDB parsing utilities for pLDDT extraction, normalization, and validation.

Provides functions to extract mean pLDDT from PDB B-factor columns, normalize
B-factor values from 0-1 scale (HuggingFace ESMFold) to 0-100 scale
(AlphaFold DB / original fair-esm), and validate PDB structure integrity.

The pLDDT scale mismatch between HuggingFace ESMFold (0-1) and AlphaFold DB
(0-100) is the most critical pitfall in the pipeline.  Centralizing
normalization here prevents the bug where all ESMFold structures get filtered
out by the pLDDT threshold (default >=70).
"""

import logging

logger = logging.getLogger(__name__)


def extract_mean_plddt(pdb_string: str, normalize: bool = True) -> float:
    """Extract mean pLDDT from PDB B-factor column.

    Parses ATOM lines from PDB content, extracts B-factor values from
    the fixed-width PDB columns 60:66, and returns the mean value.

    When *normalize* is ``True`` and the maximum B-factor value is <= 1.0,
    all values are multiplied by 100 to convert from HuggingFace ESMFold
    0-1 scale to the standard 0-100 scale used by AlphaFold DB and the
    original fair-esm library.

    Args:
        pdb_string: PDB file contents as a string.
        normalize: If ``True``, detect 0-1 scale and convert to 0-100.

    Returns:
        Mean pLDDT on 0-100 scale (or raw scale if *normalize* is ``False``).
        Returns 0.0 if no ATOM lines are found.
    """
    plddt_values = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                plddt_values.append(bfactor)
            except (ValueError, IndexError):
                pass

    if not plddt_values:
        logger.warning("No ATOM lines found in PDB string; returning pLDDT 0.0")
        return 0.0

    mean_val = sum(plddt_values) / len(plddt_values)

    # Detect HuggingFace ESMFold 0-1 scale and normalize to 0-100
    if normalize and max(plddt_values) <= 1.0:
        logger.debug(
            "Detected 0-1 pLDDT scale (max=%.4f), normalizing to 0-100",
            max(plddt_values),
        )
        mean_val *= 100

    return mean_val


def normalize_pdb_bfactors(pdb_string: str) -> tuple[str, bool]:
    """Normalize PDB B-factor values from 0-1 scale to 0-100 scale.

    Detects if B-factor values are on the 0-1 scale (max <= 1.0) as output
    by HuggingFace ESMFold. If so, multiplies all ATOM line B-factors by 100
    and returns the rewritten PDB string. If already on 0-100 scale (as in
    AlphaFold DB PDB files), returns the original string unchanged.

    Preserves exact PDB fixed-width formatting: B-factor occupies columns
    61-66 (0-indexed 60:66), right-justified, format ``%6.2f``.

    This ensures all PDB files in the pipeline use a consistent 0-100 scale,
    which is important because downstream tools (e.g. Foldseek) may read
    B-factor values directly.

    Args:
        pdb_string: PDB file contents as a string.

    Returns:
        Tuple of ``(pdb_string, was_normalized)``. If *was_normalized* is
        ``True``, B-factors were multiplied by 100.
    """
    # First pass: collect B-factors to detect scale
    bfactors = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                bfactors.append(bfactor)
            except (ValueError, IndexError):
                pass

    if not bfactors:
        logger.warning("No ATOM lines found; cannot normalize B-factors")
        return pdb_string, False

    if max(bfactors) > 1.0:
        # Already on 0-100 scale
        return pdb_string, False

    # Second pass: rewrite B-factors multiplied by 100
    logger.debug(
        "Normalizing B-factors from 0-1 to 0-100 scale (max=%.4f)",
        max(bfactors),
    )
    lines = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") and len(line) >= 66:
            try:
                bfactor = float(line[60:66].strip())
                new_bfactor = bfactor * 100
                # PDB format: B-factor is columns 61-66, right-justified, 6.2f
                line = line[:60] + f"{new_bfactor:6.2f}" + line[66:]
            except (ValueError, IndexError):
                pass
        lines.append(line)

    return "\n".join(lines), True


def validate_pdb_structure(pdb_string: str) -> dict:
    """Perform lightweight validation of a PDB structure string.

    Counts ATOM lines, extracts unique chain IDs, and counts unique
    residues. Used to catch empty or corrupt PDB downloads.

    Args:
        pdb_string: PDB file contents as a string.

    Returns:
        Dictionary with keys:
        - ``atom_count``: Number of ATOM lines.
        - ``chain_ids``: List of unique chain identifiers.
        - ``residue_count``: Number of unique residues (chain + resSeq).
        - ``is_valid``: ``True`` if atom_count > 0.
    """
    atom_count = 0
    chain_ids = set()
    residues = set()

    for line in pdb_string.split("\n"):
        if line.startswith("ATOM") and len(line) > 26:
            atom_count += 1
            chain_id = line[21]
            chain_ids.add(chain_id)
            # Residue identified by chain + residue sequence number (cols 22-26)
            try:
                res_key = (chain_id, line[22:26].strip())
                residues.add(res_key)
            except IndexError:
                pass

    result = {
        "atom_count": atom_count,
        "chain_ids": sorted(chain_ids),
        "residue_count": len(residues),
        "is_valid": atom_count > 0,
    }

    if not result["is_valid"]:
        logger.warning("PDB validation failed: no ATOM lines found")

    return result
