"""Integration test for run_hhsuite.py.

Skips cleanly when HH-suite databases aren't available — the
Soeding lab's GWDG mirror has been dead since at least 2026-04-29,
and no public mirror of Pfam/PDB70/UniRef30 in HH-suite format is
currently accessible. Once a mirror is restored (or we set up our
own), this test should run end-to-end.

Required env vars (skip-on-missing):
    SSIGN_HHSUITE_PFAM        path to pfam_v36 directory (.ffindex/.ffdata)
    SSIGN_HHSUITE_PDB70       path to pdb70_from_mmcif directory
    SSIGN_HHSUITE_UNICLUST    path to uniclust30 / UniRef30 directory

Run with:
    SSIGN_HHSUITE_PFAM=/path/to/pfam \\
    SSIGN_HHSUITE_PDB70=/path/to/pdb70 \\
    SSIGN_HHSUITE_UNICLUST=/path/to/uniclust30 \\
    pytest -m integration tests/integration/test_run_hhsuite_integration.py
"""

import os
import shutil

import pytest


pytestmark = pytest.mark.integration


def _skip_unless_hhsuite_ready():
    if shutil.which("hhblits") is None or shutil.which("hhsearch") is None:
        pytest.skip(
            "hh-suite (hhblits/hhsearch) not on PATH; "
            "install via `micromamba install -c bioconda hhsuite`"
        )

    pfam = os.environ.get("SSIGN_HHSUITE_PFAM")
    pdb70 = os.environ.get("SSIGN_HHSUITE_PDB70")
    uniclust = os.environ.get("SSIGN_HHSUITE_UNICLUST")
    # HH-suite DBs are passed by *prefix* (e.g. .../pdb70 → loads
    # pdb70_a3m.ffdata, pdb70_hhm.ffdata, ...). Validate by checking the
    # required `_a3m.ffdata` companion file rather than the bare prefix.
    missing = [
        name
        for name, path in [
            ("SSIGN_HHSUITE_PFAM", pfam),
            ("SSIGN_HHSUITE_PDB70", pdb70),
            ("SSIGN_HHSUITE_UNICLUST", uniclust),
        ]
        if not path or not os.path.exists(f"{path}_a3m.ffdata")
    ]
    if missing:
        pytest.skip(
            f"HH-suite DBs unavailable ({', '.join(missing)}). "
            f"Set <var>=<DB-prefix-without-suffix> for each DB."
        )
    return pfam, pdb70, uniclust


class TestRunHhsuite:
    def test_runs_on_fixture_proteins(
        self, tmp_dir, t1ss_fixture_proteins
    ):
        pfam, pdb70, uniclust = _skip_unless_hhsuite_ready()

        import sys
        sys.path.insert(
            0,
            os.path.abspath(
                os.path.join(
                    os.path.dirname(__file__), "..", "..",
                    "src", "ssign_app", "scripts",
                )
            ),
        )
        from run_hhsuite import run_hhsuite_parallel
        from ssign_lib.fasta_io import read_fasta

        sequences = read_fasta(t1ss_fixture_proteins)
        # Subset to BIMENO_04457 only — full run is ~5 min/protein on
        # CPU, so 9 proteins = 45 min. One protein keeps the test ~5 min.
        target = "BIMENO_04457"
        if target not in sequences:
            pytest.skip(f"{target} missing from fixture")
        sequences = {target: sequences[target]}

        output_dir = os.path.join(tmp_dir, "hhsuite_out")
        os.makedirs(output_dir)

        results = run_hhsuite_parallel(
            sequences=sequences,
            pfam_db=pfam,
            pdb70_db=pdb70,
            uniclust_db=uniclust,
            output_dir=output_dir,
            max_workers=1,
            cpu_per_job=4,
        )

        assert isinstance(results, dict)
        assert target in results, (
            f"hhsearch returned no result row for {target} — wrapper "
            f"may have dropped it on a per-protein failure."
        )
        # Wrapper returns dict[locus_tag] -> {"pfam_top1_*", "pdb_top1_*"}
        # — at minimum we expect some hit info structure to be present
        # for the autotransporter against Pfam.
        target_result = results[target]
        assert isinstance(target_result, dict)
