"""Tests for plm_effector.predict_api.

`write_predictions_tsv` is already covered by the existing
`tests/unit/test_run_plm_effector.py::TestWritePredictionsTsv` class.
This file fills in the rest:

- `_resolve_device`: torch.device coercion + CUDA-availability gate.
- `predict`: argument validation paths (effector_type / missing input /
  missing weights). The orchestration body (extract_all_features →
  run_ensemble → write) requires GPU + 15 GB of weights and is
  exercised by `tests/integration/test_run_plm_effector_integration.py`.
"""

import importlib.util
import os
import sys
from types import SimpleNamespace

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src", "ssign_app", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

# predict_api imports torch lazily — this module-level import is fine
# without torch installed. Tests that exercise _resolve_device or the
# full orchestration mark themselves needs_torch below.
from plm_effector.predict_api import (  # noqa: E402
    _VALID_EFFECTOR_TYPES,
    _resolve_device,
    predict,
)

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None,
    reason="torch not installed (covered by tests/integration/)",
)

# ---------------------------------------------------------------------------
# _resolve_device
# ---------------------------------------------------------------------------
# Mock torch lazily — we only need a stand-in `device` constructor and
# `cuda.is_available`. The real torch is fine to import (it's in the dev
# extras), but constructing torch.device + flipping CUDA availability
# without an actual GPU needs monkeypatching.


@needs_torch
class TestResolveDevice:
    def test_cpu_returns_cpu_device(self):
        import torch

        result = _resolve_device("cpu")
        assert result == torch.device("cpu")

    def test_passthrough_torch_device(self):
        import torch

        device = torch.device("cpu")
        assert _resolve_device(device) is device

    def test_unrecognised_string_raises_value_error(self):
        with pytest.raises(ValueError, match="Unrecognised device"):
            _resolve_device("tpu")

    def test_cuda_when_unavailable_raises_runtime_error(self, monkeypatch):
        import torch

        monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
        with pytest.raises(RuntimeError, match="no CUDA device"):
            _resolve_device("cuda")

    def test_default_none_treated_as_cuda(self, monkeypatch):
        import torch

        # Same gate as "cuda" — default value picks GPU when available
        monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
        with pytest.raises(RuntimeError, match="no CUDA device"):
            _resolve_device(None)

    def test_cuda_when_available_returns_cuda_device(self, monkeypatch):
        import torch

        monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
        result = _resolve_device("cuda")
        assert result.type == "cuda"


# ---------------------------------------------------------------------------
# predict() — argument validation only (orchestration is integration-tested)
# ---------------------------------------------------------------------------


class TestPredictValidation:
    def test_invalid_effector_type_raises_value_error(self, tmp_dir):
        # Create real input + weights paths so we hit the effector_type check
        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)
        with pytest.raises(ValueError, match="effector_type must be one of"):
            predict(
                proteins_fasta=input_fasta,
                weights_dir=weights,
                effector_type="T9SE",
                out_path=os.path.join(tmp_dir, "out.tsv"),
            )

    def test_missing_input_raises_file_not_found(self, tmp_dir):
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)
        with pytest.raises(FileNotFoundError, match="Input FASTA not found"):
            predict(
                proteins_fasta=os.path.join(tmp_dir, "does_not_exist.faa"),
                weights_dir=weights,
                effector_type="T1SE",
                out_path=os.path.join(tmp_dir, "out.tsv"),
            )

    def test_missing_weights_dir_raises_file_not_found(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        with pytest.raises(FileNotFoundError, match="weights directory not found"):
            predict(
                proteins_fasta=input_fasta,
                weights_dir=os.path.join(tmp_dir, "no_weights"),
                effector_type="T1SE",
                out_path=os.path.join(tmp_dir, "out.tsv"),
            )

    @needs_torch
    @pytest.mark.parametrize("effector_type", _VALID_EFFECTOR_TYPES)
    def test_every_canonical_effector_type_accepted(self, tmp_dir, effector_type, monkeypatch):
        """All five canonical effector types must pass validation. We mock the
        downstream extract → ensemble pipeline and verify only that validation
        let the call through to write_predictions_tsv."""
        import numpy as np
        import plm_effector.predict_api as predict_api

        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)

        # Stub the lazy imports inside predict()
        fake_extract = SimpleNamespace(
            extract_all_features=lambda **_kw: {"features": "irrelevant"},
        )
        # T4SE has 8 base models, others have 6
        n_base = 8 if effector_type == "T4SE" else 6
        fake_ensemble = SimpleNamespace(
            run_ensemble=lambda **_kw: (
                np.array([">P1"]),
                np.zeros((1, n_base)),
                np.array([0.5]),
                np.array([True]),
            ),
        )
        monkeypatch.setitem(
            sys.modules,
            "plm_effector.feature_extraction",
            fake_extract,
        )
        monkeypatch.setitem(
            sys.modules,
            "plm_effector.ensemble",
            fake_ensemble,
        )
        # _resolve_device pulls torch — keep it deterministic
        import torch

        monkeypatch.setattr(predict_api, "_resolve_device", lambda _d: torch.device("cpu"))

        out = os.path.join(tmp_dir, "out.tsv")
        n_positive = predict(
            proteins_fasta=input_fasta,
            weights_dir=weights,
            effector_type=effector_type,
            out_path=out,
            device="cpu",
        )
        assert n_positive == 1
        # Output written with the expected effector_type tag in the last column
        with open(out) as f:
            lines = f.readlines()
        assert lines[1].rstrip("\n").split("\t")[-1] == effector_type
