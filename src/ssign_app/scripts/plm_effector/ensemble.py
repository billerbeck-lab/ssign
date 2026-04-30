"""Two-stage ensemble prediction for PLM-Effector.

Refactored from upstream's `ensemble_step12_4predict.py` (CC-BY 3.0, see
package __init__.py). Upstream shipped this as a script that loads
per-terminal `.npz` files from disk; this module consumes the in-memory
feature dict produced by `feature_extraction.extract_all_features` and
returns per-protein final probabilities.

Upstream bug fixed in the refactor:
- `ensemble_step12_4predict.py:172` referenced `args.data_path`, which was
  never registered with argparse — would raise AttributeError on first
  run. Here the weights directory is a function argument, removing the
  ambiguity.
"""

from __future__ import annotations

import gc
import os

import numpy as np
import torch

from .inference import loadmodel_4predict, loadmodel_4test
from .utils import pool_features, set_seed

# Per-effector XGBoost decision thresholds. Values come directly from
# Zheng et al. 2026's Table 4 / supplementary — best F1 threshold
# calibrated on their held-out evaluation set.
_BEST_THRESHOLDS = {
    "T1SE": 0.5,
    "T2SE": 0.7,
    "T3SE": 0.7,
    "T4SE": 0.6,
    "T6SE": 0.5,
}

# Which effector/model combinations need mean-pooled embeddings vs raw
# per-residue embeddings. The rules come from the upstream code — each
# neural network architecture was chosen per (effector, PLM, terminal).
_POOL_RULES = {
    # (effector_type, model_index, terminal) : bool  (True = mean-pool)
    ("T1SE", 1, "Nterminal"): True,  # model1 (esm1, Nter, MLP)
    ("T1SE", 2, "Nterminal"): True,  # model2 (esm2_t33, Nter, MLP)
    ("T3SE", 2, "Nterminal"): True,  # model2 (esm2_t33, Nter, MLP)
    ("T4SE", 6, "Cterminal"): True,  # model6 (ProtT5, Cter, MLP)
    ("T6SE", 4, "Cterminal"): True,  # model4 (esm1, Cter, MLP)
}


def _maybe_pool(embedding, attention_masks, effector_type, model_idx, terminal):
    """Apply mean-pooling only when the upstream architecture expects it."""
    if _POOL_RULES.get((effector_type, model_idx, terminal), False):
        return pool_features(embedding, attention_masks, pooling="mean").numpy()
    return embedding


def _n_base_models(effector_type: str) -> int:
    """T4SE has 8 base models (includes ProtBert); the others have 6."""
    return 8 if effector_type == "T4SE" else 6


def _iter_model_specs(effector_type: str):
    """Yield `(model_idx, pretrained_type, terminal)` tuples in the upstream order.

    Matches the docstring of upstream `ensemble_step12_4predict.py`:
        model1: esm1 Nterminal
        model2: esm2_t33 Nterminal
        model3: ProtT5 Nterminal
        model4: esm1 Cterminal
        model5: esm2_t33 Cterminal
        model6: ProtT5 Cterminal
        model7: ProtBert Nterminal (T4SE only)
        model8: ProtBert Cterminal (T4SE only)
    """
    core_specs = [
        (1, "esm1", "Nterminal"),
        (2, "esm2_t33", "Nterminal"),
        (3, "ProtT5", "Nterminal"),
        (4, "esm1", "Cterminal"),
        (5, "esm2_t33", "Cterminal"),
        (6, "ProtT5", "Cterminal"),
    ]
    yield from core_specs
    if effector_type == "T4SE":
        yield (7, "ProtBert", "Nterminal")
        yield (8, "ProtBert", "Cterminal")


def run_ensemble(
    features: dict,
    weights_dir: str,
    effector_type: str,
    device,
    n_folds: int = 5,
) -> tuple:
    """Run the two-stage ensemble on extracted features.

    Args:
        features: nested dict from
            `feature_extraction.extract_all_features`.
        weights_dir: directory containing
            `trained_models/{effector_type}_model{i}_fold{f}.pth` and
            `trained_models/{effector_type}_XGB_stackingmeta_model.json`.
        effector_type: one of `T1SE`, `T2SE`, `T3SE`, `T4SE`, `T6SE`.
        device: torch device.
        n_folds: number of per-fold checkpoints to average. Upstream uses 5.

    Returns:
        `(seq_ids, stacked_probs, final_probs, passes_threshold)` as numpy
        arrays. `stacked_probs` has shape `(n_proteins, n_base_models)`.
    """
    set_seed(42)
    n_models = _n_base_models(effector_type)
    trained_models_dir = os.path.join(weights_dir, "trained_models")
    xgb_path = os.path.join(
        trained_models_dir, f"{effector_type}_XGB_stackingmeta_model.json"
    )
    if not os.path.exists(xgb_path):
        raise FileNotFoundError(
            f"XGBoost stacking model not found: {xgb_path}\n"
            f"  Expected PLM-Effector weights at "
            f"{weights_dir}/trained_models/{effector_type}_XGB_stackingmeta_model.json"
        )

    oof_probs = {f"model{i + 1}": [] for i in range(n_models)}
    seq_ids_by_model = {}

    for model_idx, pretrained_type, terminal in _iter_model_specs(effector_type):
        feat = features[pretrained_type][terminal]
        seq_ids_by_model[model_idx] = feat["seq_ids"]

        x = _maybe_pool(
            feat["embedding"],
            feat["attention_masks"],
            effector_type,
            model_idx,
            terminal,
        )
        x_tensor = torch.from_numpy(x).float()

        for fold in range(n_folds):
            ckpt_name = f"{effector_type}_model{model_idx}_fold{fold}.pth"
            if model_idx == 1:
                # Upstream uses loadmodel_4test for model1 (returns both preds
                # and probs); we only need probs but keep the call to preserve
                # the random-state side-effects that may matter for exact
                # reproducibility with the paper results.
                _, probs = loadmodel_4test(
                    trained_models_dir, ckpt_name, x_tensor, device
                )
            else:
                probs = loadmodel_4predict(
                    trained_models_dir, ckpt_name, x_tensor, device
                )
            oof_probs[f"model{model_idx}"].append(probs)

        del x_tensor
        if device.type == "cuda":
            torch.cuda.empty_cache()
        gc.collect()

    # Verify all base models agree on the protein ordering
    ref_ids = seq_ids_by_model[1]
    for i, ids in seq_ids_by_model.items():
        if not np.array_equal(ref_ids, ids):
            raise ValueError(
                f"Protein ordering disagrees between model1 and model{i}; "
                f"this indicates a bug in feature extraction."
            )

    # Average across folds to get one probability per base model per protein
    stacked = []
    for i in range(n_models):
        fold_probs = np.stack(oof_probs[f"model{i + 1}"], axis=0)
        stacked.append(fold_probs.mean(axis=0).reshape(-1, 1))
    x_stacking = np.hstack(stacked)

    # Second stage: XGBoost stacking meta-model
    from xgboost import XGBClassifier

    meta = XGBClassifier()
    meta.load_model(xgb_path)
    final_probs = meta.predict_proba(x_stacking)[:, 1]

    threshold = _BEST_THRESHOLDS[effector_type]
    passes = final_probs >= threshold

    return ref_ids, x_stacking, final_probs, passes
