"""Feature extraction via pretrained protein language models.

Refactored from upstream's `transformers_pretrainedmodel-features_extract_4predict.py`
(https://github.com/zhengdd0422/PLM-Effector, CC-BY 3.0 — see package
__init__.py). Upstream shipped the logic as a script with argparse and a
file-based `.npz` intermediate protocol; this module exposes a callable
function that returns in-memory numpy arrays instead. The heavy lifting
(tokenise → batch forward pass) is unchanged.
"""

from __future__ import annotations

import gc
import os

import numpy as np
import torch

from .utils import (
    batch_extract_features,
    read_fasta_for_prediction,
    read_fasta_for_prediction_terminal,
)

_MODEL_SUBDIRS = {
    "ProtBert": "prot_bert",
    "ProtT5": "prot_t5_xl_uniref50",
    "esm1": "esm1b_t33_650M_UR50S",
    "esm2_t33": "esm2_t33_650M_UR50D",
}

_MODEL_MAX_LENGTH = {"ProtBert": 512, "ProtT5": 512, "esm1": 1024, "esm2_t33": 1024}
_MODEL_REAL_SEQUENCE_LEN = {
    "ProtBert": 510,
    "ProtT5": 511,
    "esm1": 1022,
    "esm2_t33": 1022,
}


def _load_model_and_tokenizer(pretrained_type: str, weights_dir: str, device):
    """Load a HuggingFace pretrained PLM from `weights_dir/transformers_pretrained/<subdir>/`."""
    from transformers import AutoModel, AutoTokenizer, T5EncoderModel, T5Tokenizer

    subdir = _MODEL_SUBDIRS[pretrained_type]
    model_path = os.path.join(weights_dir, "transformers_pretrained", subdir)
    if not os.path.isdir(model_path):
        raise FileNotFoundError(
            f"Pretrained weights directory not found: {model_path}\n"
            f"  Expected PLM-Effector's HuggingFace weights at "
            f"{weights_dir}/transformers_pretrained/{subdir}/"
        )

    if pretrained_type == "ProtBert":
        tokenizer = AutoTokenizer.from_pretrained(model_path, do_lower_case=False)
        model = AutoModel.from_pretrained(model_path)
    elif pretrained_type == "ProtT5":
        tokenizer = T5Tokenizer.from_pretrained(model_path, do_lower_case=False)
        model = T5EncoderModel.from_pretrained(model_path)
    else:
        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModel.from_pretrained(model_path)

    model = model.to(device)
    model.eval()
    return model, tokenizer


def extract_terminal_features(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
) -> dict:
    """Extract N- and C-terminal features for every protein in the FASTA.

    Returns a dict keyed by terminal (`"Nterminal"` / `"Cterminal"`) with
    values `{"embedding": ndarray, "attention_masks": ndarray, "seq_ids": ndarray}`
    matching the shape that the upstream `.npz` files used.
    """
    model, tokenizer = _load_model_and_tokenizer(pretrained_type, weights_dir, device)
    max_length = _MODEL_MAX_LENGTH[pretrained_type]
    real_seq_len = _MODEL_REAL_SEQUENCE_LEN[pretrained_type]

    out: dict = {}
    try:
        for terminal in ("Nterminal", "Cterminal"):
            if terminal == "Cterminal":
                ids, sequences = read_fasta_for_prediction_terminal(
                    proteins_fasta,
                    model_type=pretrained_type,
                    terminal=terminal,
                    maxlen=real_seq_len,
                )
                tokenizer.padding_side = "left"
            else:
                ids, sequences = read_fasta_for_prediction(
                    proteins_fasta, model_type=pretrained_type
                )
                tokenizer.padding_side = "right"

            embeddings, attention_masks = batch_extract_features(
                sequences,
                pretrained_type,
                model,
                tokenizer,
                device,
                max_length=max_length,
                batch_size=batch_size,
            )

            out[terminal] = {
                "embedding": embeddings,
                "attention_masks": attention_masks,
                "seq_ids": np.array(ids),
            }

            if device.type == "cuda":
                torch.cuda.empty_cache()
            gc.collect()
    finally:
        del model, tokenizer
        if device.type == "cuda":
            torch.cuda.empty_cache()
        gc.collect()

    return out


def extract_all_features(
    proteins_fasta: str,
    effector_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
) -> dict:
    """Run feature extraction for every PLM that `effector_type` needs.

    Returns a nested dict: `features[pretrained_type][terminal] = {"embedding", "attention_masks", "seq_ids"}`.
    """
    if effector_type == "T4SE":
        pretrained_types = ["esm1", "esm2_t33", "ProtBert", "ProtT5"]
    else:
        pretrained_types = ["esm1", "esm2_t33", "ProtT5"]

    features: dict = {}
    for pretrained_type in pretrained_types:
        features[pretrained_type] = extract_terminal_features(
            proteins_fasta=proteins_fasta,
            pretrained_type=pretrained_type,
            weights_dir=weights_dir,
            device=device,
            batch_size=batch_size,
        )
    return features
