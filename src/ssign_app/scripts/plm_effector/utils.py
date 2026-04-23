"""Prediction-time utilities for PLM-Effector.

Vendored subset of the upstream `utils.py` (Zheng et al. 2026,
https://github.com/zhengdd0422/PLM-Effector). See `__init__.py` for
licence and attribution. Training-only helpers (WeightedSampler,
compute_class_weights, custom_collate, load_model_without_dataparallel,
load_val_numpy, load_test_numpy, load_test_numpy_nopool,
process_Cterminal_sequence) have been dropped because they are not
needed for inference. Upstream `torch_geometric`, `scipy.sparse`, and
`Bio.PDB` imports were vestigial in the prediction path and have been
removed.
"""

from __future__ import annotations

import random
import re

import numpy as np
import torch


_BERT_LIKE = {"Bert", "BioBERT", "ProtBert", "ProtT5"}


def _normalise_sequence(seq: str, model_type: str) -> str:
    """Prepare a single protein sequence for a pretrained tokenizer.

    BERT-family tokenisers (ProtBert, ProtT5) expect space-separated
    residues; ESM tokenisers expect the raw string. All models want the
    rare residues U/Z/O/B replaced by X.
    """
    spaced = " ".join(seq) if model_type in _BERT_LIKE else seq
    return re.sub(r"[UZOB]", "X", spaced)


def read_fasta_for_prediction(fasta_path: str, model_type: str = "esm1"):
    """Return `(ids, sequences)` parsed from a FASTA file.

    IDs include the leading `>` as emitted by the upstream parser, which
    callers strip. Residues are normalised for the chosen pretrained
    tokenizer.
    """
    ids: list[str] = []
    sequences: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                ids.append(line)
            else:
                sequences.append(_normalise_sequence(line, model_type))
    return ids, sequences


def read_fasta_for_prediction_terminal(
    fasta_path: str,
    model_type: str = "esm1",
    terminal: str = "Cterminal",
    maxlen: int = 1022,
):
    """Like `read_fasta_for_prediction` but truncated to N- or C-terminal region.

    Long proteins get cut to the first `maxlen` residues (N-terminal) or the
    last `maxlen` residues (C-terminal). Short proteins pass through
    unchanged. PLM-Effector trains one model per terminal.
    """
    ids: list[str] = []
    sequences: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                ids.append(line)
            else:
                truncated = line[-maxlen:] if terminal == "Cterminal" else line[:maxlen]
                sequences.append(_normalise_sequence(truncated, model_type))
    return ids, sequences


def batch_extract_features(
    sequences,
    pretrained_type: str,
    model,
    tokenizer,
    device,
    max_length: int = 512,
    batch_size: int = 10,
):
    """Run a batched forward pass through a pretrained PLM.

    Returns `(features, attention_masks)` as numpy arrays shaped
    `(n_sequences, max_length, hidden_dim)` and
    `(n_sequences, max_length)` respectively.
    """
    features: list[torch.Tensor] = []
    attention_masks: list[torch.Tensor] = []

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i : i + batch_size]
        inputs = tokenizer(
            batch,
            return_tensors="pt",
            padding="max_length",
            truncation=True,
            max_length=max_length,
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            if pretrained_type == "ProtT5":
                outputs = model(
                    input_ids=inputs["input_ids"],
                    attention_mask=inputs["attention_mask"],
                )
            else:
                outputs = model(**inputs)
        features.append(outputs.last_hidden_state.cpu())
        attention_masks.append(inputs["attention_mask"].cpu())
        del inputs, outputs
        if device.type == "cuda":
            torch.cuda.empty_cache()

    all_features = torch.cat(features, dim=0).numpy()
    all_attention_masks = torch.cat(attention_masks, dim=0).numpy()
    return all_features, all_attention_masks


def pool_features(features, attention_masks, pooling: str = "mean") -> torch.Tensor:
    """Pool per-residue features into a single per-sequence vector.

    Supports `"mean"` (default) and `"max"` pooling, masking padding
    positions via the attention mask.
    """
    if not torch.is_tensor(features):
        features = torch.tensor(features, dtype=torch.float32)
    if not torch.is_tensor(attention_masks):
        attention_masks = torch.tensor(attention_masks, dtype=torch.float32)

    attention_masks = attention_masks.unsqueeze(-1)

    if pooling == "mean":
        return (features * attention_masks).sum(dim=1) / attention_masks.sum(dim=1)
    if pooling == "max":
        return (features * attention_masks).max(dim=1).values
    raise ValueError(f"Unsupported pooling type {pooling!r}; use 'mean' or 'max'.")


def set_seed(seed: int) -> None:
    """Seed Python, NumPy, and PyTorch RNGs for reproducible inference."""
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
