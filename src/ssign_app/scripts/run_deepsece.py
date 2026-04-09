#!/usr/bin/env python3
"""Run DeepSecE and parse output to standard format.

DeepSecE predicts if proteins are secreted and by which secretion system type
(T1SS, T2SS, T3SS, T4SS, T6SS, or non-secreted). Uses a fine-tuned ESM-1b
transformer model via the Python API (no CLI exists).

Checkpoint auto-downloaded on first run to ~/.ssign/models/deepsece_checkpoint.pt

Column mapping from pipeline/scripts/parse_deepsece.py.
Output columns: locus_tag, dse_ss_type, dse_max_prob, plus per-type probabilities.
"""

import argparse
import csv
import logging
import os
import sys
import tempfile
import time
import urllib.request

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# DeepSecE class labels (index order from the model)
PREDICTED_LABELS = ["-", "I", "II", "III", "IV", "VI"]
SS_MAP = {
    "-": "Non-secreted", "I": "T1SS", "II": "T2SS",
    "III": "T3SS", "IV": "T4SS", "VI": "T6SS",
}

CHECKPOINT_URL = "https://tool2-mml.sjtu.edu.cn/DeepSecE/checkpoint.pt"
# No official mirror exists; the SJTU server is the only source.
CHECKPOINT_URLS = [CHECKPOINT_URL]
DEFAULT_CHECKPOINT_DIR = os.path.join(os.path.expanduser("~"), ".ssign", "models")
DEFAULT_CHECKPOINT = os.path.join(DEFAULT_CHECKPOINT_DIR, "deepsece_checkpoint.pt")
# Expected checkpoint size: ~2.5 GB (includes ESM-1b weights + classifier).
# Reject files smaller than 100 MB as truncated.
MIN_CHECKPOINT_BYTES = 100 * 1024 * 1024
DOWNLOAD_TIMEOUT_SEC = 120
DOWNLOAD_MAX_RETRIES = 3

# Column name mapping from raw DeepSecE output to standardized names
_COLUMN_MAP = {
    "protein_id": "locus_tag",
    "deepsece_prediction": "deepsece_prediction",
    "deepsece_ss_type": "dse_ss_type",
    "max_prob": "dse_max_prob",
    "nonsec_prob": "dse_nonsec_prob",
    "T1_prob": "dse_T1_prob",
    "T2_prob": "dse_T2_prob",
    "T3_prob": "dse_T3_prob",
    "T4_prob": "dse_T4_prob",
    "T6_prob": "dse_T6_prob",
}


def _validate_checkpoint(path):
    """Check that a downloaded checkpoint file is not truncated."""
    if not os.path.exists(path):
        return False
    size = os.path.getsize(path)
    if size < MIN_CHECKPOINT_BYTES:
        logger.warning(
            f"Checkpoint file is only {size:,} bytes — likely truncated "
            f"(expected ≥{MIN_CHECKPOINT_BYTES:,} bytes). Removing."
        )
        os.remove(path)
        return False
    return True


def _download_with_retries(url, dest):
    """Download *url* to *dest* with retries and exponential back-off.

    Returns True on success, False on failure.
    """
    for attempt in range(1, DOWNLOAD_MAX_RETRIES + 1):
        logger.info(f"  Attempt {attempt}/{DOWNLOAD_MAX_RETRIES}: {url}")
        partial = dest + ".part"
        try:
            def _progress(block_num, block_size, total_size):
                if total_size > 0:
                    pct = min(100, block_num * block_size * 100 // total_size)
                    if block_num % 50 == 0:
                        logger.info(f"  Download: {pct}%")

            # Set a global socket timeout so hung connections don't block forever
            urllib.request.urlretrieve(
                url, partial, reporthook=_progress,
            )

            # Atomically move into place only if the size looks right
            if os.path.exists(partial):
                size = os.path.getsize(partial)
                if size < MIN_CHECKPOINT_BYTES:
                    logger.warning(
                        f"  Downloaded file is only {size:,} bytes "
                        f"(expected ≥{MIN_CHECKPOINT_BYTES:,}). Retrying..."
                    )
                    os.remove(partial)
                    raise RuntimeError("Truncated download")
                os.replace(partial, dest)
                return True

        except Exception as e:
            logger.warning(f"  Download attempt {attempt} failed: {e}")
            if os.path.exists(partial):
                os.remove(partial)
            if attempt < DOWNLOAD_MAX_RETRIES:
                wait = 2 ** attempt  # 2 s, 4 s, 8 s
                logger.info(f"  Waiting {wait}s before retry...")
                time.sleep(wait)

    return False


def _ensure_checkpoint(checkpoint_path=None):
    """Ensure the DeepSecE checkpoint file exists, downloading if needed."""
    if checkpoint_path and os.path.exists(checkpoint_path):
        if _validate_checkpoint(checkpoint_path):
            return checkpoint_path

    path = checkpoint_path or DEFAULT_CHECKPOINT
    if os.path.exists(path) and _validate_checkpoint(path):
        return path

    # Download checkpoint — try each URL in turn
    os.makedirs(os.path.dirname(path), exist_ok=True)
    logger.info(f"Downloading DeepSecE checkpoint (~2.5 GB) to {path} ...")

    import socket
    old_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(DOWNLOAD_TIMEOUT_SEC)

    try:
        for url in CHECKPOINT_URLS:
            if _download_with_retries(url, path):
                logger.info(f"  Checkpoint saved to {path}")
                return path
    finally:
        socket.setdefaulttimeout(old_timeout)

    # All URLs / retries exhausted
    logger.error(
        "All download attempts failed. The SJTU server "
        "(tool2-mml.sjtu.edu.cn) may be temporarily unavailable."
    )
    logger.error(
        "Please download the checkpoint manually and place it at:\n"
        f"  {path}\n"
        "\n"
        "Using wget (recommended — it handles retries automatically):\n"
        f"  wget -c --tries=5 --timeout=60 {CHECKPOINT_URL} -O \"{path}\"\n"
        "\n"
        "Or using curl:\n"
        f"  curl -L --retry 5 --connect-timeout 60 -o \"{path}\" {CHECKPOINT_URL}"
    )
    raise RuntimeError(
        f"Could not download DeepSecE checkpoint after "
        f"{DOWNLOAD_MAX_RETRIES} attempts per URL.\n"
        f"  Common causes:\n"
        f"    - SJTU server (tool2-mml.sjtu.edu.cn) is temporarily down\n"
        f"    - Network/firewall blocking the download\n"
        f"    - Slow connection causing timeouts (currently {DOWNLOAD_TIMEOUT_SEC}s)\n"
        f"  How to fix:\n"
        f"    - Download manually with wget: wget -c --tries=5 {CHECKPOINT_URL} -O \"{path}\"\n"
        f"    - Then re-run this script (it will find the cached checkpoint)"
    )


def run_deepsece(input_fasta, output_dir, checkpoint_path=None, batch_size=1):
    """Run DeepSecE prediction using the Python API directly.

    Returns path to the output CSV.
    """
    import argparse
    import numpy as np

    # FRAGILE: torch import requires PyTorch installation (large dependency)
    # If this breaks: pip install torch (CPU-only: pip install torch --index-url https://download.pytorch.org/whl/cpu)
    try:
        import torch
        import torch.serialization
    except ImportError as e:
        raise RuntimeError(
            f"PyTorch not installed: {e}\n"
            f"  Common causes:\n"
            f"    - torch package is not installed in this environment\n"
            f"  How to fix:\n"
            f"    - pip install torch\n"
            f"    - For CPU-only (smaller): pip install torch --index-url https://download.pytorch.org/whl/cpu"
        ) from e

    # PyTorch 2.6+ changed default to weights_only=True, but ESM checkpoints
    # contain argparse.Namespace which isn't in the safe globals list.
    # Allow it before any model loading happens.
    try:
        torch.serialization.add_safe_globals([argparse.Namespace])
    except AttributeError:
        pass  # Older PyTorch without add_safe_globals

    # Memory optimization: monkey-patch torch.load to use mmap=True
    # ESM1b file is 7.3GB on disk — mmap avoids loading it all into RAM at once.
    # Critical for systems with <=8GB RAM.
    _orig_torch_load = torch.load

    def _mmap_torch_load(*args, **kwargs):
        kwargs.setdefault("mmap", True)
        try:
            return _orig_torch_load(*args, **kwargs)
        except RuntimeError:
            # ESM-1b uses old pickle format which doesn't support mmap
            kwargs.pop("mmap", None)
            return _orig_torch_load(*args, **kwargs)

    torch.load = _mmap_torch_load

    # FRAGILE: DeepSecE model import requires the DeepSecE package
    # If this breaks: pip install DeepSecE
    try:
        from DeepSecE.model import EffectorTransformer
    except ImportError as e:
        torch.load = _orig_torch_load
        raise RuntimeError(
            f"DeepSecE package not installed: {e}\n"
            f"  Common causes:\n"
            f"    - DeepSecE is not installed in this environment\n"
            f"  How to fix:\n"
            f"    - pip install DeepSecE\n"
            f"    - Or: pip install git+https://github.com/SijinHuang/DeepSecE.git"
        ) from e

    # FRAGILE: ESM (Evolutionary Scale Modeling) import requires fair-esm package
    # If this breaks: pip install fair-esm
    try:
        from esm import Alphabet, FastaBatchedDataset
    except ImportError as e:
        torch.load = _orig_torch_load
        raise RuntimeError(
            f"ESM (fair-esm) package not installed: {e}\n"
            f"  Common causes:\n"
            f"    - fair-esm is not installed in this environment\n"
            f"  How to fix:\n"
            f"    - pip install fair-esm"
        ) from e

    torch.load = _orig_torch_load

    from torch.utils.data import DataLoader

    # Ensure checkpoint exists
    checkpoint = _ensure_checkpoint(checkpoint_path)

    # Device setup
    if torch.cuda.is_available():
        device = torch.device("cuda")
        logger.info(f"Running DeepSecE on GPU: {torch.cuda.get_device_name(0)}")
    else:
        device = torch.device("cpu")
        logger.info("Running DeepSecE on CPU (this may be slow for large datasets)")

    # Build and load model
    # ESM1b (~7.3GB on disk) loaded via mmap to avoid RAM spike.
    # DeepSecE checkpoint (~2.6GB) loaded on top.
    logger.info("Loading DeepSecE model (mmap mode for low-memory systems)...")
    import gc

    # FRAGILE: Model construction loads ESM-1b weights (~7.3 GB) — can fail with MemoryError
    # If this breaks: ensure at least 10 GB free RAM, or use a machine with more memory
    torch.load = _mmap_torch_load
    try:
        model = EffectorTransformer(
            emb_dim=1280,
            repr_layer=33,
            hid_dim=256,
            num_layers=1,
            heads=4,
            dropout_rate=0.4,
            num_classes=6,
            return_attn=False,
        )
    except MemoryError as e:
        torch.load = _orig_torch_load
        raise RuntimeError(
            f"Out of memory loading DeepSecE model: {e}\n"
            f"  Common causes:\n"
            f"    - Not enough RAM (ESM-1b requires ~10 GB free)\n"
            f"  How to fix:\n"
            f"    - Close other applications to free RAM\n"
            f"    - Use a machine with at least 16 GB RAM\n"
            f"    - Set batch_size=1 to reduce memory usage during inference"
        ) from e
    finally:
        torch.load = _orig_torch_load

    # Load checkpoint — use mmap + map_location to avoid device mismatch
    logger.info(f"Loading checkpoint from {checkpoint}...")
    # FRAGILE: Checkpoint loading can fail with MemoryError or corrupted file
    # If this breaks: re-download the checkpoint or ensure sufficient RAM
    try:
        state_dict = torch.load(
            checkpoint, map_location=device, weights_only=True, mmap=True
        )
        model.load_state_dict(state_dict, strict=False)
        del state_dict
        gc.collect()
    except MemoryError as e:
        raise RuntimeError(
            f"Out of memory loading DeepSecE checkpoint: {e}\n"
            f"  Common causes:\n"
            f"    - Not enough RAM (checkpoint is ~2.5 GB)\n"
            f"  How to fix:\n"
            f"    - Ensure at least 16 GB total RAM available\n"
            f"    - Close other applications to free memory"
        ) from e

    model.to(device)
    model.eval()

    # Load sequences
    logger.info(f"Loading sequences from {input_fasta}...")
    dataset = FastaBatchedDataset.from_file(input_fasta)
    alphabet = Alphabet.from_architecture("roberta_large")
    loader = DataLoader(
        dataset,
        collate_fn=alphabet.get_batch_converter(),
        batch_size=batch_size,
        num_workers=0,  # 0 for Windows/WSL compatibility
    )

    # Run inference
    logger.info(f"Running predictions on {len(dataset)} proteins...")
    all_names = []
    all_probs = []
    all_preds = []
    all_lengths = []

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(loader):
            toks = toks.to(device)
            out = model(strs, toks)
            prob = torch.softmax(out, dim=1)
            _, pred = torch.max(prob, 1)

            all_probs.append(prob.cpu().numpy())
            all_preds.append(pred.cpu().numpy())

            for i, s in enumerate(strs):
                name = labels[i].split()[0]
                all_names.append(name)
                all_lengths.append(len(s))

            if (batch_idx + 1) % 100 == 0:
                logger.info(f"  Processed {batch_idx + 1} batches...")

    all_probs = np.concatenate(all_probs)
    all_preds = np.concatenate(all_preds)

    # Write output CSV
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "deepsece_predictions.csv")

    n_secreted = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "protein_id", "deepsece_prediction", "deepsece_ss_type",
            "nonsec_prob", "T1_prob", "T2_prob", "T3_prob", "T4_prob", "T6_prob",
            "max_prob", "length",
        ])
        for i, pid in enumerate(all_names):
            pred_label = PREDICTED_LABELS[all_preds[i]]
            ss_type = SS_MAP[pred_label]
            max_prob = float(all_probs[i][all_preds[i]])
            if pred_label != "-":
                n_secreted += 1

            writer.writerow([
                pid, pred_label, ss_type,
                f"{all_probs[i][0]:.4f}", f"{all_probs[i][1]:.4f}",
                f"{all_probs[i][2]:.4f}", f"{all_probs[i][3]:.4f}",
                f"{all_probs[i][4]:.4f}", f"{all_probs[i][5]:.4f}",
                f"{max_prob:.4f}", all_lengths[i],
            ])

    logger.info(f"DeepSecE: {n_secreted}/{len(all_names)} predicted as SS substrates")
    return out_path


def parse_deepsece_output(results_path):
    """Parse DeepSecE output into standardized format."""
    entries = []

    for sep in [',', '\t']:
        try:
            with open(results_path) as f:
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    entry = {}
                    for raw_col, std_col in _COLUMN_MAP.items():
                        entry[std_col] = row.get(raw_col, '')
                    entries.append(entry)
                if entries:
                    return entries
        except Exception:
            continue

    return entries


def main():
    parser = argparse.ArgumentParser(description="Run DeepSecE")
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--checkpoint", default="", help="Path to checkpoint.pt")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        results_path = run_deepsece(
            args.input, tmpdir,
            checkpoint_path=args.checkpoint if args.checkpoint else None,
        )
        entries = parse_deepsece_output(results_path)

    logger.info(f"Parsed {len(entries)} DeepSecE predictions for {args.sample}")

    fieldnames = list(_COLUMN_MAP.values())
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for e in entries:
            writer.writerow(e)


if __name__ == '__main__':
    main()
