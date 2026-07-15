# `models/`: trained model weights used by ssign

This directory is a **placeholder**. Model weights are too large to version in
Git, so they are auto-downloaded on first run (or baked into the container
image).

## What goes here at runtime

| Model                | Size     | Used by                                    |
| -------------------- | -------- | ------------------------------------------ |
| DeepSecE checkpoint  | ~2.5 GB  | `run_deepsece.py` (fine-tuned ESM-1b)      |
| ESM-1b               | ~7.3 GB  | DeepSecE backbone                          |
| ESM-2                | ~2.5 GB  | DeepLocPro backbone                        |
| ProtT5 (half)        | ~2.4 GB  | pLM-BLAST query embedding (extended tier)  |
| DeepLocPro weights   | ~2 GB    | `run_deeplocpro.py` (user-acquired, DTU)   |
| SignalP 6.0 weights  | ~1.5 GB  | `run_signalp.py` (user-acquired, DTU)      |

Total auto-downloaded on first run: ~14 GB (everything except the two
user-acquired DTU rows). Base tier needs all but ProtT5 (~12 GB); ProtT5
adds at the extended tier. Container users get every weight pre-baked in
the image, so nothing downloads at run time.

## Getting the weights

Each predictor fetches its own weights the first time it runs: DeepLocPro and
DeepSecE pull their ESM backbones through `esm.pretrained` (torch-hub cache),
pLM-BLAST pulls ProtT5 through `huggingface_hub`, and the DeepSecE checkpoint
comes down via `run_deepsece._ensure_checkpoint()`. No separate fetch step is
needed. The container image bakes all of them at build time (see
`containers/ssign.def`) so container runs are fully offline.

The **DeepSecE checkpoint** is currently hosted on an unreliable SJTU
server; mirroring to a Zenodo deposit before v1.0.0 release is part of the
longevity mitigation stack (see project plan). At release time the fetcher
flips Zenodo to the primary source with the upstream as fallback.

## DTU academic-licensed models

DTU confirmed on 2026-05-07 that SignalP 6.0 cannot be redistributed; users
obtain it directly from the [DTU portal](https://services.healthtech.dtu.dk/)
(free academic licence). DeepLocPro is pending separate clarification with
Ole, the DeepLocPro maintainer; treated as user-acquires-it for now.

ssign is offline-first, so the canonical path uses local DTU installs.
Users without a DTU licence can opt into the DTU webserver fallback with
`--signalp-mode remote --deeplocpro-mode remote`; this is a convenience
path for trial runs and licence-free use, and depends on DTU continuing
to host the service. For reproducible / paper-grade work, install
locally. See [`docs/how-to/install.md`](../docs/how-to/install.md).

## Integrity checking (post-v1.0.0)

Once the Zenodo mirror lands, every weight file will have an SHA-256
checksum recorded alongside it:

```bash
cd models && sha256sum -c checksums.sha256
```

For now, integrity relies on HTTPS + each upstream's own size validation.

## What must **not** go here

- User-trained or user-fine-tuned models (those are out-of-scope for ssign).
- Raw training data (that's in `data/` if anywhere).
- Model outputs (those are pipeline `results/`).
