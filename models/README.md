# `models/`: trained model weights used by ssign

This directory is a **placeholder**. Model weights are too large to version in
Git, so they are hosted externally and fetched at install time.

## What goes here at runtime

| Model                | Size     | Used by                                   |
| -------------------- | -------- | ----------------------------------------- |
| DeepSecE checkpoint  | ~500 MB  | `run_deepsece.py`                         |
| ESM-1b               | ~7 GB    | DeepSecE feature extraction               |
| ESM-2                | ~3 GB    | PLM-Effector feature extraction           |
| ProtT5               | ~2.5 GB  | PLM-Effector feature extraction           |
| PLM-Effector weights | ~200 MB  | `run_plm_effector.py`                     |
| DeepLocPro weights   | ~2 GB    | `run_deeplocpro.py`                       |
| SignalP 6.0 weights  | variable | `run_signalp.py` (subject to DTU license) |

## Getting the weights

ssign ships a fetch script that pulls every weight file from its canonical
source:

```bash
bash scripts/fetch_weights.sh
```

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
`--signalp-mode remote --deeplocpro-mode remote`. See
[`docs/how-to/install.md`](../docs/how-to/install.md).

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
