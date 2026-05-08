# `colab/`: Google Colab demo notebook

This directory will hold the **"try ssign in your browser"** demo notebook
shipped alongside the paper: a zero-install, single-genome trial that lets
anyone with a Google account test ssign without setting up their own machine.

## Planned contents

- `ssign_demo.ipynb`: end-to-end notebook:
  1. Clones the ssign repo into the Colab session.
  2. `pip install -e .` in the session's Python environment.
  3. Downloads model weights + a minimal database subset from Zenodo.
  4. Uploads (or fetches) a user-supplied genome.
  5. Runs the pipeline on the free T4 GPU.
  6. Displays the resulting HTML report inline.

- `README.md`: this file (scope, limitations, badge code).

## Limitations

Colab is great for a demo, **not** a production service:

- Session timeout ~12 h on the free tier; 24 h on Pro.
- GPU is unpredictable (T4 usually; occasionally no GPU).
- Storage is ephemeral; every session re-downloads weights + database.
- 390 GB BLAST NR is out of scope; demo uses a small curated reference set.
- Single-genome only; multi-genome batches need local install.
- Subject to Google's ToS; no sustained public web-service hosting.

## Why ship this at all

Lowest-friction way to let reviewers, educators, and curious researchers
actually run ssign on their own input without installing anything. Also
useful as a reproducibility artifact attached to the paper.

## Status

Notebook will be created in Phase 8 of the publication roadmap, once the
pipeline is stable and model weights are on Zenodo. For now this directory
reserves the layout.
