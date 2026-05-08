# `figures/`: publication figures and their regeneration scripts

This directory holds final figures used in the ssign manuscript and the scripts
that regenerate them from primary data. Anyone should be able to clone the
repo, run a single script per figure, and reproduce the exact image in the
paper.

## Convention

- One subdirectory per figure: `figures/fig1_pipeline_overview/`,
  `figures/fig3_benchmark/`, etc.
- Each subdirectory contains:
  - The final image(s): `fig1.pdf`, `fig1.png` (both formats).
  - `regenerate.py` (or `regenerate.sh`): self-contained script to produce
    the image from data.
  - `README.md`: one paragraph describing what the figure shows, inputs, and
    any non-obvious choices (colour scheme, scaling, etc.).

## Running

```bash
cd figures/fig1_pipeline_overview/
python regenerate.py
```

Scripts should read inputs from `../data/` or `../results/` (or fetch from
Zenodo); never hard-code paths to Teo's laptop.

## Status

Contents will be populated during Phase 8 of the publication roadmap. For now
this directory exists only as a placeholder to reserve the layout.
