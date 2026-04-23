"""ssign's vendored copy of PLM-Effector (Zheng et al. 2026).

Original source: https://github.com/zhengdd0422/PLM-Effector
Original author: Dandan Zheng (`zhengdanizh@ipbcams.ac.cn`)
Licence (for this directory only): Creative Commons Attribution 3.0. See
`LICENSE` in this directory for the full text.

Citation (required by CC-BY 3.0):
    Zheng, D. et al. "PLM-effector: unleashing the potential of protein
    language models for bacterial secreted protein prediction." Briefings
    in Bioinformatics 27(2), bbag143 (2026).

Why ssign vendors this code instead of `pip install`-ing it:
- There is no PyPI package for PLM-Effector upstream.
- Upstream install is conda-pinned to Python 3.9 + CUDA 11.3, which
  conflicts with ssign's Python 3.10+ baseline.
- Upstream `run_pipeline.py` contains hardcoded machine-specific paths.

Modifications made to the upstream code:
- Dropped `torch_geometric`, `scipy.sparse`, and `Bio.PDB` imports that
  were vestigial in the prediction path (they were only used in
  training-time utilities we do not need).
- Dropped training-only helpers (WeightedSampler, WeightedSampler_4combine,
  compute_class_weights, custom_collate).
- Refactored `run_pipeline.py` and the two subprocess-orchestrated scripts
  into a single callable `predict()` function exposed from this package.
- Replaced the file-based `.npz` intermediate protocol with in-memory
  tensor passing between feature extraction and ensemble stages.

Top-level API (shape, wired up in a follow-up commit — see 3.1.e.2):
    from plm_effector import predict
    predict(
        proteins_fasta="/path/to/input.faa",
        weights_dir="/path/to/plm_effector_weights",
        effector_type="T1SE",
        out_path="/path/to/output.tsv",
        device="cuda",
    )

This commit (3.1.e.1) lands the vendored skeleton only: the utils,
neural-net definitions, and inference helpers. The `predict()` entry
point, the refactored feature-extraction and ensemble stages, and the
ssign `run_plm_effector.py` wrapper are added in 3.1.e.2 once this
skeleton compiles cleanly.
"""

__all__: list[str] = []
