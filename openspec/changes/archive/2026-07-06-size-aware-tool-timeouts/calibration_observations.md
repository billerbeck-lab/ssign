# Calibration observations — pooled run 3165431 (2026-07-02, RTX6000)

Real per-step wall-clock from the 52-genome single-job run, on the pooled sets
(**160,831 proteins** for whole-genome predictions; **819 substrates** for
annotation). Feed these into the official `calibration/fit.py` refit. Reference
machine for `coefficients.json` is CX3-A40; this run was **RTX6000**, so a
per-machine rate factor applies (RTX6000 ≈ A40 here — DLP/DSE landed within 4%).

| tool | regime | size | actual_s | predicted_s (A40 fit) | pred/actual | verdict |
|---|---|---|---|---|---|---|
| deeplocpro | whole_genome | 160831 | 20024 | 19558 | 0.98 | fit excellent — keep |
| deepsece | whole_genome | 160831 | 18600 | 17886 | 0.96 | fit excellent — keep |
| signalp | whole_genome | 160831 | **36076** (run 3169556) | 758 (a=0, mean) → 19558 (DLP proxy) | 0.54 | **measured**: real a ≈ 0.224 s/protein (~1.8× DLP); proxy under-predicted but 3×-margin cap (16.3h) held |
| plm_blast | substrates | 819 | 27342 | 31784 | 1.16 | conservative — keep (also has 24h floor) |
| interproscan | substrates | 819 | 4860 | 5687 | 1.17 | good — keep |
| eggnog | substrates | 819 | 1333 | 16096 | **12.08** | over-predicts 12× (fixed-cost-dominated small-N fit); safe (loose cap) but refit needed. real ≈ 1.63 s/sub |
| protparam | substrates | 819 | 0.4 | ~0 | — | fine |

## Actions taken
- **signalp.whole_genome**: manually proxied to `deeplocpro.whole_genome` (a=0.1214,
  b=28.9), `method="proxy_deeplocpro"`, `low_confidence=True`. Cap went 4h → 16.3h.
  SignalP was killed incomplete so its true rate is unknown; DLP (heavier transformer,
  same 160k input) took 5.56h, so SignalP is expected ~5.5-6h. The rerun will measure it.

## For the official refit (not done here — calibration repo is separate)
- **eggnog.substrates**: add the 819→1333s point; the current slope (a=19.26) is
  inflated by small-N (≤~30 substrate) calibration runs where the DIAMOND DB-load
  fixed cost dominates. Real marginal rate ≈ 1.6 s/sub.
- **signalp.whole_genome**: replace the DLP proxy with the rerun's measured SignalP time.
- **plm_effector.whole_genome**: still `a=0, n=3` (degenerate), but PLM-Effector is off
  by default (`skip_plm_effector=True`) so it never runs whole-genome in practice; low
  priority. It remains the unit-test example of a degenerate-fit-protected-by-floor.
