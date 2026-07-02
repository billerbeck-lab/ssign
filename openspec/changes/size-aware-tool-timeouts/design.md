## Context

External tool calls are capped by a flat `TOOL_TIMEOUT_S = 14400` (4h) constant, applied at two layers:
1. **Runner**: `core/runner.py` calls `run_script("run_X.py", args, timeout=14400|7200)` at ~10 sites; `_run_script_streaming` kills the whole wrapper process on expiry.
2. **Wrapper**: several wrappers (`run_deeplocpro`, `run_signalp`, `run_interproscan`, `run_eggnog`, `run_blastp`, `run_bakta`, `run_ortholog_grouping`) additionally wrap the actual tool in `subprocess.run(cmd, timeout=TOOL_TIMEOUT_S)`.

The 52-genome benchmark (`MultiGenomeRunner`, `--enrichment-stats`) pools all genomes' proteomes into one 160,831-protein whole-genome prediction call. DeepLocPro/DeepSecE/SignalP each needed ≈5–5.5h; the runner's 4h `run_script` cap killed them (`"Timeout after 14400s"`), cross-validate got no output, and every genome emitted 0 substrates. It is a time wall, not memory (tools stream in batches).

ssign already ships a fitted per-tool runtime model: `runtime/effort_model.py` (`effort(tool, size, regime, coeffs) -> Effort(seconds, low_confidence, method, n, loo_pct)`, `resolve_regime(tool, whole_genome)`) backed by `runtime/coefficients.json` (linear `a·size + b` in reference-machine seconds, per tool/regime). It is currently consumed only by the GUI ETA (`Home.py`) and `calibration/replay.py`, not by the runner. This change reuses that model to size the timeout.

## Goals / Non-Goals

**Goals:**
- Replace flat per-tool timeouts with a size-aware cap derived from the existing effort model, for all 11 modelled tools.
- Keep the current 4h as a floor so small/neighborhood runs are bit-for-bit unchanged.
- Make the runner and wrapper timeout layers agree (one computed value threaded to both).
- Make a timeout an attributable, loud failure rather than a silent 0-substrate cascade.

**Non-Goals:**
- No change to `MultiGenomeRunner` pooling structure (per-genome / chunking rework is explicitly rejected).
- No PBS `--walltime` logic in ssign (set manually at CX3 submit).
- No new runtime dependency; no online machine-rate inference on the hot path (reference-machine effort + a margin is sufficient).
- DeepSecE `batch_size=1` → VRAM-auto batching is a separate speed optimization, not required for correctness.

## Decisions

**1. Derive the timeout from reference-machine effort × margin, not from a live machine rate.**
`scaled_timeout(tool, size, regime, *, margin=2.0, floor=TOOL_TIMEOUT_S)` returns `int(max(floor, ceil(margin × effort(...).seconds)))`, or `floor` when `effort` is None. Alternative considered: wire the online `Estimator` (which divides effort by an inferred per-machine rate) into the runner. Rejected for now: it adds state threading through the step loop, and at the point the first GPU tool runs there is no inferred GPU rate yet, so it would fall back to the reference prior anyway. A static margin on the reference estimate is simpler and, for the target hardware (RTX6000/A40 are within ~1.5× of the reference A40), safe. The online estimator remains available for a future tightening.

**2. Confidence-aware margin.** Well-calibrated fits (`low_confidence = False`; DLP/DSE whole_genome at 6.5–7% LOO error) use `margin = 2`. `low_confidence = True` fits use a wider margin so a thin/noisy fit cannot under-cap a genuinely-running tool. Alternative considered: a single flat margin. Rejected because several relevant fits are weak (SignalP whole_genome LOO 663% n=4, plm_effector whole_genome n=3, EggNOG substrates LOO 137%). The floor is the backstop for the degenerate case where the fit predicts far too low (SignalP whole_genome predicts ~13 min but the tool is genuinely fast, so the floor is never the killer).

**3. Compute in the runner, thread to the wrapper.** The runner is the enforcing layer (it killed the failed run) and it knows the tool + input FASTA at each call site, so it computes `scaled_timeout(tool, count_sequences(input), regime)` and passes it to `run_script(..., timeout=...)`. It also passes `--timeout N` to the wrapper so the wrapper's inner `subprocess.run` uses the same value (default `TOOL_TIMEOUT_S` when standalone). Alternative considered: each wrapper computes its own timeout independently. Rejected: duplicates the effort-model read and the size/regime resolution in every wrapper, and risks the two layers drifting. Single owner (runner), value threaded down.

**4. `TOOL_TIMEOUT_S` becomes the floor, not the cap.** Its numeric value (14400) is unchanged; only its role changes. `plm_blast` keeps its existing 24h floor (`PLMBLAST_TIMEOUT_S`). Per-protein caps (HH-suite `HHBLITS/HHSEARCH_TIMEOUT_S`) are a different structure and are left as-is.

**5. Fail-loud.** When the runner observes a step timeout, it emits a message naming the tool and the cap before the downstream step fails on the missing output, so the true cause is visible in the log rather than only the cascade.

## Risks / Trade-offs

- **A confidently-wrong low fit under-caps a slow tool** → Mitigated by the confidence-aware margin (wider for `low_confidence`) and the floor. The only fits that both exceed the floor and are noisy are EggNOG (substrates) and possibly IPS; these get the wider margin, and validation measures their real time to refit.
- **Reference machine ≠ run machine** → The 2× margin absorbs the RTX6000/A40 gap; a much slower GPU could still exceed 2× of the reference estimate, but that is out of the current hardware envelope and would be caught by validation timing.
- **`coefficients.json` now read on the hot path** → It is a small bundled JSON already loaded elsewhere; read once per run and cache. A missing/corrupt file must degrade gracefully to `floor` (treat as unmodelled), not crash the run.
- **EggNOG at panel scale** → predicted ≈4.4h at 819 pooled substrates already brushes the 4h floor; the scaled cap (wider low-confidence margin) covers it, and validation confirms.

## Migration Plan

Additive and backward-compatible: small runs keep the 4h floor. Rollout is the code change plus a validation run of the 52-genome single job with `--walltime 24:00:00`. Rollback = revert to the hardcoded timeouts (no data migration). The validation run also records real vs predicted per-tool wall-clock to feed the effort-model refit.

## Open Questions

- Exact low-confidence margin: fixed (e.g. 3× or 4×) vs derived from `loo_pct`. Default to a fixed wider value unless the LOO-scaled form is clearly better; settle during implementation with the measured EggNOG/IPS times.
- Whether to also thread the timeout into the `t5ass_whole` pooled annotation pass (same tools, substrate-scoped) or accept the floor there for now.
