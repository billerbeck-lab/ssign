## 1. Timeout helper

- [x] 1.1 Add `scaled_timeout(tool, size, regime, *, margin=2.0, floor=TOOL_TIMEOUT_S) -> int` (new `src/ssign_app/runtime/timeouts.py`), returning `int(max(floor, ceil(margin * effort(tool,size,regime,coeffs).seconds)))`, and `floor` when `effort` is None. Load `coefficients.json` once (cache); degrade to `floor` if the file is missing/corrupt.
- [x] 1.2 Confidence-aware margin: apply a wider margin when the effort fit has `low_confidence=True` (fixed wider value; leave the LOO-scaled variant as an open question resolved in 5.x).
- [x] 1.3 Unit tests: floor honored, margin applied, `ceil`/`int` return type, unmodelled tool → floor, low-confidence widening, degenerate low fit still ≥ floor (SignalP whole_genome case).

## 2. Runner wiring

- [x] 2.1 Added `PipelineRunner._scaled_tool_timeout(tool, size, whole_genome, floor)` + `_substrate_count()`; wired the 8 sizeable modelled sites: DLP/DSE/SignalP/plm_effector (size=`count_sequences(input_proteins)`), IPS/EggNOG/plm_blast/protparam (size=`_substrate_count()`). **Left flat (documented):** bakta + extract_proteins (size = output proteome, not known pre-run; per-genome bounded), macsyfinder (fast, per-genome, no `run_script` timeout), blastp/hh_suite/orthologs-blast (unmodelled → helper would floor anyway).
- [x] 2.2 `regime` via `resolve_regime(tool, whole_genome=<config flag>)`; `size` from the FASTA/substrate set the runner passes (pool runner sees the pooled 819 substrates for segment D).
- [x] 2.3 `plm_blast` passes `floor=PLMBLAST_TIMEOUT_S` (24h).

## 3. Wrapper timeout plumbing

- [x] 3.1 Added `--timeout` (default `TOOL_TIMEOUT_S`, or `PLMBLAST_TIMEOUT_S` for pLM-BLAST) threaded to the inner `subprocess.run` in the 5 wired-and-modelled wrappers: `run_deeplocpro`, `run_signalp`, `run_interproscan`, `run_eggnog`, `run_plm_blast` (both embed + search sites). Also updated the hardcoded "timed out after 4 hours" messages to `{timeout}s`. blastp/bakta/orthologs wrappers left untouched (unmodelled / left flat).
- [x] 3.2 Runner appends `--timeout <scaled>` for DLP/SignalP/IPS/EggNOG/plm_blast so inner and outer caps agree.
- [x] 3.3 DeepSecE + plm_effector + protparam have no inner `subprocess.run` timeout → runner-side `run_script` timeout only; confirmed those scale via the runner change.

## 4. Fail-loud

- [x] 4.1 Both `run_script` timeout handlers (non-streaming + streaming) now `logger.warning` naming the script + cap ("hit its Ns timeout and was killed — a hard cap, not a tool crash") and return `Timeout after {timeout}s (killed {script})`.

## 5. Validation

- [~] 5.1 Run 3165431 (2026-07-02): 51/52 genomes emitted substrates, but SignalP hit the 4h floor (its `a=0` fit) and was killed → T5SS undercounted. DLP/DSE + all annotation succeeded (fix validated). **A clean rerun is still needed** after the 5.4 SignalP patch.
- [x] 5.2 Real timings recorded in `calibration_observations.md`: DLP 5.56h (pred 5.43h), DSE 5.17h (4.97h), pLM-BLAST 7.6h (8.83h), IPS 81m (95m), EggNOG 22m (4.47h) — all on the pooled 160k/819 sets. DLP/DSE fits validated to <4%.
- [x] 5.3 EggNOG pooled-substrate real = 22min (819 subs); fit predicted 4.47h → **12× over** (safe/loose, fixed-cost-dominated small-N fit). Flagged for official refit, not blocking.
- [~] 5.4 SignalP `whole_genome` (`method=mean`, a=0, killed the run) **patched to the measured DLP rate** (a=0.1214, `method=proxy_deeplocpro`, low_confidence) → 16.3h cap; unit test updated. EggNOG + plm_effector fits recorded for the official `calibration/fit.py` refit (separate repo). SignalP's true rate lands from the rerun.

## 6. Docs

- [x] 6.1 CLAUDE.md Key Parameters: `TOOL_TIMEOUT_S` is now the floor; caps scale as `max(floor, 2× effort prediction)`, low-confidence wider; noted which tools stay flat.
- [x] 6.2 Added a header note to `scripts/cx3/submit_batched_overnight.sh` (near the `WALLTIME` default) that a full panel needs `--walltime 24:00:00`.

## 7. Stretch (separable, not blocking)

- [ ] 7.1 DeepSecE `batch_size=1` → VRAM-auto batching by reusing `auto_batch_size_from_vram` (the slowest predictor at batch=1); guard memory and keep batch=1 as the CPU/low-VRAM fallback.
