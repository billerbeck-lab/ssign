## Why

The circular-shift enrichment test measures whether DeepLocPro and DeepSecE positive calls cluster at secretion-system loci, but it ignores SignalP. SignalP detects the N-terminal signal peptide that marks Sec-dependent export: autotransporter (T5aSS/T5cSS) passengers and T2SS/T5SS substrates carry one, while T1SS/T3SS/T6SS effectors bypass Sec. A SignalP enrichment track therefore adds a biologically meaningful, directly measurable signal (strong, expected enrichment where Sec is used; a clean negative contrast where it is not) that DLP/DSE alone do not capture.

## What Changes

- Add **SignalP** as a third circular-shift enrichment predictor alongside DeepLocPro and DeepSecE, for **all** secretion-system types (window types use the ±3-gene window mask; T5aSS/T5cSS use self-detection). Emits a SignalP fold + permutation p + BH q per type, in the same per-tool multiple-testing family as DLP/DSE.
- **Force whole-genome SignalP whenever `--enrichment-stats` is on**, the same way the flag already forces whole-genome DLP/DSE. The whole-genome SignalP machinery already exists (`Config.sp_whole_genome` + `runner._step_signalp`); it is simply not triggered by enrichment today. Runs **local** (this tool does not use webservers).
- `enrichment_testing.py` gains a `--signalp` input and a SignalP positivity vector; `run_enrichment` emits the SignalP track; SignalP joins `ENRICH_TOOLS`.
- The per-tool enrichment figure (`run_enrichment_figure.py`) gains a third SignalP bar per type (each tool an independent score).
- The combined (one-bar-per-type) figure shows DLP-or-DSE for window types and **SignalP-alone for autotransporters (T5aSS/T5cSS)** — not a 3-way OR. SignalP-alone has a lower genome background, so over the few autotransporter loci it gives a stronger, more often significant bar; a union would raise the background and weaken it (design decision 4).
- SignalP-positive is fixed as a Sec signal peptide (`VALID_SEC_SIGNAL_TYPES` = Sec/SPI + Sec/SPII), the definition `cross_validate` already uses.

## Capabilities

### New Capabilities
- `signalp-enrichment`: SignalP as a circular-shift enrichment predictor (whole-genome positivity, per-type fold/p/q, figure bar), forced on with `--enrichment-stats`.

### Modified Capabilities
<!-- No baseline spec captures the enrichment predictor set; the figures-v2 change's enrichment-stats delta is the closest, but it is a separate in-flight change. SignalP enrichment is introduced as a new capability rather than modifying a not-yet-archived delta. -->

## Impact

- **Code**: `src/ssign_app/scripts/enrichment_testing.py` (new `--signalp` input, positivity vector, SignalP track, `ENRICH_TOOLS`), `src/ssign_app/core/runner.py` (force `sp_whole_genome` under enrichment; pass the whole-genome SignalP TSV to the test; pool SignalP track), `src/ssign_app/scripts/run_enrichment_figure.py` (third bar), `src/ssign_app/scripts/ssign_lib/constants.py` (`ENRICH_TOOLS`), `plot_style.py` (SignalP colour).
- **Runtime**: one extra whole-genome SignalP pass per genome when enrichment is on (local; acceptable per Teo).
- **Docs/tests**: CLAUDE.md enrichment section, README; enrichment_testing + figure + runner pooling tests.
- **Out of scope**: the figures-v2 figure set; the combined-figure design (unchanged except the optional 3-way OR).
