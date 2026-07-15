# Deferred work

Tracks items skipped during tasks. One bullet per item: what, why, trigger to revisit.

## ⭐ RESUME 2026-07-13 — "container just works" plan (approved); Phase 1 DONE

Plan file: `~/.claude/plans/expressive-cooking-pebble.md` (make extended/full
container self-contained + durable; base benefits too — one image, one runner).
Driven by a long CX3 debugging chain: every extended-container failure was
host-tool integration, not ssign logic (missing tRNAscan → wrong bakta shadow →
stale apt blastn 2.12 → Bakta ENOSPC on the 64 MiB `--writable-tmpfs` /tmp).

- **Phase 1 (scratch fix) = DONE + pushed (172f121).** `runner.resolve_scratch_dir()`
  points `$TMPDIR`+tempfile at real disk (keep $TMPDIR if ≥5 GB free, else
  outdir/.ssign_scratch) before any work_dir; wired into `run()` AND
  `_run_multi()` (multi path bypasses run()). New `--scratch-dir`. 1456 tests green.
- **Phase 2 (DEF DONE + PUSHED 2026-07-13; CX3 validation pending) = bake free toolchain.**
  Commits: `5547759` (def bake + tasks 5.4) + `fa9d523` (stripped container PBS),
  pushed on `enrichment-circular-shift-per-run`. Image built + smoke-verified:
  `~/.ssign_sif_build/ssign_v3.sif` 5.7 GB, bakta 1.12.0 / blastn 2.17 / diamond /
  real hmmsearch / tRNAscan from `/opt/conda/envs/bakta`, emapper.py from eggnog env.
  **REMAINING = the CX3 parity proof** (self-contained image, NO host bakta-deps/eggnog):
    1. laptop: `scp ~/.ssign_sif_build/ssign_v3.sif ttr25@login.cx3.hpc.ic.ac.uk:ephemeral/ssign_v3.sif`
       (fallback if no `~/ephemeral` symlink: ssh in, `echo $EPHEMERAL`, scp to that abs path)
    2. CX3: `cd ~/blastp_t5a/ssign && git pull`
    3. CX3: `qsub -v GENOME=$HOME/xantho_gbff/<substrate-yielding>.gbff,SIF=$EPHEMERAL/ssign_v3.sif ~/blastp_t5a/ssign/scripts/cx3/run_container_extended.pbs`
       then health-check `qstat -f <id> | grep -E 'job_state|Resource_List.select'` (gpu_type must be present).
    Pass = 23/23 steps with EggNOG/IPS/pLM-BLAST columns populated, matching a native run.
    SUBMITTED 2026-07-13: job **3282805.pbs-7** on Roseixanthobacter_finlandensis_VTT_E-85241
    (image confirmed at /rds/general/user/ttr25/ephemeral/ssign_v3.sif). select has
    gpu_type=RTX6000 (healthy; queued waiting on busy RTX6000, not a config bug).
    Check later: `qstat -u ttr25` (R = running); when `ssign.run.log` says "Pipeline
    complete", retrieve via cx3-retrieve skill. Watch emapper-vs-bakta-diamond compat.
    Watch: emapper picks up bakta env's diamond 2.2.3 on PATH (not eggnog's) — confirm
    compatible with emapperdb-5.0.2 when the real DB is mounted (smoke couldn't test it).
- (superseded detail) Phase 2 def AUTHORED: micromamba added to `%post`; `/opt/conda/envs/bakta`
  Def AUTHORED: micromamba added to `%post`; `/opt/conda/envs/bakta`
  (`-c conda-forge -c bioconda bakta=1.12.0` → blast 2.17 + real HMMER + amrfinder +
  diamond + tRNAscan/aragorn/infernal/pilercr) + isolated `/opt/conda/envs/eggnog`
  (eggnog-mapper=2.1.13, biopython 1.76 pin isolated); solved EARLY in %post
  (fail-fast before torch); dropped apt `ncbi-blast+`; `%environment` prepends both
  env bins; env specs frozen to `/opt/conda/*.lock.yml`. **Scope tightened vs plan:**
  (a) **IPS NOT baked** — it ships software+member-DBs as one 24 GB dir that must be
  host-mounted for the DBs anyway, so that mount already brings `interproscan.sh`
  (baking = redundant). (b) **ESM bake DEFERRED** — DeepSecE's ESM1b alone is 7.3 GB
  (not 2.5), PLM-E pulls 3 more, both off by default; default-on DeepLocPro is
  host-mounted + pulls its own ESM. Needs checkpoint verification (DTU pkg) before
  baking the right one; host torch-cache mount already covers it. So image ~7 GB,
  not 12-18. **BUILDING NOW:** bg PID 660926 → `~/.ssign_sif_build/build_v3.log`,
  out `ssign_v3.sif` (v2.sif = last known-good, untouched). Monitor task b5xq5m36p
  waits for exit. After green build: smoke (bakta 1.12.0 / blastn 2.17 / emapper.py
  all from /opt/conda), then commit def + tasks.md, then CX3 extended re-validate
  with ONLY DBs+DTU mounted (no host bakta-deps/eggnog). Tasks #18-22.
- **Phase 3 = `scripts/ssign-run` one-line wrapper: DRAFTED + dry-run verified 2026-07-13.**
  Portable launcher (no PBS): auto-detects DB sub-trees under `--db-root`, binds DTU
  envs, scratch→/tmp, torch cache, auto `--nv`. Assembled apptainer cmd matches the
  validated PBS. Remaining: README doc + refactor PBS to call it (after 5.3 passes).
- **DTU-easy = DONE 2026-07-13 (task 5.6).** Local DTU is now the hard default: runner
  errors (not silent webserver) when no local SignalP/DeepLocPro; `--<tool>-mode remote`
  is explicit opt-in. New `scripts/ssign-setup-dtu` installs both into `~/.conda/envs/*`
  (DeepLocPro fully auto from GitHub; SignalP needs the user's licensed tarball only).
  1457 tests green. Facts: memory `reference_dtu_install`. KEY finding: DeepLocPro is
  NOT licence-gated (public GitHub) — only SignalP is. Easiest-path = image + fetch-DBs
  + ssign-setup-dtu + ssign-run.
- **"Just works" glue DONE 2026-07-13 (task 5.8 + ssign-run polish).** `ssign fetch-databases`
  subcommand (runs the bundled fetch script from the image, no host tools); README
  "Container quickstart" 4-step block; `ssign-run` auto-detects the host SignalP env
  (scans conda roots) so it needs ZERO DTU flags after ssign-setup-dtu. 1460 tests green.
  **Only release-time item left for full "just works": GHCR + Zenodo publish** (needs a
  token, user-run; I can stage the `apptainer push` / upload commands on request).
- **v4 image = DeepLocPro BAKED (task 5.7, 2026-07-13).** `~/.ssign_sif_build/ssign_v4.sif`
  9.7 GB: DeepLocPro (pin 32c0b37c, --no-deps into ssign env) + its ESM2-650M (2.5 GB)
  baked. **Verified in-image offline**: real prediction runs under torch 2.12 (testprot→
  Outer Membrane 0.82). So SignalP 6 is now the ONLY host-provided predictor. Image is
  NON-COMMERCIAL (CC BY-NC-SA DeepLocPro + AGPL EggNOG; Teo OK'd, academic tool).
  **CX3 v4 validation supersedes the queued v3 job 3282805** (kill it): scp ssign_v4.sif
  to CX3 ephemeral, `git pull`, qsub run_container_extended.pbs (now defaults to v4, mounts
  only SignalP + DBs). ssign-run/PBS updated: torch-cache mount opt-in (`--torch-cache`),
  no host DeepLocPro mount.
- Stays host-provided forever: DTU tools (SignalP/DeepLocPro) + reference DBs.
- The mount-based `run_container_extended.pbs` (gpu_type in directive, PREPEND
  /usr/local/bin:bakta-deps:eggnog, /tmp real-disk bind) is a validation HARNESS,
  superseded once Phase 2 bakes the tools. Last container run got to Bakta
  ENOSPC (job 3278618) → Phase 1 fixes that; needs a rerun to confirm.

## ⭐ RESUME AFTER COMPACT (2026-07-08)

**State:** full-tier (tier-3) is code-complete + reviewed + pushed on branch `enrichment-circular-shift-per-run`
(HEAD e62bc60). HH-suite dir→prefix bug fixed (98b2f5f); BLASTp switched NR→Swiss-Prot (74829d6) + review
cleanup (e62bc60). Full unit suite 1434 green. Xanthobacter 74-genome extended panel = DONE + clean +
inspected (T1/T4/T5a/T5b/T5c/T6, zero T3SS = real Xanthobacter biology; result in ~/Desktop/cx3_runs).
**IN FLIGHT (Teo-driven CX3, node login-ai):**
1. **Full-tier DB install = COMPLETE** (2026-07-08). Swiss-Prot + Bakta-full both fetched; `ssign doctor
   --tier full` = 8/8 databases OK (Bakta `bakta/db` full 84G, Swiss-Prot, UniRef30 all green), 2/2 weights,
   31/31 py. db-light still parked at `$EPHEMERAL/ssign-databases/bakta_db_light_backup` — reclaim with
   `rm -rf` once smoke-test-2 passes. (doctor's only 2 FAILs when run by hand = hhsearch/hhblits not on the
   login shell's PATH; cosmetic — real runs self-supply it via the fix below.)
2. **HH-suite PATH fix (NEW, this commit).** `hhblits` lives in `~/.conda/envs/hhsuite/bin` but `~/.bashrc`
   does NOT export it, so PBS runs never had it on PATH (why smoke-test-1's HH-suite failed before the prefix
   fix could even run). `run_batched_multi.pbs` + `run_k12_validation.pbs` now auto-detect + PATH-add it.
3. **Full-tier smoke-test-2 = DONE + PASS for HH-suite/BLASTp** (run 3252429, 2 Roseixanthobacter, HEAD
   aa1c936, in `~/Desktop/cx3_runs/batched_RTX6000_20260708_194046_3252429`). HH-suite → OK (PATH fix +
   prefix fix both validated), BLASTp → OK (Swiss-Prot), 16/16 both genomes, 24 substrates. BUT **EggNOG
   FAILED** = task #12: diamond exit 1 / empty stderr = OOM. `run_eggnog.py` sizes `--dbmem` (~44GB) +
   `--block_size 8` (~48GB) ≈ 92GB from full-node `effective_ram_gb()`, but runs in a 6-tool parallel group
   → overcommits the 120GB node → diamond OOM-killed. Nondeterministic (passed smoke-1 by luck). Optional so
   run still hit 16/16, but a PANEL hits it every genome / could OOM a core tool. Fix: budget EggNOG RAM ÷
   parallel-group concurrency (same pattern resolve_threads uses for threads). GATES the full-tier panel.
4. **#11 + #12 both FIXED (2026-07-09, committed + 4-agent simplify-reviewed, 1444 unit tests green).**
   #12 EggNOG OOM (commit 1e20f66): new `parallel_share_ram_gb()` in resources.py (RAM analogue of
   `parallel_share_cpus`, via shared `_parallel_group_size()`); EggNOG's 3 autodetect fns budget --dbmem/
   --block_size/--index_chunks off the group RAM share → on 120GB/6-tool group share=20GB → --dbmem off,
   no overcommit. #11 Bakta `db*` glob (next commit): fetch_databases.sh keys skip-guard on exact variant
   subdir (db/ vs db-light/) so extended→full re-downloads; dependency_manifest resolve_path returns
   `dirname(max(matches))` so db/ (full) beats db-light/. **VALIDATED: full-tier smoke-test-3 PASSED
   2026-07-09** (run 3256040, 2 xantho genomes): EggNOG-mapper -> OK (the fix worked, no OOM), HH-suite +
   BLASTp OK, 16/16 both genomes. All 4 full-tier fixes confirmed end-to-end. **NEXT: retrieve smoke-test-3
   (confirm emapper --block_size 2 + pull HH-suite/BLASTp strings for the consensus tool-weighting design);
   reclaim bakta_db_light_backup; run a real full-tier panel; close full-tier-cx3-wiring openspec (validation
   tasks 5.1-5.3 now satisfiable).**
**OPEN:** #8 (un-gate DeepSecE for T3SS? — recommend keep gate) undecided. db-light still parked at
bakta_db_light_backup (reclaim after smoke-test-3 confirms the resolve fix uses db/).
**NOTE:** CX3 has 822GB NR unused-by-default (opt-in via --blastp-db, or delete to reclaim).
Next likely: rerun smoke test → if EggNOG green, a real full-tier panel; close full-tier-cx3-wiring
openspec change + signalp-t5ss (task #4).

## 2026-07-12 — full-tier reruns retrieved (CLEAN); container extended = Bakta gap fixed; HH-suite timeout OPEN

Retrieved both full-tier reruns (in `~/Desktop/cx3_runs/`): run 3261357 = **74 genomes** (Xantho+Roseixantho, 17/17 each), run 3261459 = **57 genomes** (benchmark, 17/17 each). Pooled enrichment is biologically clean: Sec-dependent types high on SignalP (T5aSS self DLP 28×/SignalP 7.6×***, T5bSS DLP 37×, T2SS all sig), Sec-independent a clean SignalP negative (T1SS 0.5×, T3SS 0.39×, T6SS 0.2×, all n.s.) with DLP/DSE still high; T5cSS rescued by SignalP in benchmark (4.8×, q<0.001). Figures render on house style.

**OPEN ITEMS:**
- **HH-suite pooled step timed out at the flat 4h floor (14400s) in BOTH reruns.** Optional homology annotation only; substrate calls + enrichment unaffected, but the HH-suite annotation column is empty. Not in the size-aware timeout scaling. Trigger to revisit: if HH-suite annotations are wanted in the paper → add hhsuite to `runtime/effort_model.py` scaling or raise its floor; else accept (EggNOG/InterProScan/consensus cover function).
- **Container extended validation (reproducible-install 5.3):** run 3268527 failed at Bakta (missing tRNAscan-SE + external binaries, not baked). Fixed `run_container_extended.pbs` to mount host `bakta-deps` (commit 0099ca9); RERUN PENDING. Archival image should bake Bakta's toolchain + ESM ("mount now, bake for archival").
- **v2 functional figures (03–06) NOT yet verified** — weren't in the retrieval tarball (filter grabbed enrichment only). Re-pull `figures/` to confirm the v2 consensus rendered.

## 2026-07-09 — consensus-annotation-v2 LANDED (tool-weighted voter); CX3 full-tier reruns NEXT

**State:** consensus-annotation-v2 openspec change APPLIED (15/15 tasks). `annotation_consensus.py` +
`functional_vocab.py` rewritten to the tool-weighted voter (port of `docs/development/v2lib.py`). Gate
(`docs/development/regrade.py`) PASSED: port-fidelity 672/672, Xantho top-cat 39% (<42%), fallback-noise
0 (was 37), 55-gold outside-vocab 13→3, full-tier fold-name false calls fixed. 1449 unit tests green.
Golden fixture regenerated (surgical, GBFF-only). Change is ready for `/opsx:archive` after the reruns confirm.

**SUBMITTED 2026-07-09 (Teo-driven CX3): BOTH genome sets running at full tier** so the paper's functional
figures use the v2 consensus. Two pooled jobs, `--tier full --enrichment-stats --walltime 72:00:00`: Xantho 74
(`~/xantho_gbff/*.gbff`) + benchmark 75 (`$HOME/blastp_t5a/ssign/validation_sweeps/benchmark/inputs_gb/*.gbff`).
**RESUME = retrieve via cx3-retrieve → check v2 functional figures (03–06) + `_enrichment_stats` → then
`/opsx:archive consensus-annotation-v2`.** Runtime caveat: benchmark 75 at full tier is multi-day
(per-substrate HH-suite over UniRef30) and **may exceed 72h**; if it times out, decide split vs raise-walltime
vs per-genome enrichment. (Xantho full-tier ≈ extended for the *consensus* figures: BLASTp/Swiss-Prot is empty
on environmental Roseixanthobacter, so full only adds down-weighted HH-suite.)

Deferred from the v2 apply (revisit triggers noted):
- **Numeric-suffix machinery miss.** `_MACH` uses strict word boundaries, so `Hcp1`/`VgrG1`/`Tae2` (digit
  suffix) don't match → those components fall to `Unclassified` instead of `Apparatus-associated`. Faithful to
  the validated prototype; not a mislabel (honest Unclassified), just a lost Apparatus signal. Revisit if the
  figures show Apparatus under-counting, or fold into the curated-override follow-on.
- **Garbage-in residuals** (Tde1→Phage, VasX→Oxidoreductase, apxIIIA→Protease): tools carry no correct call.
  Explicitly out of scope here → the curated effector-family override follow-on change.
- **TOOL_WEIGHT + tier cutoffs live inline in `annotation_consensus.py`**, not `ssign_lib/constants.py`
  (the stated thresholds SSOT). Kept inline: no other consumer, and moving adds an import layer. Revisit if a
  future webserver UI needs them user-tunable (see memory feedback_user_controllable_thresholds).

## 2026-07-08 — Full-tier smoke-test bugs FIXED (HH-suite prefix + BLASTp→Swiss-Prot)

Both full-tier smoke-test bugs fixed (full suite 1434 green). **#9 HH-suite:** `_resolve_ffindex_prefix`
in run_hhsuite.py converts the resolved DB dir → ffindex prefix (`.../UniRef30_2023_02`) so hhblits gets
a prefix not a dir (was IsADirectoryError in the SSD stager); +5 tests. Commit 98b2f5f. **#10 BLASTp:**
Teo chose Swiss-Prot over NR. Full-tier BLASTp now defaults to **Swiss-Prot** (curated, ~300 MB, minutes)
not NR (~800 GB, >2h for 24 substrates → days on a panel). Changes: `fetch_blast_swissprot` in
fetch_databases.sh (run_full fetches swissprot, NOT nr — nr fn kept for manual opt-in); manifest NR entry
→ Swiss-Prot (SSIGN_BLAST_SWISSPROT / blast_swissprot / swissprot.pdb) so doctor checks the right DB;
runner resolves blastp_db from SSIGN_BLAST_SWISSPROT → `<dir>/swissprot` (NR no longer auto-resolved, opt-in
via --blastp-db); PBS exports SSIGN_BLAST_SWISSPROT; timeout 7200→TOOL_TIMEOUT_S (4h). Docs synced across
env_vars, run_on_hpc, install, README, design_decisions (+ new §3.5), cx3-submit skill. Tests swapped
NR→Swiss-Prot. **Committed 74829d6 + e62bc60** (2-agent simplify review passed: `-db swissprot` and the
ffindex-prefix derivation both verified correct upstream; review-flagged stale NR strings in cli/Home/
constants/run_blastp/pbs/preflight all fixed in e62bc60). **NOTE:** Teo's CX3 install has the
822GB NR on disk from the earlier fetch; it's now unused-by-default (harmless; usable via --blastp-db, or
delete to reclaim). A fresh --tier full fetch would pull Swiss-Prot (~300MB) + skip nr. **RE-VALIDATE:**
re-run the full-tier smoke test after `git pull` to confirm both fixes (HH-suite completes, Swiss-Prot
BLASTp completes in minutes) before a real full-tier panel.

## 2026-07-08 — Retrieved + inspected: Xantho panel (CLEAN) + full-tier smoke test (2 bugs)

Both runs retrieved to `~/Desktop/cx3_runs/`. **Xanthobacter panel (extended tier, run 3232160): CLEAN
+ good science.** 74/74 genomes at 17/17, no failures, pooled enrichment = 32 cross-genome rows, ~593
substrates. Pooled enrichment is biologically correct: T1SS strongly enriched (DLP 8.9× q3e-4, DSE 5.5×,
SignalP 0.5× NOT — Sec-bypassing C-terminal signal ✓), T5aSS-self very strong (DLP 27.7×, SignalP 7.6×
✓ Sec-dependent), T5bSS strong (DLP 36.9×, COMBINED 2.4×), T6SS sig via DSE (5.0× q0.017, SignalP 0.16×
NOT ✓ injected). **NO T3SS row + no T2SS** → Xanthobacter (environmental) has no injectisome, confirms
the "0-T3SS fits Xanthobacter not Xanthomonas" suspicion from the §3.3 genus flag. This is a valid,
publication-relevant result. (Real per-tool wallclocks for calibration/runs.jsonl still TODO from the log.)
**Full-tier smoke test (run 3238415, 2 Roseixanthobacter genomes, ~24 substrates): did its job, caught 2
real bugs.** My NR + UniRef30 RESOLUTION wiring WORKED (both DBs found, tools started against them). Both
new tools then failed (optional → run still completed 17/17 without them): (1) **HH-suite crashed** =
task #9 (dir→ffindex-prefix bug, resources.py:453). (2) **BLASTp timed out** at the hardcoded 7200s =
task #10 (24 substrates >2h vs 822GB NR; scales badly — real panel could be a day+). IPS/EggNOG/pLM-BLAST/
physicochem all OK at full tier. So full-tier plumbing is ~validated; the 2 heavy homology tools need
fixes (#9 code, #10 timeout + a design call on NR-vs-smaller-DB) before a real full-tier panel.

## 2026-07-07 — OPEN DECISION: un-gate DeepSecE for T3SS? (task #8, Teo undecided)

Teo floated: when Flagellum is purposefully included in a run, also allow DeepSecE to call
flagellum/T3SS (or even un-gate in general). Undecided; parked as task #8 (blocked on the Xantho
retrieve for real numbers). Mechanism (confirmed): DeepSecE has NO flagellum class
(`run_deepsece.py` PREDICTED_LABELS = non-secreted/T1/T2/T3/T4/T6SE); flagellum is a T3SS homolog,
so flagellar proteins funnel into T3SE. Current gate = `_dse_flag` unconditional
(`cross_validate_predictions.py:149`). **Claude's recommendation: keep the gate, don't couple it** -
because (1) DeepSecE can't emit "flagellum", un-gating only re-admits undifferentiated T3SS labels;
(2) including Flagellum ALREADY finds flagellar substrates via DLP-extracellular + proximity
(DeepSecE adds only mislabeling); (3) genomes with flagellum + a real injectisome (common) get their
real T3SS list contaminated; (4) coupling to Flagellum-inclusion is arbitrary (contamination is
T3SS-proximity-driven); (5) T3SS already has 2 reliable signals (MacSyFinder + DLP 71% ext on real
effectors). Honest alt if we want DeepSecE's flagellar signal: RELABEL a DeepSecE-T3SE call inside a
flagellum window as flagellar-associated (deliberate feature, arguably out of ssign's remit).
Added the "no flagellum class" mechanism to design_decisions §3.3 (commit 15a9a48). *Trigger:* after
Xantho run retrieved, recount DeepSecE-T3SS vs flagellar annotation, then Teo decides. **Also
outstanding:** genus label in §3.3 ("Xanthomonas") likely wrong - 0-injectisome fits Xanthobacter.

## 2026-07-07 — T3SS/flagellum doc drift fixed; one spec discrepancy left for Teo

Verified (for Teo): ssign DOES distinguish T3SS injectisome from flagellum. TXSScan 1.1.4 ships
separate `Flagellum.xml` + `T3SS.xml` models (T3SS built from injectisome `sct*` incl. the OM secretin
`sctC` the flagellum lacks); ssign reads the type from `model_fqn`. Flagellum is in the default excluded
set (dropped at substrate-calling, `system_filtering.py:276`); T3SS is NOT excluded (detected +
substrate-called). DeepSecE can't tell them apart, so it's flagged unconditionally for T3SS
(`t3ss-detection` spec). Fixed stale "T3SS excluded by default" claims across 7 docs (commit 0e3da7a):
design_decisions §2.1+§3.3, pipeline_overview, configure, troubleshooting, cli.md, output_files, README;
also corrected §3.3's wrong "flagged only if MacSyFinder validates T3SS" (it's unconditional) and the
default excluded VALUE (real set = Flagellum Tad T4aP T4bP MSH ComM Archaeal-T4P, docs had shown 2-3).
**DEFERRED (spec-vs-code, Teo's call):** `openspec/specs/t3ss-detection/spec.md:10` still says the default
`excluded_systems` SHALL be `[Flagellum, Tad]`, but the constant + the authoritative `secretion-system-scope`
spec have all 7 appendages. The `t3ss-detection` req over-specifies the full list (its real concern is only
"SHALL NOT contain T3SS"). Fix = reword line 10 to defer the full list to secretion-system-scope. Didn't
edit unprompted (synced contract spec). *Trigger:* next time the t3ss-detection spec is touched, or Teo OKs.

## 2026-07-06 — fetch_databases.sh pre-flight (fail-fast on missing tools)

Added `_preflight_tier <tier>` (checks every external tool the tier needs BEFORE any download,
instead of lazily when the fetch reaches that DB). Motivating case: `update_blastdb.pl` (BLAST NR)
was checked only at the LAST step, so a missing NR tool wasted the hours spent fetching Bakta +
HH-suite first. Hoisted the 5 install hints to `_HINT_*` vars so pre-flight + per-fetch checks can't
drift; kept the per-fetch `_require_command` calls as defense-in-depth. Honors `--dry-run`
(informational). Tests: `tests/unit/test_fetch_databases_preflight.py` (fail-fast, extended doesn't
require NR tool, full checks NR up front). Root cause it prevents = the live 2026-07-06 false start
(amrfinder/update_blastdb.pl not on PATH). **Deferred (out of scope, pre-existing):** `--dry-run` of
any bakta-fetching tier exits 1 — `fetch_bakta`'s `bakta_subdir=$(compgen -G "$dir/db*" | head -1)`
fails under pipefail when the DB dir isn't created in a preview. Confirmed present in HEAD (not mine).
Fix would be `|| true` on that assignment, but dry-run may have further preview-only warts downstream,
so a clean `--dry-run` exit is its own task. *Trigger:* if anyone relies on `--dry-run` exit code.

## 2026-07-06 — TIER-3 (full) run setup: 3 wrapper/runner gaps found (planning, no code yet)

Tier-3 (`--tier full`) adds `blastp` (vs NR) + `hhsuite` (vs UniRef30) on the substrate
set only. Auditing the CX3 wrapper + runner against `--tier full` surfaced 3 gaps that
would silently fail a full-tier run if we only changed the tier flag:
1. **`run_batched_multi.pbs:89` hardcodes `--tier extended`.** SSIGN_EXTRA_ARGS is appended
   AFTER it, so `SSIGN_EXTRA_ARGS="--tier full"` works (argparse last-wins), but that path
   is fragile (spaces in the PBS `-v` value). Clean fix: a `--tier` flag on
   `submit_batched_overnight.sh` → `TIER` env → `--tier $TIER` in the PBS (default extended).
2. **PBS wrapper never exports `SSIGN_HHSUITE_UNICLUST` (UniRef30).** It exports pfam+pdb70
   only (lines 72-73). The runner DOES read `SSIGN_HHSUITE_UNICLUST` (runner.py:403,432), so
   HH-suite just needs the export added: `[ -d "$DBROOT/hhsuite/uniref30" ] && export
   SSIGN_HHSUITE_UNICLUST="$DBROOT/hhsuite/uniref30"`.
3. **BLASTp NR does NOT auto-resolve — real inconsistency.** The manifest defines
   `SSIGN_BLAST_NR` (subpath `blast_nr`, sentinel `nr.pdb`) but the runner's config
   resolution (Steps A+B, runner.py:399-448) omits `blastp_db`, so exporting `SSIGN_BLAST_NR`
   does nothing and the run dies at runner.py:2099 "BLASTp requires a local database". Two
   fixes: (a) quick — append `--blastp-db $DBROOT/blast_nr/nr` in the wrapper; (b) proper —
   add `("blastp_db","SSIGN_BLAST_NR")` to the runner's Step-A tuple so NR behaves like every
   other DB (then wrapper just exports `SSIGN_BLAST_NR`). Recommend (b); it's a genuine gap.
**Runtime caveat (load-bearing):** blastp + hhsuite are NOT in the effort model, so
size-aware timeouts don't scale them — fixed caps only (BLASTp = TOOL_TIMEOUT_S 4h floor;
HH-suite = HHBLITS_TIMEOUT_S 1h + HHSEARCH_TIMEOUT_S 30min PER PROTEIN). HH-suite runs
per-substrate via ThreadPoolExecutor (cpu_per_job=2 → ~32 parallel on 64c). Benchmark pooled
substrate set was ~819; a 74-genome Xantho panel likely ~1000-1200 substrates → HH-suite is
the tail, plan --walltime 24h+ (realistic hhblits ~2-10 min/protein vs UniRef30). DB fetch
(`fetch_databases.sh --tier full`) is ~700 GB total: NR ~390 GB (needs `update_blastdb.pl`
from ncbi-blast+), UniRef30 ~25 GB (GWDG only), Bakta full upgrade. *Trigger:* when Teo OKs
the tier-3 run — propose the wrapper+runner change via /opsx:propose, fetch DBs on $EPHEMERAL,
verify with `ssign doctor --tier full`, then submit.

**IMPLEMENTED 2026-07-06 (openspec `full-tier-cx3-wiring`, 13/16 tasks, 1426 unit tests green).** Runner
auto-resolves NR from `SSIGN_BLAST_NR` via the shared `find_db_by_env_var` helper (dir→`/nr` prefix);
`submit_batched_overnight.sh --tier <t>` + PBS `TIER` var + UniRef30/NR exports. 4-agent simplify review
clean; only fix applied = use `find_db_by_env_var` at both runner lookup sites instead of inline `next()`.
Third-party review confirmed `<dir>/nr` is the correct BLAST+ `-db` value (v5 NR `nr.pal`/`nr.pdb` alias).
Remaining tasks 5.1-5.3 = CX3-user-driven validation (fetch + smoke test + record wallclock).
**CX3 DB STATE (checked 2026-07-06):** extended DBs all present (bakta db-light, eggnog, IPS, ECOD30,
plm_effector, taxdump); ALL 4 full-tier DBs MISSING (NR + hhsuite pfam/pdb70/uniref30). `--tier full`
fetch pulls ~555 GB NEW incl an unwanted Bakta-full swap (db-light works; doctor accepts it).
Deferred review items (out of scope, Teo's call): (1) sibling PBS scripts `run_serial_genomes.pbs` +
`run_k12_validation.pbs` still hardcode `--tier extended` — no tier flag; only the batched path is
full-tier-wired. (2) `fetch_databases.sh` `fetch_blast_nr` doesn't pass `--blastdb_version 5` to
update_blastdb.pl (relies on the modern default; fails loud on a legacy v4 NR, never silent-wrong).
**FETCH STARTED 2026-07-06 (PID relaunched after PATH fix).** Root cause of the first two false starts:
`amrfinder` + `update_blastdb.pl` (both in `~/.conda/envs/bakta-deps/bin`) were NOT on PATH — the
`export PATH=~/.conda/...` before `source .venv/bin/activate` didn't survive the activate. Fix:
`export PATH="$HOME/.conda/envs/bakta-deps/bin:$PATH"` AFTER activation, verify with `command -v`, then
launch. **BAKTA-FULL SKIPPED:** `fetch_bakta full`'s _already_done guard keys on the `bakta/` parent,
sees existing `db-light`, and skips the ~84 GB full `db` download. So this fetch delivers full-tier
NR + HH-suite (pfam/pdb70/uniref30) but Bakta stays db-light (good side effect: no db/db-light resolver
ambiguity). **Teo CONFIRMED he wants true Bakta-full (2026-07-06)** — deferred until AFTER the Xantho run
finishes (can't move db-light while it's in use): stage `bakta/db-light` aside (mv, don't delete), re-run
the fetch to pull `bakta/db`, verify `ssign doctor`, then the full run uses the full DB. Once `db` (full)
exists alongside `db-light`, the resolver's `db*/version.json` sentinel is ambiguous — so remove/rename
db-light after, leaving only `db`. NR is the ~390 GB tail (fetched LAST, after HH-suite). PATH fix that
unblocked the fetch: export `bakta-deps/bin` AFTER `source .venv/bin/activate` (activate resets PATH);
cx3-submit skill Databases prelude corrected to match + given a `command -v` pre-launch gate.

## 2026-07-02 — t5-hitchhiker-enrichment LANDED (openspec, 21/21, 1422 unit tests green)

T5aSS/T5cSS autotransporters now emit TWO enrichment results (`mode` column): `self`
(autotransporter self-detection, COMBINED=DLP-or-SignalP) + `window` **hitchhiker**
(secreted-predicted ±3 neighbours that may piggyback the T5 pore, COMBINED=DLP-or-DSE).
Enrichment figure draws them as two adjacent `(self)`/`(hitch)` x-groups; annotation
figures split T5a/c by `substrate_source`. Deleted fig-02 (autotransporter scatter),
renumbered curated set `01`–`06`, removed `fig_autotransporter` toggle, removed fig-01
single-genome bar labels. Hitchhiker concept documented (design_decisions § 5.2).
Null key gained `mode` (`enrich_null_key(ss, tool, mode)`) to stop self/window collision.
Simplify pass (4 agents): promoted shared `constants.enrich_combined_uses_signalp` +
`T5_SELF_TAG`/`T5_HITCH_TAG`/`ENRICH_MODE_*`/`SUBSTRATE_SOURCE_T5_SELF` constants (killed
a stat-vs-figure drift trap); correctness review found no bugs.
Committed + pushed `3490de3` (partial-staged runner.py to exclude the parallel
size-aware-timeouts work). **CX3-VALIDATED 2026-07-03** (5-genome pooled job 3165435 +
single E. coli 3165443): pooled T5aSS emits self (DLP 7.0×✱✱, COMBINED 5.1×✱✱✱) AND
hitchhiker window (DLP 3.1×, COMBINED 2.0×, DLP-or-DSE confirmed — SignalP-window obs 13
ignored); figures render the two adjacent `(self)`/`(hitch)` groups (self=SignalP hue,
hitch=DLP-or-DSE hue), single-genome path works (E. coli T5aSS self DLP 7.3×✱, hitch
4.4× faded), annotation `01` shows the split with muted-hitch colour. Hitchhiker signal
is real (E. coli's 7 T5aSS neighbours). **ONLY REMAINING:** archival order — archive
`figures-v2` + `signalp-enrichment-track` BEFORE this change (its run-figures/enrichment
deltas supersede theirs). *Trigger:* finish + archive the two predecessors, then archive
`t5-hitchhiker-enrichment`.

## ⭐⭐ 2026-07-06 — BENCHMARK DONE: run 3169556 clean, recall 38/85 (44.7%).

**FINAL RESULT (run 3169556, 2026-07-03, final code):** 52/52 genomes 17/17; SignalP succeeded (real 36076s=10.0h
on 160k, ~1.8× DLP; 16.3h proxy cap held). Recall **38/85 = 44.7%** (T1 16/18, T2 1/8, T3 5/19, T4 0/8, T5 10/19,
T6 6/13), **0 flips** vs the prior baseline → the ca18322 T5SS gate is recall-neutral on gold (drops only FP T5).
All 85 molecules reconciled. Recall recomputed via `RerunIndex.emitted_overlap`; recall figures regenerated
(scripts 27/32/52 → `data/phase2/figures/{,summary/}`); pooled auto-figures in the tarball. Gallery artifact:
https://claude.ai/code/artifact/9109aac0-04c5-4215-91ab-f4778a189e07 . New rerun reshaped into
`rerun/<unit>/results/`; old 2026-06-25 preserved as `rerun.stale_20260625/`.
**⚠️ FULLASM CAVEAT:** the 4 full-assembly T1SS genomes (NZ_CP031766.1/CBDBTK/JABJZG/JBCGCZ) are still served from
the OLD `rerun_fullasm/` (2026-06-25); the override is LOAD-BEARING (+3: Serralysin/apxIA/hlyA only recoverable
from the full assembly — fragment inputs cut the gene). T1SS is unaffected by recent changes so 38/85 is very
likely final, but for single-code-version rigor re-run those 4 on new code from `inputs_gb_fullasm/` (small, ~1-2h).
**NEXT / status (2026-07-06):** (b) `size-aware-tool-timeouts` ARCHIVED (specs synced; DeepSecE VRAM-auto batch
added, task 7.1, commit f353593); (c) benchmark docs DONE = `validation_sweeps/benchmark/docs/phase2_recall_results.md`
(commit 6b5df1a; tasks 4.4/4.5 closed). (a) fullasm 4-genome rerun command handed to Teo (uses `inputs_gb_fullasm/`
on CX3, --walltime 12h) — pending his run. **Xanthobacter panel run — SUBMITTED + RUNNING (2026-07-06).**
74 Bakta `.gbff` from `/home/teo/Desktop/Billerbeck - SS Identification/all_Xantho_incl_SpringerLab/Annotations/bakta/gbff/`,
rsync'd to CX3 `~/xantho_gbff/`, single full tier-2 MultiGenome job, `--enrichment-stats --walltime 72:00:00`, default
Bakta re-annotation (no `--use-input-annotations`). Estimate ~40h wall (~400k pooled proteins; SignalP tail ~25h real vs
~40h proxy cap). Run dir on CX3: `~/runs/batched_RTX6000_<ts>_<jobid>/`. *Next:* retrieve + inspect when it lands
(~2026-07-08) — confirm 74/74 at 17/17, record real per-tool wallclocks for the effort-model calibration.

## (history) 2026-07-03 — run 3165431 = 16/17 (SignalP-only fail); SignalP coeff PATCHED; RERUN READY.

**Outcome of run 3165431 (2026-07-02):** timeout fix VALIDATED — DLP 5.56h (pred 5.43h), DSE 5.17h (4.97h),
EggNOG/IPS/pLM-BLAST all OK on the pooled 160k/819 sets. **Only SignalP failed:** its `whole_genome` fit was
`method=mean` (a=0, n=4, size<=5572) → floored at 4h → killed at 14400s. Every genome 16/17 (SignalP the only
miss); substrates DID emit (DLP succeeded) so non-T5 recall is valid but T5SS is undercounted + SignalP enrichment
empty → **not usable as the final benchmark; must rerun.** Real timings + fit validation saved in
`openspec/changes/size-aware-tool-timeouts/calibration_observations.md`. **FIX applied + committed:** SignalP
`whole_genome` proxied to DLP rate (a=0.1214, method=proxy_deeplocpro, low_confidence) → scaled cap 4h→16.3h;
SignalP true time still unmeasured (killed) but concurrent DLP=5.56h implies ~6h → 16.3h is safe. EggNOG 12×-over
+ plm_effector a=0 recorded for the official calibration refit (safe/loose, non-blocking). Local run copy +
extract: `validation_sweeps/benchmark/bench_rerun2.tgz` + scratchpad `rerun2/`.
**RERUN COMMAND (after git pull on CX3):** `bash scripts/cx3/submit_batched_overnight.sh --enrichment-stats
--walltime 72:00:00 $(sed "s#^#$HOME/blastp_t5a/ssign/validation_sweeps/benchmark/inputs_gb/#" validation_sweeps/benchmark/data/phase2/rerun_inputs.txt)`.
Expect ~15h wall (segment B ~6h + segment D pLM-BLAST tail 7.6h + A ~1.7h). Then reshape + recompute recall.

## (history) 2026-07-02 — CX3 VALIDATION RUN HANDED OFF (single-job, timeout-fixed).

Committed+pushed `4a510bd` (size-aware-tool-timeouts) on top of figures `3490de3`, branch
`enrichment-circular-shift-per-run`. Handed Teo the **single-job** CX3 submit (the path that failed before,
so it validates the fix AND gives gold-list recall in one go):
`bash scripts/cx3/submit_batched_overnight.sh --enrichment-stats --walltime 72:00:00 $(sed "s#^#$HOME/blastp_t5a/ssign/validation_sweeps/benchmark/inputs_gb/#" validation_sweeps/benchmark/data/phase2/rerun_inputs.txt)`
after `cd ~/blastp_t5a/ssign && git pull`. Config verified = gold-list baseline: DLP/DSE/SignalP + IPS/EggNOG/
pLM-BLAST/ProtParam ON, PLM-E/BLASTp/HH OFF, Bakta reannotate (use_input_annotations=False). Lands in
`~/runs/batched_RTX6000_<ts>_<jobid>/<genome>/<genome>_results{,_raw}.csv` (FLAT, no results/ subdir).
**WHEN IT LANDS:** (1) tar the per-genome `_results{,_raw}.csv` + enrichment + ssign.run.log; scp down;
(2) **RESHAPE** into `validation_sweeps/benchmark/rerun/<unit>/results/` (the layout RerunIndex reads —
insert the `results/` level); (3) recompute recall under ca18322 T5SS gate (expect T5SS 10/19 -> lower);
(4) group-5 validation: confirm `[_pool]` predictions OK (not `Timeout after 14400s`) + record real-vs-predicted
per-tool wallclock. Health markers: `[_pool] ... -> OK`, each genome `N/N steps succeeded`.

## ⚠️ 2026-07-02 — single-job 52-genome rerun FAILED (pooled-prediction timeout); redo via HARNESS

The ca18322 benchmark rerun submitted via `submit_batched_overnight.sh` (single MultiGenomeRunner job,
`batched_RTX6000_20260701_154401_3156011`) is a **total loss: 0 secreted proteins emitted for all 52
genomes**. Cause: MultiGenomeRunner **pools every genome's whole proteome into one prediction call** —
here a 53 MB / **160,831-protein** FASTA (`_pool/pooled_whole_genome.faa`) — and `--enrichment-stats`
forces whole-genome DLP+DSE+SignalP on it. All three hit `TOOL_TIMEOUT_S=14400` (4h) and FAILED, so every
genome's cross-validate died with "No DeepLocPro output". The ~5h wall was 3 tools grinding to the timeout.
- **Fix / redo:** use the per-genome benchmark **harness** `scripts/cx3/SUBMIT_rerun.sh` (12 batch jobs,
  each genome a separate ssign run, whole-genome preds ~13 min/genome, no pooling). Same full-standards
  config (`REANNOTATE=1,INCLUDE_T3SS=1,WHOLE_GENOME=1,ENRICH=1,ANNOT=1`). This is the path that produced
  the working 2026-06-25 results and lands in `rerun/<unit>/results/` (the layout RerunIndex reads). Also
  submit `SUBMIT_rerun_fullasm.sh` for the 4 full-assembly T1SS genomes (needs `inputs_gb_fullasm/` on CX3).
- **DEFERRED (real bug, design decision pending Teo):** `MultiGenomeRunner` segment B (predictions
  DLP/DSE/SignalP/PLM-E) pools ALL genomes into ONE call; with `--enrichment-stats` that pool is the whole
  proteome, so at 52 genomes it's 160k proteins in one `subprocess.run(timeout=4h)` per tool → all three
  killed at 4h. Verified specifics:
    * Scope: only the **segment-B predictors in whole-genome mode** hit this (DLP, DSE, SignalP; PLM-E
      structurally but off fleet-wide). **Segment-D annotation** (EggNOG/IPS/pLM-BLAST/BLASTp/HH) pools over
      **substrates only** (small), so it's safe at panel scale (pLM-BLAST already has a 24h cap).
    * It's a **TIME wall, not memory**: all three stream in batches (DLP withholds mega-proteins via
      `partition_by_length`; SignalP `--mode fast`; DSE `DataLoader`). Raising the timeout WOULD let it finish
      (no OOM) but gives a monolithic multi-hour call with no checkpointing.
    * **Efficiency find:** DSE runs `batch_size=1` (run_deepsece.py:212), never overridden — slowest of the
      three. PLM-E already has `auto_batch_size_from_vram` + `--chunk-size`; DSE should reuse that.
  - **Recommended fix = MINIMAL (revised 2026-07-02 after measuring throughput; Teo pushed back on the big
    rework, correctly).** The 4h cap was just ~10–30% too short, not a structural pooling flaw. Measured on the
    successful harness K-12 run (4,314 proteins, same 3 tools, same 1-GPU sharing): SignalP 1.5min, DSE 8.4min,
    DLP 9min (concurrent → ~9min group wall). Scaled to the 160,831-protein pool (one amortized model load):
    parallel group ≈ **~4.5–6h**, DLP-dominated. The tools write output only at the end, so the failed run's
    "0 output at 4h" = ~90% done, not stuck. **Fix: make the timeout size-aware** in the 3 prediction wrappers
    — `timeout = max(4h, n_proteins/rate × margin)` (wrappers already call `count_sequences`); ~400 prot/min ×2
    → 160k ≈ 13h cap (real work fits, a hung DTU call still dies). Scale, don't delete (delete loses the
    hung-tool guardrail; fine for local-only but strictly worse). **Also raise PBS walltime**: wrapper default is
    8h, a full panel needs ~8–12h → re-run `submit_batched_overnight.sh --walltime 24:00:00` (48h first time).
    Optional speedup (not required): DSE `batch_size=1` → reuse PLM-E `auto_batch_size_from_vram`.
  - **REJECTED (over-scoped):** per-genome segment-B rework / chunking `multi_runner`. Not needed — only segment
    B (whole-genome, 160k) hits the cap; segment D pools just **819** substrates (measured), nowhere near it.
  - **IMPLEMENTED 2026-07-02 (groups 1-4 + 6; `/opsx:apply size-aware-tool-timeouts`).** New `runtime/timeouts.py`
    `scaled_timeout(tool,size,regime)=max(floor, margin×effort)` (margin 2, low-confidence 3); runner
    `_scaled_tool_timeout`+`_substrate_count` wired at 8 sites (DLP/DSE/SignalP/plm_effector + IPS/EggNOG/plm_blast/
    protparam); `--timeout` threaded to 5 wrappers; both `run_script` timeout handlers fail loud. bakta/macsyfinder/
    extract_proteins stay flat (size unknown pre-run / never pooled). Verified: DLP 160k→10.9h, EggNOG 819→13.4h,
    small runs→4h floor; 199 dedicated tests green (scaled_timeout 8 + 5 wrapper suites + runtime 191). **Deferred:**
    group 5 (CX3 validation re-run — needs the user) + group 7 (DSE VRAM batch, stretch). **⚠️ COMMIT NOTE:**
    `core/runner.py` is shared — it also carries the parallel session's uncommitted `fig_*` config lines, so a
    path-scoped commit of runner.py would sweep those in; coordinate before committing (let figures session commit first).
  - **PROPOSAL WRITTEN 2026-07-02:** openspec `size-aware-tool-timeouts` (4/4 artifacts, validates). Generalized
    per Teo to ALL 11 modelled tools (not just whole-genome predictors): `scaled_timeout = max(4h floor,
    2×effort-model prediction)`, wider margin for low-confidence fits, computed in the runner + threaded to
    wrappers via `--timeout`. Walltime stays manual (out of scope). Next: `/opsx:apply` to implement.
  *Decision pending:* implement the fix now (then re-run single job with `--walltime 24:00:00`) vs harness-now-fix-later.

## 2026-06-30 — figures-v2 (parallel figures session, separate from the gold-list work below)

Done + tested this session (full unit suite green, 1392): enrichment figure → single combined
DLP+DSE fold/significance bar chart (`_enrichment_fold.png` / `pooled_enrichment_fold.png`);
per-run restructure (T5a/c pulled out of cargo figs `01`/`03`/`04`, new `05` autotransporter
self-detection scatter); functional-category candidates `06`–`13` (COG/KEGG/EggNOG/curated-consensus
× overall+by-SS-type) in `generate_figures.py` + new `ssign_lib/functional_vocab.py`; genome-group
`0N_pooled_*` figures; `outer_membrane_prob` plumbed end-to-end to the integrated CSV.

LANDED / DEFERRED items from figures-v2:
- **Appendage-exclusion default — LANDED 2026-06-30 (Teo approved).** `DEFAULT_EXCLUDED_SYSTEMS` now
  `[Flagellum, Tad, T4aP, T4bP, MSH, ComM, Archaeal-T4P]`, wired via the constant across
  cli/Home/runner/system_filtering/validate_macsyfinder_systems. **HEADS-UP for the gold-list/recall
  session:** default substrate scope shrank — proteins whose only neighbour is a T4 pilus / MSH /
  ComM / archaeal pilus are no longer substrate-called. They were never real secreted effectors, so
  recall on the gold list should be unaffected (it removes false positives), but if you recompute
  found/total against a fresh run, expect fewer total substrate calls.
- **Consensus-function keyword grouping — AUDITED 2026-07-09** (was REVISIT flag from 2026-06-30).
  Full audit in `docs/development/consensus_annotation_audit.md` (graded against Xantho 593 + benchmark
  55-gold accuracy sheet + lit). **Grade C-.** Key defects: 3 over-broad rules (Autotransporter/
  Structural/Transporter) win 60% of calls; `t[1-9]ss` in the "Secretion system" rule hijacks effectors
  → Apparatus (Tae amidases, many T3SS effectors); classify returns ALL matches → 809 vote-splits;
  Hypothetical + title-case fallback are winnable (28/593 noise); `fli[a-z]` matches "flippase"/"conflict";
  13/55 gold land outside vocab (RTX toxins→Protease, TPS adhesins→Nuclease, S-layer→Protease).
  **Proposed v2 in the report** (one-vote-per-tool + specificity rank; Hypothetical/Other as floors;
  drop fallback minting; machinery-by-component-identity not keyword; lit-grounded 2-level taxonomy).
  Grading script: `docs/development/grade_consensus.py`. **NEXT (gated on Teo review): implement v2 (task #15).**
- **SignalP-as-enrichment-predictor — LANDED 2026-06-30, VALIDATED 2026-07-01 (openspec `signalp-enrichment-track`, 15/16).**
  SignalP is now a third circular-shift predictor for all SS types: per-tool figure gets a 3rd bar;
  combined figure uses **DLP-or-SignalP for ALL T5 (a/b/c)** (reversed from the initial SignalP-alone
  decision, per Teo), DLP-or-DSE for window types, DLP-only for T3SS.
  `--enrichment-stats` forces whole-genome SignalP (local). SignalP-positive = `VALID_SEC_SIGNAL_TYPES`
  via the new shared `constants.is_sec_signal_peptide` (also adopted by t5ss_handler). Unit + integration
  tests green (1410). **Task 5.1 CX3 validation DONE** (single PAO1 + 5-genome pooled, 2026-07-01,
  `~/Desktop/cx3_runs/batched_RTX6000_20260701_1602*`): SignalP row on every type; T3SS has no DSE row;
  T5 combined = union re-test (T5aSS 5.1×***, T5bSS 2.3×**, T5cSS 5.7× rescued by SignalP alone, DLP=DSE=0);
  faded-not-grey n.s. bars; titles say "significance"; `n_null` = 5679 (single, exact rotations) vs 10000
  (pooled MC). **Ready to `/opsx:archive signalp-enrichment-track`.**
- **T5SS substrate gate LANDED 2026-07-01 (openspec `signalp-t5ss-substrate-call`, code+tests done).**
  A trace found the T5 DLP rule (bug-fix #6) was **inert**: `t5ss_handler` emits every detected T5
  component and `system_filtering` kept them all unconditionally; cross_validate's `is_secreted` is
  never read for inclusion. Added a real gate in `system_filtering._t5_self_has_evidence`: a T5 component
  is a substrate iff DeepLocPro localizes it (T5a/c ext-or-OM; T5b OM-only, preserving #6) OR SignalP
  predicts a Sec signal (Sec-only; Tat excluded, literature-confirmed). Net: T5 counts DROP (evidence-less
  components removed). 1413 unit tests green. **CX3 rerun 2026-07-01 confirmed the gate WORKS and DROPS:**
  enrichment self-mask detected 18 T5aSS + 3 T5cSS components; only 12 + 2 retained as substrates, so
  the gate **dropped 7 autotransporters (6 T5a + 1 T5c)** that had neither DLP OM/ext ≥0.8 nor a Sec
  signal. (Correction: an earlier note in this session said "0 dropped" — that was wrong; it was based
  on "0 leaks among *retained* rows", which can't see dropped rows.) All retained T5SS-self rows carry
  evidence, 0 leaks; none of the drops leaked back via proximity (the 7 T5aSS-proximity rows are
  disjoint, strongly-extracellular *neighbors*, not the dropped components). **DEFERRED:**
  (1) regenerate the `t5ass_minimal` golden fixture (its T5 substrate `dlp_extra=0.60`<0.8,
  `signalp=OTHER` is now correctly dropped), so the frozen results changed (integration test skips
  locally w/o tools; regen per REGENERATE.md on a full install). (2) 52-genome benchmark for the drop
  *magnitude* per genome (the 5-genome run already shows the drop exists; the stalled 52-genome CX3
  batch would quantify it at scale). *Trigger:* full-install run for the golden regen; 52-genome CX3 batch.
- **Figure-02 SignalP fill uses a looser rule than the enrichment track (REVISIT, surfaced by simplify
  2026-06-30).** `generate_figures._signalp_positive` is a denylist (`_SIGNALP_NEGATIVE`) that counts
  TAT/TATLIPO/PILIN as "has signal", while the enrichment SignalP track + t5ss_handler use the strict
  allowlist SP/LIPO only (`is_sec_signal_peptide`). So a TAT-signal autotransporter shows filled in
  figure 02 but counts negative in the enrichment stat. *Why deferred:* changing it shifts figures-v2
  output (separate change); needs Teo's call on whether figure-02 fill means "any signal peptide" or
  "Sec signal peptide". *Trigger:* decide alongside the figures-v2 figure-02 review; if aligning, point
  `_signalp_positive` at `is_sec_signal_peptide`.
- **Functional-candidate pruning (figures-v2 task 6.3).** `06`–`13` are over-produced on purpose. Teo
  reviews a real run, names the keep-set; delete the rejected figure functions in a follow-up.
- **Pooled figures ignore per-figure toggles (minor).** `multi_runner.py`'s `generate_figures --mode pooled`
  call doesn't forward the user's `--no-*` flags, so `0N_pooled_*` always renders the full curated set
  even if a per-genome figure was disabled. *Why deferred:* arguably-correct (the group view should be
  complete) and low-value; forwarding needs the config.fig_*→--no-* mapping runner.py already has.
  *Trigger:* only if a user reports the inconsistency.

4-agent simplify pass on the full figures-v2 diff (2026-06-30) caught a shipped crasher: the enrichment
figure used `ordered_ss_types` without importing it (only the opt-in integration test covered the script,
so the default unit run was green). Fixed + added unit-level render coverage
(`test_run_enrichment_figure.py`). Also fixed: `functional_vocab` import fallback, a triplicated
missing-value check (now `functional_vocab.is_missing`), and a return-code check in `pool_and_plot_enrichment`.

## ⭐ RESUME HERE (2026-06-30) — gold v6 (85 rows, recall 38/85 = 44.7%); EXHAUSTIVE BLIND RE-READ DONE, list FINAL; figs in a SEPARATE session

**EXHAUSTIVE BLIND RE-READ DONE 2026-06-30 (Teo: "do the exhaustive blind reread").** 6 blind agents (round-robin
batched so each spans every SS type; inputs via `scripts/76_make_blind_batches.py`) independently re-verified all 87
rows from primary evidence -- per row: (A) the UniProt IS the named effector (not a machinery/accessory paralog),
(B) the locus is that gene, (C) the cited paper EXPERIMENTALLY shows secretion by this ss_type. **Result: 83/87
confirmed; ZERO wrong-gene anchors, ZERO fabricated/mismatched accessions, ZERO length errors** -- the 4 earlier
wrong-gene catches (Tae4/Tle1/YspE/Ssp2) were the last of that class. The ONLY substantive finding: 4 T3SS rows are
TRANSLOCON/apparatus components (secreted BY the T3SS but not delivered host-effector cargo). Teo's scope call =
drop pure-machinery, keep dual translocator/effectors:
- **DROPPED** T3SS_28 EspA (needle-filament sheath, PF03433) + T3SS_16 HrpK1 (translocon-pore helper, PF16937) --
  both correctly anchored, both found=no, pure structural translocon -> out of scope for an effector gold list.
- **KEPT** T3SS_01 CopN (YopN/LcrE gatekeeper w/ documented tubulin/cell-cycle effector activity) + T3SS_21 BipC
  (SipC/IpaC-family translocator w/ documented actin-nucleation effector activity) -- genuine dual-role effectors.
Also swapped 2 weak citations to verified PRIMARY papers (blind agent confirmed DOI + verbatim quote by fetching):
T4SS_02 Ceg14 -> Burstein 2009 PLoS Pathog (CyaA translocation, lpg0437 in Table 3); T2SS_03 LapA -> Ball 2002 Mol
Microbiol (Hxc-system discovery, hxc mutant loses LapA from supernatant). Cross-strain ortholog cites (T1/T3
Tobe-2006-table rows, T5SS_03 iga) left as-is -- identity + mechanism are correct, only the paper is an ortholog's.
**gold v6 = 85 rows, recall 38/85 = 44.7%** (T1 16/18, T2 1/8, T3 5/19, T4 0/8, T5 10/19, T6 6/13). Run order
65->73->65 (idempotent); blind verdicts saved in scratchpad `blind_verdict_{1..6}.tsv`. **THE GOLD LIST IS NOW FINAL.**

### (history) 2026-06-30 — gold v5 (87 rows); machinery-anchor recheck

**MACHINERY-ANCHOR RECHECK DONE 2026-06-30 (Teo: "recheck the gold list for placeholder/machinery anchoring").**
`scripts/75_machinery_anchor_audit.py` checks, per row, what the protein at the row's locus is ANNOTATED as
(RefSeq gene name + Bakta rerun product) and flags apparatus/accessory/structural/regulatory classes + near-
sequential locus clusters (the Ysa placeholder pattern). **Result: 0/87 machinery-anchored, 0 placeholder
clusters.** Detector self-tested to fire on all 4 known-bad products (YscL/TagJ/TagF/OmpA) and pass real effectors.
**Dropped T3SS_30 NleD** (Teo): real verified effector, correctly anchored, but Bakta missed its ORF (IS-dense, the
NleD Bakta-gene-miss) so it had no product and was found=no for an ANNOTATION reason -> dropped so recall measures
prediction+proximity on cleanly-annotated genes. **Recall now 38/87 = 44%.**

(That re-read is now DONE -- see the ⭐ block above. This subsection is kept only for the machinery-anchor history.)

## (history) 2026-06-30 — gold (88 rows); DOMAIN/NAME audit (dropped YspE+Ssp2 machinery anchors)

**DOMAIN/NAME AUDIT DONE 2026-06-30 (Teo: "whenever I ask you find something wrong — was there nothing else to check?").**
Built `scripts/74_domain_audit.py`: for every row, fetch the UniProt entry and flag (a) RED_DOMAIN = a machinery/
structural/mobile-element Pfam (the systematic version of the T6SS_17/18 catch), (b) NAME_MISMATCH = the UniProt
entry doesn't name the effector, (c) length/organism mismatch. Plus integrity (no dup uniprot/locus/coords) and a
blank-field sweep. KEY structural lesson Teo flagged: **a BLANK field dodges every check that keys on it** -- 3 rows
had NO uniprot, so they skipped all identity sweeps.
- **2 MORE wrong-gene anchors found + DROPPED** (both blank-uniprot, both anchored to MACHINERY, both real
  experimentally-validated effectors but mis-anchored; Teo's call = drop, not re-anchor):
    - **T3SS_23 YspE** -> YE_RS17825 = a HrpE/YscL-family T3SS APPARATUS protein (in the ysa operon), not the
      secreted Ysp. (It was found=yes because ssign emitted the apparatus protein, so the drop slightly DE-inflated
      recall -- more honest.)
    - **T6SS_21 Ssp2** -> SMDB11_RS11225 = T6SS ACCESSORY protein TagJ. Real Ssp2 = SMDB11_RS24150/SMA2264/
      WP_089196535.1 (proteomics MW 17867 confirms), but Bakta missed that ORF -> would be found=no anyway.
- **T6SS_13 YezP KEPT**: real effector, correctly anchored at the T6SS-cluster edge (next to TssM), Wang 2015 PLoS
  Pathog experimentally secretes it ("establish that YezP is a substrate secreted by T6SS-4"); just no UniProt entry.
- **T2SS_05 plaA FALSE ALARM (walked back)**: UniProt mislabels lpg2343's gene as sseJ and misassigns "plaA" to
  lpg2837(=PlaC). RefSeq + Flieger 2002 confirm lpg2343 = PlaA = the row. Row is CORRECT. (Lesson: check the
  gene-order index before "fixing" a flag.)
- **Are the source papers experimental? YES** (Teo asked): 85/90 citation quotes carry a direct experimental
  secretion signal (immunoblot/secretome/translocation); the ~5 exceptions are autotransporters whose secretion is
  experimentally established but whose quote describes structure (Ag43/NadA/YadA), a proteomic-secretome row with a
  name-listing quote (T2SS_01), and the YspE truncation. No row rests on a prediction/review/database paper.
- **Recall now 38/88 = 43%** (was 39/90): T1 16/18, T2 1/8, T3 5/22, T4 0/8, T5 10/19, T6 6/13. corrections.tsv = 336 rows.
- **SEPARATE-BENCHMARK FINDING (NOT the 88-row gold list -- flag for the effector-recovery-benchmark / 337-corpus):**
  `data/source_corpus/T3SS_verified.tsv` anchors the whole dispersed-Ysa-effector block (YspE/F/I/K/L/M) to a
  SEQUENTIAL placeholder run YE3556-YE3564 = apparatus/regulatory genes (ysaP pilotin, response regulator, ysrS
  sensor kinase). Agent verified the REAL loci are dispersed (yspI=YE2444, yspK=YE2447, yspM=YE3614). Exact YspE YE
  number is only in Matsumoto&Young 2006 Table 1 (paywalled). The 337-corpus Ysa block needs re-anchoring or dropping;
  does NOT affect the 88-row gold list (only YspE reached it, now dropped). **Trigger: before trusting 337-corpus T3SS recall.**

## (history) RESUME 2026-06-29 (eve) — gold v4 (90 rows); ALL-90 reachable recompute DONE; clean found×reachable 2x2

**ALL-90 REACHABLE RECOMPUTE DONE 2026-06-29 eve (Teo: "do the full reachable recompute, I don't trust this now").**
Replaced the curated `machinery_resolved.tsv` as the source for `nearest_machinery_gene/locus`, `distance_to_machinery_genes`,
`reachable_within_3` on ALL 90 rows with the machinery TXSScan actually DETECTED in the rerun (gene-order distance to the
nearest COGNATE component; `scripts/73_recompute_geometry.py`, renamed from the t1ss_t5ss version). Motivation: the curated
table is a partial component list that put effectors 100-3000 genes from "their" machinery while ssign emitted them via
proximity (CopN 520->1, BipC 11->3, TplE 1454->1, Tle1 265->1, TseM 453->2, YezP 3114->1). **found_by_ssign untouched ->
recall stays 39/90.**
- **Clean found × reachable 2x2 now:** found=yes/reach=yes 31, found=yes/reach=self 7 (self-secretors), **found=yes/reach=no 1
  (VirA only)**, found=no/reach=yes 9, found=no/reach=no 40, found=no/reach=self 2 (plpD, eae). `found ⊆ reachable` holds bar VirA.
- **VirA (T3SS_20) is the one genuine found-but-far:** a T3SS effector on the Shigella virulence plasmid, 24 genes from its
  T3SS apparatus but **1 gene from a T5aSS autotransporter**; ssign emitted it via DLP-extracellular 0.9994 + that non-cognate
  neighbor (dse_type_match=False). Real "recovered, but not by cognate proximity" datapoint, not a near-miss.
- **BUG found + fixed (mine):** the simplify "consistency" tweak that made `detected_components` stop at the first blank line
  was WRONG -- results.csv splits secretion systems into TWO subsections, `# Secretion Systems (with secreted proteins)` and
  `# Secretion Systems (other)`, blank-separated, and a system with no emitted substrate (often the T6SS) lives only in
  `(other)`. The terminator dropped `(other)`, so e.g. Ssp2 came back "no machinery" with 13 T6SS comps on its contig.
  Reverted to read-to-EOF + record_type filter (commented why it canNOT share `_read_emitted_loci`'s terminating reader).
- **Miss-reason categories for the figure (verified honest, status col in `geometry_recompute.tsv`):**
  (1) recovered = found (39); (2) reachable-but-missed = found=no/reach=yes 9 (cognate machinery detected + effector <=3 away,
  but the secretion-prediction/localization gate failed -- the TUNABLE misses); (3) unreachable, machinery detected, effector
  far = found=no/reach=no `recomputed` (NleD 2944 from LEE, SipA 5, BepA 220, VasX 78 -- proximity-filter limit); (4)
  unreachable, cognate NOT detected / cross-replicon = the 14 `no_cognate_machinery` (3 cross-contig: VceC/BtpB Brucella chr1
  effector vs VirB on chr2, tagA V.cholerae chr2 vs T2SS chr1; 11 absent: cag T4SS/Dot-Icm not modeled as pT4SSt, Citrobacter
  LEE undetected genome-wide [same IS-dense genome that lost NleD to a Bakta miss], HrpK1 Hrp T3SS, Serralysin) -- DETECTION
  limit / replicon split; (5) self-secretor not emitted = reach=self/found=no (plpD, eae). T3SS_30 NleD = `effector_not_in_order`
  (Bakta gene-miss, ORF absent from rerun) -> nearest="(effector ORF absent from rerun)", reach=no.
- Scripts converge as 65 -> 73 -> 65 (gold md5 stable, "reachable changed vs stored: 0"); all idempotent. corrections.tsv = 338 rows.

## (history) RESUME 2026-06-29 (pm) — gold list v4 (90 rows); high-churn re-review + RTX/T5SS/T6SS geometry DONE; recall 39/90

**HIGH-CHURN RE-REVIEW + GEOMETRY RECOMPUTE DONE 2026-06-29 pm (Teo: "confirm the sheet; re-review high-churn rows + the 32 geometry rows").**
- **2 wrong-gene errors found + fixed (recall 37->39/90).** Blind 3-agent re-review of the 30 rows with >=2 correction
  types, cross-checked by my own deterministic Pfam lookup (authoritative over the vote). An earlier `swap_up` had
  anchored two rows onto the WRONG gene:
    - **T6SS_18 Tle1 (EAEC 042):** old D3GUV9/EC042_RS24190 = Pfam PF00691 OmpA = the TagL PG-binding accessory
      (PDB 5M38/7BBA), NOT Tle1. Real Tle1 = **B7LG63 / EC042_RS24220** (560aa, PF09994 Tle1-like_cat), RefSeq
      "T6SS effector phospholipase Tle1-EAEC", beside its Tli1 immunity genes. Re-anchored (`reanchor` in scripts/65).
    - **T6SS_17 Tae4 (S. Tm 14028s):** old A0A0F6AX88/STM14_RS02020 = Pfam PF09867 TagF_N = a T6SS regulator, NOT the
      Tae4 amidase. Real Tae4 = **A0A0F6AX79 / STM14_RS01975** (161aa), RefSeq "T6SS amidase effector Tae4", beside Tai4.
      (agents split weak/ok/fix_id; the Pfam settles it.) Re-anchored.
  Both real effectors EMIT in the rerun (Tae4 DLP 0.99, Tle1 DLP 0.84, both DLP+DSE, nearby_ss=T6SSi) -> **found no->yes
  for both, so the wrong-gene anchoring had UNDER-stated recall by 2.** Verified the emissions are real (in the emitted
  secreted-proteins list), not a matcher artifact.
- **Deterministic cross-field check (scripts/72) on all 30 high-churn rows: 0 mechanical errors** — UniProt length matches
  the coord-encoded protein (28/30 ratio exactly 1.000, the 2 prtA off by a few aa = signal peptide); coords land on a real
  gene; found fresh. The identity (uniprot<->locus<->coords) is airtight; only the 2 Pfam-level wrong-gene errors slipped
  the length check (coords had been set to match the WRONG uniprot, so length was self-consistent).
- **Other re-review flags were pre-adjudicated, kept (3-agent reconciliation in gold_review2/sweep3_highchurn/):**
  5 fix_quote (T1SS_R02 apxIA, T1SS_R15 rtxA, T5SS_02 Ag43, T5SS_16 nadA, T5SS_17 yadA) + 1 weak (T5SS_03 iga, cross-species
  ref). All are rows where the protein IS a bona fide secreted substrate but the citation_quote describes structure/
  nomenclature rather than the export event; scripts/65 already carries `note` dispositions keeping them. **SURFACE TO TEO:
  T5SS_02 (Ag43) — all 3 agents say the current quote isn't verbatim in the cited 1999 paper; prior pass kept it as a valid
  autotransport statement. Decide whether to hunt a verbatim secretion sentence.** T4SS_02 Ceg14 = ok (2/3; 1 agent hit a
  2024 paywall).
- **Geometry recompute (scripts/73) — 34 rows now sourced from rerun-DETECTED components (gene-order distance), not the
  partial machinery_resolved table:** 13 RTX-T1SS (11 distances confirmed identical, apxIIA now concrete d453, Serralysin no
  T1SS on contig), 10 T5bSS (the big fix: "(self-secreting) d0" -> real TpsB distance; **8/10 reachable**, fhaB d4 + lspA1
  d153 outside; matches my earlier hand-derived NOTES values exactly), 7 T5a/T5cSS self-secretors relabeled honestly
  ("self", reach="self" not yes), 2 T5d/T5eSS self/no-model, + T6SS_17/18 (machinery_resolved gave d5/7=not-reachable but
  the effectors sit d1 from a detected T6SSi component and emit -> recomputed d1 reach=yes). **Geometry recompute does NOT
  touch found_by_ssign**; recall is driven by emission overlap only.
- **Recall = 39/90 = 43% (headline no-T5SS 29/71 = 41%)**: T1 16/18, T2 1/8, T3 6/23, T4 0/8, T5 10/19, T6 6/14. corrections.tsv
  = 173 audit rows. Scripts run order: `65 -> 73 -> 65` (reanchor sets coords; 73 reads them; 65 applies 73's geometry); all idempotent.
- **DEFERRED:** machinery_resolved.tsv is a CURATED PARTIAL component set; for full consistency ALL 90 rows' reachable/distance
  should be recomputed from rerun-detected components (only the 34 above were). Trigger: if reachable_within_3 is reported as a
  proximity ceiling in the paper. (Supersedes the older T5bSS-only reachability deferral.)

**NEXT: back-half of 4.4 = wire the 39/90 recall into figs 01-06 (task 1.7); then 4.5 docs. Then talk figure auto-generation (Teo queued this).**

## (history) RESUME 2026-06-29 (am) — gold list v3 FINAL (90 rows); geometry+overlap audit DONE

**GEOMETRY + OVERLAP found_by_ssign AUDIT DONE 2026-06-29 (task #19, OpenSpec 4.4a).** Built
`scripts/71_geometry_overlap_audit.py` + canonical `RerunIndex.emitted_overlap()` (the any-overlap rule Teo
fixed on 2026-06-29: any bp overlap with an EMITTED secreted protein = found; ORF-called-but-not-emitted or
Bakta-gene-miss = genuine fail). `scripts/65` now RECOMPUTES found_by_ssign for every row from its coords, so
the column can't go stale. Results: **geometry 0 true mis-anchors** (58 machinery-resolved all match; 21 RTX/
T5SS self-secretors have no resolved machinery instance; 11 effector-locus unindexed). **Molecule identity
reconciled**: 3 RefSeq↔INSDC contig aliases verified at lenratio 1.000 (PAO1 AE004091.2↔NC_002516.2,
B.pertussis NC_002929.2↔BX470248.1, +MC58) → all 90 reconcile, 0 un-reconcilable, NO naming-mismatch fails.
**6 found corrections** (all 3 matchers agree): apxIA/hlyA/ltxA no→yes (rerun_fullasm emits them, stored was
stale); TseM/Tae4/Tle1 yes→no (ORF annotated but ssign never emitted it). **Recall = 37/90 = 41% (headline
no-T5SS 27/71 = 38%)**: T1 16/18, T2 1/8, T3 6/23, T4 0/8, T5 10/19, T6 4/14. T4SS 0% = correct consequence of
the cross-replicon re-anchors (effectors far from VirB), not a bug. corrections.tsv now 81 rows (75 + 6
found_recompute). 4 simplify agents reviewed: 0 correctness bugs (overlap `>0` proven inert, min real overlap
293 bp); folded `reason` into emitted_overlap return, dropped a dead line, unrolled a walrus. **NEXT: back-half
of 4.4 = wire this 90-row recall into figs 01-06 (task 1.7); then 4.5 docs.** Committing now.

**RECALL SANITY CHECK 2026-06-29 (Teo flagged low T2SS/T5SS).** Verified NO genes dropped: every T2SS+T5SS
gold row cleanly overlaps its annotated rerun protein (thousands of bp). The misses are real non-emissions,
interpretable by subtype:
  - **T5SS canonical self-secretors found 7/7** (T5aSS espP/flu/iga/pic/sat 5/5 via `src=T5SS-self`; T5cSS
    nadA/yadA 2/2 on OM+SignalP). Exactly as expected.
  - **T5bSS (TPS, two-partner, NOT self-secreting) 3/10.** KEY: the TpsA substrate is RIGHT NEXT to the
    detected TpsB pore (gene-distance 1) in most misses, so it is NOT a proximity problem — it is the
    DeepLocPro localization gate. 5/7 misses are within ±3 of the pore but fail DLP extracellular>=0.80:
    cdiA d1 (0.79), lspA2 d1 (0.69), bcpA d3 (0.75) extracellular-but-just-under-cutoff; hpmA d1 (0.11) +
    shlA d1 (0.09) DLP calls the big TpsA exoprotein Outer Membrane. Only 2 are true proximity misses: fhaB
    d4 (just outside window, also OM) and lspA1 d153 (H. ducreyi's two LspA share one distant LspB). **Teo's
    call: do NOT change the pipeline — this is the correct consequence of a strict cutoff; report honestly.**
    Found T5bSS (cdrA/hmw1A/hxuA) all d1 + DLP 0.91-0.98. Paper-worthy line: T5bSS recall is DLP-gated, not
    distance-gated.
  - **T5SS non-canonical 0/2** (plpD T5dSS, eae T5eSS): OM, no self path. Architecture, as Teo predicted.
  - **T2SS 1/8** (Teo: "fine"): substrates genomically scattered, no T2SS component within +/-3 (`near=`
    empty) even though DLP called most extracellular 0.80-0.997. Known proximity-filter limit; motivates the
    classifier. Only lapA (next to machinery) emitted.
  - **DEFERRED (data-correctness, does NOT affect found/recall):** the 10 T5bSS gold rows are mislabeled
    `nearest_machinery_gene=(self-secreting)`, `reachable_within_3=yes`, `distance=0`. TPS substrates are NOT
    self-secreting; the real gene-distance to the cognate TpsB pore is 1-4 (153 for lspA1). `found_by_ssign`
    is computed from emission overlap, not this column, so numbers are unaffected. Proper fix needs curating
    each TpsB pore locus (machinery_resolved has no T5SS instances). **Trigger to revisit:** if a paper figure
    reports a proximity-reachability ceiling (then T5bSS true reachable = 8/10, not 10/10), curate TpsB loci.

**SHEET COVERAGE AUDIT 2026-06-29 (Teo: "confirm the sheet is correct; anything unchecked / high-churn?").**
Mapped every verification dimension of the 90-row `gold_list_final.tsv`. All FOUR independent dimensions swept:
  1. **Identity** (organism/accession/locus/uniprot/seq): 4-agent blind review (`gold_review` Jun26, 97 rows,
     25 fixes) + focused uniprot sweep (`gold_review2/sweep1_identity` Jun28, 24 rows). Applied via scripts/65.
  2. **Citation** (primary_ref + quote): 3-agent blind review (`gold_review2/sweep2_citation` Jun28, all 90):
     66 ok, 14 fix_quote, 4 fix_ref, 4 weak, 2 ties. ALL adjudicated+applied (weak/tie rows too).
  3. **Geometry** (nearest/distance/reachable): scripts/71 recompute — 58 machinery-resolved match exactly,
     0 mis-anchors.
  4. **found_by_ssign** (recall): scripts/71 any-overlap recompute, all 90, 6 fixes.
Rechecked the two scariest residuals this turn, both CLEAN: (a) the single `no_overlap_bakta_miss` row
(T3SS_30 NleD, D2TML3) — verified GENUINE Bakta gene-miss: gold span 638308-639015 = 708bp = 235aa (matches
NleD exactly), but rerun Bakta left a ~1kb annotation gap there (skipped locus GMFMKK_00600) flanked by an
integrase + 2 IS/transposases; coordinate is correct, honest found=no. (b) weak/tie citation rows — all carry
an applied correction. **Two genuine residual gaps, NEITHER corrupts recall** (recall = found_by_ssign, which
was recomputed for all 90):
  - **DEFER: 32 rows have UNVERIFIED geometry** (13 T1SS RTX + 8 T5SS no_machinery + 11 T5SS unindexed):
    distance/reachable never recomputed (machinery_resolved.tsv has no instance for RTX or T5SS). Subsumes the
    T5bSS self-secretor mislabel above. Trigger: only matters if a figure/paper prints distance/reachability for
    these types; fix needs TpsB-pore + RTX-machinery curation.
  - **OPTIONAL: 6 triple-corrected rows** (hlyA 4 correction-types; lspA2/TseH/prtA[R11]/YopJ_YopP/lspA1 3
    each) are the highest-churn = where a compounding error would hide. Each has a documented correction; a final
    human eyeball (or a focused 3-agent re-review) is the cheapest belt-and-suspenders check. 30 rows touched by
    >=2 correction types; events/row max = 4.

## (history) RESUME 2026-06-28 — Phase A gold list v3 FINAL (90 rows); both review sweeps DONE

**SECOND-PASS REVIEW IN PROGRESS (Teo asked 2026-06-26).** After the tier-B re-anchor, a direct look at the
4 "note" rows found all 4 were real fixes (committed a30ed53): T4SS_05 reanchor (ATU_RS23890 was the 83aa
Atu6188/P08061 mislabelled virE3; real effector = Atu6191/ATU_RS23905/A0A2Z2PK47, flips reachable->no),
T5SS_03 swap_up Q51163(fragment)->Q9K0B4(1815aa MC58), T5SS_13/14 swap_ref to Ward 2004 LspB-secretion paper
+ fill LspA2=Q9ZHL0. **gold_list_final.tsv now 71 confirmed + 19 corrected + 0 noted = 90.**

That the notes hid real defects (only caught because one agent noticed) prompted a full second-pass review of
the dimensions the blind review did NOT check. **Two agent sweeps over the 90 rows** (Agent tool, batched,
NOT Workflow). **HARNESS BUILT + COMMITTED 02162ed** under `gold_review2/`:
  - `scripts/66` precompute (`identity_signals.tsv`): cleared 66/90; flagged T3SS_29 (Q6T6T6 ratio 1.32),
    T3SS_30 (Q6RYA7 ratio 0.85) + 22 blank-accession rows; 0 dups.
  - `scripts/67` builds batches: sweep1_identity = the 24 flagged rows (3 batches s1_1..s1_3),
    sweep2_citation = all 90 grouped by ss_type (11 batches). Each sweep has its own anti-hallucination
    SCHEMA.md. `scripts/68` reconciles per sweep (generic; deterministic plurality + TIE, 3 agents/batch).
    simplify pass folded modal() into bench_io. `scripts/69` deterministically re-validates every
    agent-proposed accession (length+organism vs UniProt) — the machine-checked safety net for adjudication.
  - **SWEEP 1 (identity) DONE 2026-06-28.** 9 agents (3/batch, 3/3 coverage) -> reconcile (`68 sweep1_identity`)
    -> validate (`69`). Outcome folded into `scripts/65` + re-run: **14 accession ADDs** (blank rows, gene-name
    confirmed, length-matched), **7 REANCHORS**, **1 organism fix** (T6SS_13), **3 none_exists** kept blank
    (T3SS_23, T6SS_21, + T6SS_13 accession). gold_list_final = **55 confirmed + 35 corrected = 90** (per-type
    unchanged T1 18 / T2 8 / T3 23 / T4 8 / T5 19 / T6 14). corrections.tsv = 48 rows. Committed 6b74fd5.
  - **!!! SYSTEMIC MIS-ANCHORING FOUND.** 7 of 23 T3SS rows had `effector_locus` pointing to the WRONG gene
    (a transposase/hydrolase/regulator, not the named effector); agents caught it via NCBI feature-table.
    Reanchored to the true effector locus (length confirmed vs index + UniProt): T3SS_16 (hrpS->HrpK1
    PSPTO_RS07380/G3XDB3), T3SS_27 (ydgJ->NleF E2348C_RS07680/B7UR63), T3SS_28 (HTH->EspA ROD_RS14675/D2TKE3),
    T3SS_29 (IS110->NleA ROD_RS01735/D2TJZ4), T3SS_30 (hydrolase->NleD ROD_RS02910/D2TML3), T3SS_31
    (IS3->NleB1 ROD_RS05475/D2TT37, was spuriously found=yes!), T3SS_32 (transglycosylase->EspJ
    ROD_RS24095/D2TRY1). **OPEN QUESTION for Teo: run a deterministic locus-identity audit on ALL 90 rows**
    (does each gold effector_locus actually encode the named gene?) before finalizing -- the build (scripts/63)
    clearly mis-mapped some loci; non-flagged rows with a coincidentally length-matching wrong gene would NOT
    have been caught by the precompute. Decide: locus-audit first, or sweep-2-citation first, or both.
  - **LOCUS-IDENTITY AUDIT DONE 2026-06-28 (Teo chose this before sweep 2).** `scripts/70` (+ bench_io
    `locus_index`/`norm_locus` helpers) deterministically maps each accession's UniProt ordered-locus-name AND
    the gold effector_locus through the gene-order index; a `match` = same gene. Result: **54 match, 0 MISMATCH**
    -> the systemic mis-anchoring was confined to the 7 T3SS rows already fixed; nothing new. Residual 36 are
    UNVERIFIABLE-BY-INDEX, no mis-anchor signal: 21 no_oln (reviewed Swiss-Prot toxins, distinctive lengths,
    gene-name+length already checked), 11 T5SS gold_locus_unindexed + 1 espP (phase1 index doesn't cover those
    genomes; graded from rerun; large distinctive autotransporters/TPS), 3 no_acc (none_exists). Output
    `locus_identity_audit.tsv`. Low residual risk; optional follow-up = cross-check the 12 T5SS vs the rerun
    Bakta annotation. Committed 3edbb5a.
  - **SWEEP 2 (citation, all 90) DONE 2026-06-28.** Re-ran `scripts/67` (batches carry sweep-1 accessions) -> 33
    agents (11 batches x 3, full 3/3 coverage) -> `68 sweep2_citation`: **66 ok / 32 adjudicate** (14 fix_quote,
    4 fix_ref, 4 weak, 2 ties, 8 majority-ok). Signal was overwhelmingly QUOTE quality, not membership: many quotes
    were garbled (placeholder `[Bartonella]`, dropped species names, `uropathogenicCFT073` typo, an EMPTY T4SS_11
    quote, a cross-wired FhaC-on-a-Proteus-HpmA row) or described function/structure not secretion. Folded into
    `scripts/65` (new `fix_quote` action; STATUS+=fix_quote): **9 quote-only fixes + 3 ref swaps (YopJ->Brodsky2008,
    hlyA->Koronakis1991, prtA->Ghigo1992, all verbatim-verified via Europe PMC this session) + 6 existing
    dispositions extended with citation_quote + 4 weak-but-valid kept as notes** (rtxA T1SS, nadA/yadA T5cSS TAAs,
    flu T5aSS -- bona fide substrates, current/structure quote kept, NO drops). **Caught a latent bug:** T6SS_03
    was duplicated in DISPOSITIONS (sweep-1 swap_up shadowed the first-pass mBio swap_ref -> ref correction silently
    lost); consolidated. gold_list_final v3 = **39 confirmed + 47 corrected + 4 confirmed_note = 90** (per-type
    unchanged T1 18/T2 8/T3 23/T4 8/T5 19/T6 14). corrections.tsv = 75 rows. simplify: factored the shared Bowen-2003
    prtA quote into `PRTA_RTX_QUOTE` (used by T1SS_06 + T1SS_R12). Committed 18d48a4.
  - Sweep-2 weak-but-kept rows to remember (citation describes structure, protein IS a valid substrate): T3SS_17
    VopS (left as-is, AMPylating T3SS effector), T1SS_R15 rtxA, T5SS_16 nadA, T5SS_17 yadA, T5SS_02 flu.

**THEN (still pending):** 1.7 regen figures from gold_list_final + 4.4 recall rebuild + 4.5 docs.

## (history) Phase A gold list FINALIZED (90 rows) — tier-B re-anchor

**Phase A strategy changed (Teo, 2026-06-26):** instead of re-verifying all 337 effectors, collapse to
**ONE representative substrate per secretion-system instance** (benchmark scores "found a system + a
substrate", not how many). The closest substrate to the instance's machinery is kept (mirrors ssign's
+/-N proximity). This 1:1 list is the **new trusted GOLD LIST** for the paper, replacing the 337 table.

**Built (committed):** `scripts/60` (337-row ledger, superseded), `scripts/61` (mechanical DOI+UniProt
check: 0 dead DOIs / 0 dead accessions across 337), `scripts/62` (early representative selector, superseded),
**`scripts/63_build_gold_list.py` -> `data/phase2/verification_phase_a/gold_list.tsv` = 97 trusted
system<->substrate pairs** (T1 18, T2 11, T3 25, T4 9, T5 19, T6 15). Fully enriched: host, gene,
uniprot, genome+coords, nearest machinery + gene-distance, reachability, found_by_ssign,
`stage_replicons`, DOI, quote.

**KEY FIX (2026-06-26):** an earlier 79-row version restricted to the 90 machinery instances, which
SILENTLY DROPPED 13 verified T1SS effectors (hlyA/apxIA/ltxA/lktA/rtxA/apxII/apxIII RTX toxins, prtA-C +
Serralysin proteases, hasA) -- their systems are T1SS_R## RESCUE instances, not in the original
machinery audit. Teo caught it (his fig 06 showed T1SS n=18, found 13; my 79-row list had T1SS n=5).
Rebuilt `scripts/63` from the SAME answer key `52_system_recall` uses (instance set + `nearest_dist`
from `ceiling_per_effector` incl. rescue; cleaned effectors + `ssign_call` from `actual_per_effector`
via clean_dataset; provenance from positives_all). Per instance keep the CLOSEST effector, preferring
one that preserves the instance found/reachable status. **Now reconciles with fig 06 EXACTLY** for
T1-T6 in its 3-segment logic (found/reach_miss/unreach): T1 13/3/2, T2 2/0/9, T3 8/3/14, T4 0/4/5,
T6 7/3/5. T5SS uses the REAL rerun emission (10/19 found) not 53's inferred 15 (the 4.2 over-count);
T5d/e (plpD/eae) flagged no-TXSScan-model = unreachable. The "discrepancy" was NOT a reachability bug
(my per-substrate reachable matched the answer key 56/56); it was the 90-only restriction + me comparing
reachable-vs-found columns.

**STAGING RULE (Teo, important):** multi-replicon hosts must be staged WHOLE for ssign, never split
into separate runs (the V. cholerae chr2 / Ralstonia megaplasmid bug class). `gold_list.stage_replicons`
lists the replicons to co-stage. **7 genomes need it:** V. cholerae (T2SS_02, T6SS_01), Ralstonia
(T3SS_13), Yersinia ent. (T3SS_23), B. melitensis (T4SS_06), B. abortus (T4SS_08), B. thailandensis
(T6SS_10). 4 are cross-replicon (machinery+substrate on different molecules -> substrate unreachable by
proximity, flagged). Reachable within +/-3 genes: 39/79.

**Multi-agent blind review (DONE 2026-06-26):** 97 gold rows in 10 batches, **4 independent agents each
via the Agent tool** (NOT workflow; the overnight failure was the workflow `args`-as-string char-split
bug). Each agent read only `gold_review/batches/<id>.tsv` + `SCHEMA.md`, hard-failed if missing, used
PubMed MCP + UniProt/Crossref, wrote `gold_review/verdicts/agent{k}_<id>.tsv` (16-col schema). 40/40
passes, 4/4 coverage (one agent4_T2SS_1 dropped T2SS_12 -> 3/4 there; one T4SS_1 agent3 hit a usage-policy
block, re-run clean).

**Reconciliation (`scripts/64_reconcile_gold_review.py` -> `gold_review/reconciliation.tsv`):** majority-vote
final verdict per row + agreement label + row_id-drift/coverage check. **71 confirmed / 25 fix / 1 split /
0 drop.** The split T4SS_03 + 3 dead-accession fixes (T2SS_10/12, T6SS_12) were DROPPED per Teo.

**Corrections applied (`scripts/65_apply_corrections.py` -> `data/phase2/verification_phase_a/gold_list_final.tsv`
+ `gold_review/corrections.tsv`):** disposition table is the human-adjudicated outcome; raw gold_list.tsv
untouched. **gold_list_final.tsv = 90 rows, ALL SCORABLE (0 holds): 71 confirmed + 15 corrected + 4 noted;
7 dropped.** Per-type kept: T1 18, T2 8, T3 23, T4 8, T5 19, T6 14. Tier-A fixes = 9 organism relabels +
2 uniprot swaps (T3SS_17/21) + 1 ref swap (T6SS_03).

**Tier-B re-anchor pass DONE (2026-06-26, Teo: "high quality, drop if you can't"):** the 7 held rows resolved.
`scripts/65` gained a `reanchor` action that recomputes coords + nearest-machinery + gene-distance +
reachability from the SAME gene-order index/machinery tables 62/63 use (no hand-typed numbers); found_by_ssign
read from the rerun by coordinate join (RerunIndex), like the T5SS rows.
  - **T1SS_04 -> relabel** (NOT mis-anchored): BX470248.1 ORGANISM line = *B. pertussis Tohama I*; cya/BP0760/
    P0DKX7 + CyaB neighbour all consistent. Was a mislabel. Scorable (found yes).
  - **T2SS_05 -> reanchor** plaA = LPG_RS11775/lpg2343/Q5ZT22 (old LPG_RS04595 was ravI/lpg0926, a Dot/Icm
    effector). Far from Lsp machinery (nearest LspM, 472 genes) -> reachable no, found no.
  - **T4SS_08 -> reanchor** BtpB = BAB1_0756/BAB_RS19540/Q2YN91 (TIR) on chr I (old was BAB1_0782/DUF2069).
    VirB machinery on chr II -> cross-replicon, not reachable. Both chrs staged.
  - **T4SS_09 -> reanchor** VceC = BAB1_1058/BAB_RS20990/Q2YQ34 on chr I (old was chr II). Cross-replicon, not
    reachable. Both chrs staged.
  - **T2SS_08 -> DROP:** instance machinery is Xanthomonas Xps (NC_007508.1); cited substrate is V. cholerae
    lipase. Substrate organism != machinery organism -> cannot re-anchor coherently.
  - **T3SS_14 -> DROP:** listed locus RS_RS21300/RSp0871 is *hrpD* (machinery), not RipJ -- anchored onto the
    machinery (spurious dist 1); deleted RipJ accession's only successor is a different assembly, not GMI1000.
  - **T3SS_13 -> DROP** (Teo: "drop the one with the dead uniprot"): PopC's cited A0ABY6NJB5 was deleted from
    UniProt. A live Reviewed replacement EXISTS if reinstated (PopC = RSp0875/RS_RS21320/Q9RBS2/POPC_RALN1 on
    the GMI1000 megaplasmid NC_003296.1, adjacent to HrcC, reachable=yes) -- recorded in the drop basis.

**>>> NEXT:** **1.7** regen figures from `gold_list_final.tsv` (record before/after delta) + the held **4.4**
recall rebuild on the finalized 90-row gold list, then **4.5** docs + close tasks 70/71.

**DEFERRED (simplify):** `scripts/65 Reanchor.__init__`, `62 ordinal_index`, `63 gene_order` are three
near-identical gene-order-index loaders. Fold into one `bench_io.gene_order_index()` -> locus ->
[(rec, ordinal, start, end, strand)]. Deferred because it touches 62/63 (already produced committed
artifacts); do it with a re-run + reconcile check, not right before a commit. Trigger: next time 62/63/65
are edited for other reasons.

## (superseded) RESUME (2026-06-25 late) — Phase B reconcile 4.1/4.2/4.3 DONE, 4.4 HELD for Phase A

**Next action = Phase A** (`benchmark-final-validation` tasks 1.1-1.7): the deep re-audit of the 337
secreted-protein effectors + every DOI. UNSTARTED, needs Teo's explicit go-ahead. Everything else in
Phase B is done or deliberately held. Why held: 4.4 (recall figure rebuild) waits on Phase A because
the re-audit will change the effector set, so rebuilding recall now would be throwaway (Teo's call
2026-06-25).

**Phase B (4.x) status:** 4.1 DONE (RTX toxins confirmed emit from fullasm rerun), 4.2 DONE (T5SS
annotation), 4.3 DONE (enrichment), **4.4 RECONCILED but figure-rebuild HELD**, 4.5 pending (docs).
Reconciliation result (preserve): coordinate-joined all 496 testable effectors to the tier-2 rerun =
**471/481 agree (98%)**; **+3** toxins now emit (hlyA/apxIA/ltxA, the fragment-bug fix), **-6 in-figure**
the rerun no longer emits (RipJ, NleB [T3SS]; Tae4_Stm, Tle1_Sci1, BopE, TseM [T6SS], clean-Bakta
borderline shifts), 15 unjoinable (genomes outside the 57-panel). Full detail in `benchmark-final-
validation` tasks.md 4.4. The join tool is `scripts/rerun_coords.py` (reusable when 4.4 resumes).

**State:** the `benchmark-final-validation` tier-2 CX3 rerun is COMPLETE + healthy. 57 distinct
genomes, clean-slate Bakta re-annotation, T3SS-on, `--enrichment-stats`, annotation-on. Results
extracted at `validation_sweeps/benchmark/rerun/<genome>/results/` (+ `rerun_fullasm/<g>/` for the 4
RTX/Serralysin genomes re-run on full assemblies). 52/57 emitted secreted + enrichment; 3671 secreted
proteins fleet-wide.

**Branch:** `enrichment-circular-shift-per-run` (NOT main, all recent work here). Only tracked dirty
file is NOTES.md.

**UNTRACKED run data — do NOT `git add`** (none of it is gitignored yet): `rerun/`,
`rerun_results.tgz`, `rerun/rerun_nulls.tgz`, `cx3_t3ss/`, `cx3_enrichment/`, `*.tgz`. GB-scale.

**Done this session:** reconstructed the pooled cross-genome enrichment figure post-hoc (the benchmark
batch runs one genome per `ssign run`, so MultiGenomeRunner's pooled figure was never auto-emitted).
Re-tarred the 57 `*_enrichment_nulls.npz` from CX3, ran `rerun/_pool_enrichment.py` (throwaway, calls
the production `pool_and_plot_enrichment`). Output `rerun/_pooled/{pooled_enrichment_stats.tsv,
figures/pooled_enrichment_null_distributions.png}`. Result: T1/T2/T3/T5a/T5b/T6 significantly enriched
(DLP+DSE); **T3SS DLP 5.9x *** with DSE correctly empty** (validates bug-fix #4 fleet-wide); T4SS +
T5cSS n.s. (low count). Effectively answers task #70 at the pooled level.

### The 4 active OpenSpec changes (what each scopes)
| change | status | scope |
|---|---|---|
| `benchmark-final-validation` | 4/20 | THE deep audit + tier-2 rerun + reconcile. Phase A (1.1-1.7) = re-verify 100% of 337 effectors + every DOI = the secreted-protein deep audit, UNSTARTED, held for go-ahead. Phase B rerun DONE. |
| `effector-recovery-benchmark` | 45/51 | original benchmark: verified gold set, ssign-INDEPENDENT machinery answer key (the SYSTEMS audit, lit-only, NO MacSyFinder), proximity ceiling, actual recall |
| `figures-v2` | 0/20 | the auto figure-gen redesign ONLY (appendage exclusion, T5a/c split, autotransporter fig, COG/KEGG/EggNOG functional figs, fold-bar enrichment, genome-group pooled 01-05). NOT implemented, paused for poster |
| `secretion-classifier-dataset` | 24/29 | classifier training set (audit 347 predicted rows, instance assignment, T5SS sourcing, feature matrix) |

"Deep proper audit of systems + secreted proteins" = SYSTEMS in `effector-recovery-benchmark`
(machinery answer key) + EFFECTORS in `benchmark-final-validation` Phase A (the unstarted 337-row pass).

### ⚠ BLOCKER FOUND 2026-06-25: rerun used FRAGMENT inputs for 4 genomes (script 58 bug)
Reconciliation 4.1 surfaced a staging defect. `scripts/58_stage_rerun_inputs.py` always reuses
`inputs_gb/<acc>.gbff` and never consults `inputs_gb_fullasm/`. For genomes whose corpus
`refseq_genome` accession is a plasmid/WGS-contig, the rerun ran ssign on that FRAGMENT, not the
full assembly that `scripts/50_stage_full_assembly.py` had already staged in `inputs_gb_fullasm/`.
The rerun thus reproduced the exact partial-assembly staging defect the T1SS fix was meant to cure.

Affected (rerun-input CDS -> full-assembly CDS, emitted): 
- `NZ_CP031766.1` hlyA (P08715, T1SS): 175 (plasmid) -> 5220, emitted 0
- `NZ_CBDBTK010000022.1` apxIA (P55128, T1SS): 33 (contig) -> 2192, emitted 0
- `NZ_JABJZG010000001.1` ltxA (P16462, T1SS): 241 (contig) -> 2164, emitted 0
- `NZ_JBCGCZ010000007.1` Serralysin (P23694, T1SS): 256 (contig) -> 4947, emitted 4 (toxin not among them)
NOT affected: `NC_010939.1` apxIIA (inputs_gb 2140 ≈ fullasm 2158, already ~complete); `NZ_SMAM01000003.1`
(not in rerun, no current effector). `NZ_CP006957.3` lktA was already the full chromosome.

**4.1 status:** of the 4 staging-fix toxins {hlyA,apxIA,ltxA,lktA}, only **lktA emits** (P0C082 NOT P55117;
UniProt drift — clean_dataset.T1SS_STAGING_FIX still lists the stale P55117). lktA in NZ_CP006957.3 =
locus OKNMBO_02445 "Leukotoxin", DLP 0.988, DLP+DSE, emitted ✓. hlyA/apxIA/ltxA cannot be confirmed from
this rerun (wrong input). Side note: apxIIA (P15377, full asm, DLP 0.98) did NOT emit = proximity/machinery
miss, a real 4.4 recall datapoint.

**Coordinate-join method (works):** rerun raw CSV has contig+start+end+strand; join answer-key effector
(placement_start/stop or start/stop) to the max-overlap rerun protein, check if its locus_tag is in the
emitted `*_results.csv`. Contig name == RefSeq accession for single-replicon assemblies. Near-exact overlap.

**Resolution (Teo chose: fix staging + re-run 4 on CX3):**
- **Code fix DONE:** `scripts/58_stage_rerun_inputs.py` now has `stage_whole_assembly()` — prefers
  `inputs_gb_fullasm/<acc>.gbff` over the fragment for single staged units (unit-tested; idempotent;
  flags copied files as NEW for CX3 transfer). rerun_inputs.txt stems are unchanged (only file content).
- **Targeted re-run STAGED (awaiting Teo's qsub):** `scripts/cx3/rerun_fullasm_batch.txt` (the 4 stems) +
  `scripts/cx3/SUBMIT_rerun_fullasm.sh` (one job, `INPUT_DIR_GB=inputs_gb_fullasm`, `RUN_TAG=rerun_fullasm`,
  config identical to the main rerun: genbank+REANNOTATE+T3SS+whole-genome+ENRICH+ANNOT). Separate RUN_TAG
  so it never clobbers `rerun/`. Results -> `~/runs/benchmark_phase2/rerun_fullasm/<unit>/`.
- **CX3 steps for Teo:** (1) `git pull` on CX3; (2) scp the 4 full-assembly gbffs to CX3
  `inputs_gb_fullasm/` (they were NOT in the original rsync — see the scp line in SUBMIT_rerun_fullasm.sh
  header); (3) `bash validation_sweeps/benchmark/scripts/cx3/SUBMIT_rerun_fullasm.sh`; (4) retrieve to
  `validation_sweeps/benchmark/rerun_fullasm/<unit>/`. NB Diamond PATH symlink must still be in place
  (ANNOT=1 hits eggnog env) — see the Diamond section below.
- **On return:** confirm hlyA/apxIA/ltxA/Serralysin emit via the coordinate join; if the 3 RTX toxins
  emit, retire APPLY_T1SS_STAGING_FIX (flag stays as fallback) + fix its stale P55117->P0C082 lktA id.
- **Changes are uncommitted** (working tree on `enrichment-circular-shift-per-run`): script 58 + 2 cx3
  files + NOTES + openspec tasks.md. Commit when Teo is ready.

### Session progress 2026-06-25 (cont.) — 4.1 staged for CX3, 4.3 DONE, 4.4 scoped
- **4.3 DONE** (enrichment sane fleet-wide): scanned all 57 `*_enrichment_stats.tsv` (444 rows, 0 numeric
  issues; every p/q in [0,1], every fold finite). Bug-fix #4 holds in every genome (21/21 T3SS DSE rows
  blanked, 21/21 T3SS DLP populated, 7 sig). Pooled: T3SS DLP fold 5.90 q=2e-4 sig, no T3SS DSE row.
  openspec 4.3 ticked. (formal "close task 70" = a docs line in 4.5.)
- **4.1 staged for CX3** (see BLOCKER block above): awaiting `rerun_fullasm/` for hlyA/apxIA/ltxA/Serralysin.
- **4.4 scoping (DO THIS NEXT, with rerun_fullasm in hand):** the recall figures (01 effector-level, 06
  system-instance) are built from `data/phase2/actual_per_effector.<tag>.tsv` (582 rows; tags
  `panel_genbank_default` + `panel_genbank_t3ss`) filtered by `clean_dataset`, NOT from
  `positives_all.testable` (that col is a sparse 13-row T1SS leftover — do not use it). Join keys in
  actual_per_effector: `unit_id` (genome), `effector_locus` (RefSeq locus), `ssign_call`
  (emitted_secreted/not_emitted), `testable`, `reachable_n3`. hlyA already shows `ssign_call=not_emitted`
  there too -> the fragment-input problem PREDATES the rerun; clean_dataset's APPLY_T1SS_STAGING_FIX was
  the workaround. **4.4 pipeline to build:** (a) for each (unit_id, effector_locus) parse the RefSeq gbff
  (`inputs_gb/<unit>.gbff`, or `inputs_gb_fullasm/` for the 4) to map effector_locus -> (contig,start,end);
  (b) coordinate-join to the rerun proteins (`rerun[_fullasm]/<unit>/results/<unit>_results_raw.csv` has
  contig+start+end; emitted = locus in `<unit>_results.csv`); (c) per-SS-type found/not vs the old
  actual_per_effector ssign_call; (d) compare to figure 06; (e) regen poster 01-06 ONCE with all 57
  (prefer rerun_fullasm/<g> over rerun/<g>). Coordinate-join method validated (near-exact overlaps).
- **4.2 DONE** (T5SS annotation from rerun): built reusable `scripts/rerun_coords.py` (`RerunIndex`,
  coordinate-join on (contig,start,stop) since Bakta renames loci; prefers `rerun_fullasm/` once it
  lands) and rewired `annotation_correctness.py` to source T5SS from `rerun/` (killed the dead
  `/tmp/ssign_fleet_67` path). Corpus has 23 T5SS effectors; **13 emit + graded** (NOT the 15 "found"
  53_t5ss_recall estimated). Re-source reproduced the old 4 verdicts (espP/pic/yadA/cdrA) = join
  validated. **All 13 overall-correct**; only T5SS per-tool misses are Pfam on the 2 TAAs (yadA wrong,
  nadA partial). Pfam channel now `eggnog__pfam_ids` (rerun Bakta emits no `pfam_ids`; EggNOG gives the
  same short-name format). Sheet 55 rows (42 proximity + 13 T5SS); 07 = 45/54 (83%), 08 = T5SS 13/13.
  T5SS verdicts gene-keyed (`T5SS_V`) so they survive Bakta locus drift.
  - **4.4 datapoint (recall discrepancy):** 8 T5SS corpus effectors do NOT emit in the rerun despite
    `53_t5ss_recall` marking 5 of them "found" via *local MacSyFinder + inferred* emission:
    plpD, eae (no TXSScan model), **shlA, hpmA, lspA1, lspA2, cdiA, bcpA** (T5bSS TpsA substrates not
    recovered). `53_t5ss_recall`'s "found=15" overcounts real emission (=10 of its 15); fold this into
    the 4.4 T5SS recall reconciliation. (prn, fhaB: genomes not in panel.)
  - **DEFERRED (4.4 follow-up):** the 42 proximity rows in `annotation_accuracy_sheet.tsv` still carry
    fleet_67 provenance strings (the hand-authored verdicts are biology judgments, unaffected, so 07/08
    are correct); a full rerun re-source of the proximity panel via `rerun_coords` belongs in 4.4.

### Open work inventory (pick up here)
1. **Reconciliation** (benchmark-final-validation tasks 4.1-4.5): 4.1 confirm the 4 RTX T1SS toxins
   (hlyA/apxIA/ltxA/lktA = P08715/P55128/P16462/P55117) emit from rerun -> retire `clean_dataset`
   APPLY_T1SS_STAGING_FIX; 4.2 grade 15 T5SS annotation; 4.3 enrichment fleet-wide (pooled fig already
   shows OK) -> close #70; 4.4 reconcile rerun recall vs figure 06; 4.5 docs, close #70+#71. Read
   `scripts/clean_dataset.py`, `scripts/52_system_recall.py`, `analysis/annotation_correctness.py` first.
2. **Phase A** (benchmark-final-validation 1.1-1.7): the 337-effector + DOI deep re-audit. Large
   multi-agent run, HELD for Teo's go-ahead. Reuses the tasks-74-84 blind-agent method.
3. **Poster figures:** `validation_sweeps/benchmark/figures/` has only `01`-`04`, STALE (pre-rerun). The
   openspec recall set is `01`-`06` + annotation `07/08`. Regenerating from the rerun = reconciliation 4.4 + 4.2.
4. **figures-v2** (0/20): auto figure-gen redesign, paused for poster. Self-contained in its design.md.
5. **Task #97:** proper Diamond PATH fix in run_benchmark_batch.pbs + `ssign doctor` Bakta-subtool
   version-check (tactical `ln -sf` already applied on CX3, see Diamond section below).
6. **Pre-archive:** sync 66->57 / 16->15 numbers in benchmark-final-validation proposal/specs (see panel
   section below) before `/opsx:archive`.

### Key file pointers
- Rerun per-genome: `rerun/<g>/results/{<g>_results.csv emitted, <g>_results_raw.csv all+annot,
  <g>_enrichment_stats.tsv, <g>_enrichment_nulls.npz}`
- Pooled figure (this session): `rerun/_pooled/figures/pooled_enrichment_null_distributions.png`
- Panel: `data/phase2/rerun_panel_manifest.tsv` (57 distinct), `data/phase2/rerun_inputs.txt`

## CX3 gotcha: Bakta 1.12 needs Diamond >= 2.1.21 (clean-slate rerun blocker, 2026-06-24)
- **Symptom:** rerun batch jobs die at step 2 (extract_proteins) in ~30s: "Wrong Diamond version
  installed. Please install Diamond version v2.1.21 or ...". Hits any Bakta run (we almost always
  re-annotate with Bakta, NOT --use-input-annotations).
- **Cause:** REGRESSION + a PATH-shadowing trap. `bakta/utils.py` `DEPENDENCY_DIAMOND=(2,1,21)`; Bakta
  1.12 needs Diamond >= 2.1.21. `bakta-deps` Diamond was actually already 2.2.0, BUT the
  `run_benchmark_batch.pbs` `ANNOT=1` block PREPENDS `~/.conda/envs/eggnog/bin` to PATH, and eggnog's
  bundled Diamond is OLDER than 2.1.21, so Bakta resolved eggnog's old diamond first. `which -a diamond`
  in a plain shell hid this (bakta-deps was first there; only the ANNOT PATH prepend reorders it).
- **Tactical fix (CX3):** `ln -sf ~/.conda/envs/bakta-deps/bin/diamond ~/.conda/envs/eggnog/bin/diamond`
  (eggnog-mapper tolerates diamond 2.2.0). Proper fix (deferred): in run_benchmark_batch.pbs, prepend a
  diamond-only dir (symlink to bakta-deps diamond) AFTER the annot block so the right diamond always
  wins without touching the eggnog env; and have `ssign doctor` version-check Bakta's Diamond.
- **Fix:** `module load Mamba/23.11.0-0; mamba install -n bakta-deps -c bioconda 'diamond>=2.1.21' -y`.
- **Follow-up (deferred):** `ssign doctor` should version-check Bakta's sub-tools (Diamond), it only
  checks the bakta binary exists, so it reported "All checks passed" while Bakta was unrunnable.

## benchmark-final-validation: panel is 59 genomes, not the proposed 66
- **What:** `scripts/57_build_rerun_panel.py` computed the rerun panel = **59 genomes** (52 staged-to-run
  + 7 stage-from-cache), **15** dropped, **0** needing a network fetch. The proposal/design/specs still
  say "66 genomes / 15 missing T5SS / 16 dropped" (the directive's estimate).
- **Why the estimate was high:** 13 of 20 T5SS genomes already overlap staged proximity genomes (the
  directive assumed only 5), and every genuinely-missing T5SS genome is already in `refseq_cache`.
- **Trigger to revisit:** before `/opsx:archive` of this change, sync the 66→59 / 16→15 numbers into
  `proposal.md` + `specs/tier2-benchmark-rerun/spec.md` so the contract matches `rerun_panel_manifest.tsv`.
- **Also deferred:** task 2.2b (stage the 7 cache-not-staged T5SS genomes into `inputs_gb/` + collapse the
  2 same-assembly merges) before the CX3 submit dry-run (2.3).
- **57 distinct genomes** (not 59): PAO1 (`AE004091`+`NC_002516.2`) and A. fabrum C58 (`NC_003063.2`
  +`NC_003065.3`) are each one assembly split across units; run each as ONE whole-assembly input.
- **Corpus organism mislabels** found + corrected in the manifest via NCBI: `NC_006351`=B. pseudomallei
  K96243 (corpus said B. thailandensis), `NC_013716`=Citrobacter rodentium (said EPEC), `NC_007508`=
  Xanthomonas euvesicatoria (said P. aeruginosa). Worth a corpus-level fix during Phase A re-verification.
- **inputs_gb -> CX3:** don't git-track the 537 MB dir; only ~9 new/merged `.gbff` are needed on CX3
  (7 staged + PAO1 + C58), tar just those and scp.

## RESUME HERE (2026-06-24): apply OpenSpec change `figures-v2`
Teo stepped away to revamp poster figures; resume the figure overhaul after.
- **What:** `openspec/changes/figures-v2/` — APPROVED proposal, 0/20 tasks. All decisions + verified
  facts are in its `design.md` (self-contained). Next step: `/opsx:apply figures-v2`.
- **Branch:** `enrichment-circular-shift-per-run` (main is behind; all recent work is here).
- **Decisions locked (design.md):** exclude ALL non-secretion appendages by default
  (Flagellum,Tad,T4aP,T4bP,MSH,ComM,Archaeal-T4P; keep T1-T6SS/T9SS/pT4SSi/pT4SSt); COG figures use
  category NAMES not letters; apparatus-associated secreted proteins -> labeled bucket (not filtered).
- **Verified facts (don't re-investigate):** `outer_membrane_prob` is in cross_validate predictions
  (line 277) but dropped at `integrate_annotations` -> surface it; `cog_category`/`kegg_ko` already in
  the integrated CSV when EggNOG ran; `dlp_extracellular_prob` already there.
- **Already shipped this session (commit 5152c0a):** dropped `04_signalp_by_type` + `05_tool_coverage`;
  SS variant labels collapse via `display_type`; per-run set is now `01`-`05`.
- **Multi-genome enrichment bug FIXED (d7435a8):** `_pool_segment_b_inputs` now pools the whole proteome
  for DLP/DSE under `--enrichment-stats`. Validation check: batched Salmonella T3SS null_mean should be
  ~0.55 (matching single-genome), not 0.08.
- **Work order:** task groups 1-6 in tasks.md (1 data+appendage-exclusion -> 2 T5a/c split+autotransporter
  fig -> 3 functional COG/KEGG/EggNOG/consensus candidates -> 4 enrichment fold-bars+group genome-wide
  null -> 5 genome-group pooled 01-05 -> 6 docs/tests/validate). CHECK IN with Teo after groups 1-2 before
  building the 10 functional candidates.
- Untracked artifacts kept locally: `validation_sweeps/cx3_t3ss/` (the 1/4/20-genome run figures+stats).

## CX3 validation runs PENDING (set 2026-06-22) — figure-revamp + t3ss-on-by-default

Branch `enrichment-circular-shift-per-run` @ b6805bf carries BOTH changes (figure-revamp
COMPLETE+archived-pending; t3ss-on-by-default groups 1-4 done, group 5 = these runs).
Three CX3 runs to validate (all `--small --enrichment-stats`, T3SS now on by default):
- **1 genome**: `$HOME/ssign-tutorial/salmonella_typhimurium_lt2.gbff` (T3SS) — satisfies t3ss task 5.1.
- **4 genome**: `--tutorial-all` (k12/pao1/salmonella/vibrio).
- **20 genome**: T3SS-rich benchmark set, tarball at `/tmp/ssign_benchmark_20.tgz` (56MB, NOT on CX3 yet — scp up to `$HOME/ssign-benchmark-20/`). Accessions: NC_003197.2 NC_004337.2 NC_003131.1 NC_008791.1 NC_004578.1 NC_011601.1 NC_013716.1 NC_000117.1 NC_003295.1 NC_003902.1 NC_002516.2 NC_002505.1 NC_004603.1 NC_008570.1 NC_002942.5 NC_003063.2 BX470248 NC_006351.1 NC_014500.1 AE002098. Use `--walltime 24:00:00`.
WHAT TO VERIFY: (a) new figures 01-07 + pooled P01-P03 render; (b) T3SS substrate count is SANE
(not 1,808-flagellar-blowup); (c) enrichment stats have a `T3SS / DLP / window` row and NO
`T3SS / DSE` row; (d) pooled cross-genome figures at 20-genome scale.

### RESULT 2026-06-23 (runs retrieved to validation_sweeps/cx3_t3ss/, untracked):
PASS on the headline. All 3 runs completed (1-genome Salmonella 24/24; 4-genome all 17/17 + pooled;
20-genome all 20 x 17/17 + pooled). T3SS now DETECTED with SANE counts (0-12/genome, 61 pooled over
20 — NOT the 1,808 blowup). Every `T3SS / DSE` row is empty/n.a. (DSE correctly excluded); `T3SS / DLP /
window` rows present with strong enrichment. New figures 01-07 + P01-P03 + enrichment grids all render;
fig 04 SignalP-by-type reproduces Family A/B biology on real Salmonella (T1/T3/T6 = 0% SP+, T5a 7/8,
T5c 1/1). t3ss tasks 5.1/5.2 + #68 satisfied.

### BUG FIXED 2026-06-23 (multi-genome enrichment background):
Root cause: `multi_runner._pool_segment_b_inputs` pooled `neighborhood_proteins` UNCONDITIONALLY and
staged it as both `neighborhood_proteins` AND `proteins`, so the pool runner's DLP/DSE (which read
`proteins` via `_resolve_step_input_fasta(whole_genome=True)`) ran on neighborhood-only even though
`--enrichment-stats` forces `dlp_whole_genome`. Fix: when `dlp_whole_genome`/`dse_whole_genome`, also
pool each genome's whole proteome (`proteins`) into `pooled_whole_genome.faa` and stage THAT as
`proteins`; SignalP/PLM-E keep the neighborhood pool. Now batched DLP/DSE see the full proteome, the
background is genome-wide, and per-genome enrichment matches single-genome. Regression tests:
test_multi_runner.py::TestSegmentBWholeGenomePooling (2). 1379 unit tests pass. TRIGGER TO CONFIRM ON
REAL DATA: re-run the 4-genome CX3 batch with --enrichment-stats; Salmonella T3SS null_mean should now
be ~0.55 (matching the single-genome run), not 0.08. Original report below.

### BUG (original report, now FIXED above):
For the SAME genome (Salmonella), single-genome vs 4-genome runs give IDENTICAL observed + n_mask
(T3SS 4/22, T6SS 2/18) but null_mean differs ~7x (T3SS 0.551 vs 0.08; T6SS 0.451 vs 0.066). Implied
genome-wide DLP positive rate: single ~0.025/pos (~113 positives, biologically plausible ~2.5%),
multi ~0.0036/pos (~16, implausibly low). So MultiGenomeRunner's enrichment background is computed
over too few whole-genome positives, INFLATING fold + significance in batched runs (T3SS fold 7.3x->50x,
n.s.->***). Substrate CALLS are correct; only the enrichment null baseline is wrong. This is the
"segment-B prediction pooling x forced-whole-genome-DLP" interaction flagged at submit time. Trigger to
revisit: before trusting any multi-genome enrichment p/q. Likely in core/multi_runner.py segment-B pooling
or enrichment_testing background construction; needs the full run dir (raw deeplocpro.tsv per genome) on
CX3 to confirm the positive-count gap. Single-genome enrichment is correct and unaffected.
COSMETIC: figures label SS variants raw (T6SSi, a T4aP column) while enrichment collapses via
display_type (T6SS) — consider routing figure ss_type labels through display_type too.

## DeepLocPro crashes on mega-proteins — FIXED 2026-06-18 (openspec deeplocpro-mega-protein-guard)

SHIPPED: DeepLocPro wrapper now withholds sequences > DEEPLOCPRO_MAX_AA (5000, env-overridable)
from the model and emits them as "Not predicted (too long)" rows; the run no longer crashes.
Remaining: rerun BX470251 on CX3 to confirm 24/24 → 67/67, then archive the change. Original report:


Fleet genome **BX470251** (Photorhabdus laumondii TTO1, 4683 proteins) failed: DeepLocPro
exited 1 (GPU OOM, deterministic across 2 nodes), and as a CORE step it cascaded the whole
genome to failure (8/24 steps). Cause: **plu2670 is 16,367 aa** (next longest 5,457) — a giant
Tc/Mcf toxin / NRPS megasynthase. DeepSecE + SignalP handled the genome fine; only DeepLocPro
died. `run_deeplocpro.py` local path (`deeplocpro -f .. -o .. -g negative -d cuda`) has NO
length guard (the 500 cap is DTU sequence-COUNT only). This is a general defect: any toxin/
secondary-metabolite-rich genome (Photorhabdus, Xenorhabdus, giant adhesins) carries mega-
proteins that will kill a run. Matters for the publication + zero-maintenance/longevity pitch.

**Fix (a small change, needs /opsx:propose):** in the DeepLocPro wrapper, set aside sequences
over a safe length (~5000 aa) before invoking DeepLocPro, mark them unpredicted (or default
localization) in the output with a warning, so a single mega-protein can't crash a core step.
Consider the same guard for the other PLM predictors (DSE/SignalP/PLM-E) defensively.
Trigger: implement to get BX470251 to 67/67 and future-proof toxin-rich genomes. Note the
fleet output currently lives on $EPHEMERAL (home FS was ENOSPC at 9% quota — separate RCS issue).

## enrichment-stats validation + PLM-E over-prediction (2026-06-17)

Findings on the PAO1 smoke run (job 3013556), analysis in
`validation_sweeps/benchmark/analysis/enrichment_validation/` (scripts 01/02, figs 01-05):
- **Background bug**: the 200-protein null undershoots the true non-neighborhood background
  (DLP 0.5% vs ~1.3-1.65%, DSE 1.0% vs ~1.4-1.7%), inflating significance. Null-size sweep: 200
  over-calls, 1000 ≈ converged to "all". One false DSE call removed going 200->1000.
- **PLM-E over-predicts massively**: 25.3% of the PAO1 proteome called effector at native
  threshold, 18% even gated at max_prob>=0.8. T6SE-dominated (loosest threshold 0.5), T2SE
  essentially never (their weak type), 36% multi-type. Per-system enrichment: only 2/18 systems
  significant (both weak/spurious), and the real T3SS is DEPLETED (2/25). PLM-E adds no reliable
  enrichment signal. Paper (Zheng 2026 bbag143) reports specificity only on ~150 curated negatives,
  no genome-scale FPR test; recall-tuned thresholds + OR-of-5-ensembles guarantee inflation at scale.

DECISIONS (Teo, 2026-06-17), to implement via a new OpenSpec change (no active change covers this; #70):
1. n_null default 200 -> 1000; use ALL non-neighborhood proteins for the background when predictors
   ran whole-genome (free, exact).
2. PLM-E positivity gated at max_prob >= 0.8 everywhere it's a binary call (enrichment +
   cross_validate), consistent with DLP/DSE.
3. PLM-E OFF by default entirely (Teo's call after seeing the 25%/18% over-prediction).
4. Drop PLM-E from the enrichment test (subsumed by 3, but keep explicit if PLME is ever re-enabled).
Trigger: run `/opsx:propose` for these, pending the deeper PLM-E paper research (2 agents running).

Pre-fleet-launch: CX3 checkout needs `git pull` (job 3013556 predated the PLME-enrichment wiring
a7afbd9 and the SP_WHOLE_GENOME opt-in).

## effector-recovery-benchmark — CITATION-INTEGRITY: NEXT TASK (2026-06-12)

**CONFIRMED real** (not an agent hallucination): 3/3 spot-checks via CrossRef API (deterministic)
match the audit agents exactly — celA/plaA DOI = a soil-16S-ecology paper, VirA DOI = a mouse
tooth-development paper, Tle4 DOI = a protein-design paper. Genes are REAL (real loci/genomes);
it's the sourcing-DOI metadata that's wrong, + 2 genuine SS-type mislabels (BopA/BopE = Bsa T3SS
not T6SS). Our own `doi_resolves` col is blank for all 19 and `verification_status`=VERIFIED, so
"VERIFIED" never meant "DOI points to the right paper."

**CAVEAT — don't extrapolate:** the 19 audited are the most-suspicious slice (emitted-but-machinery-
far), NOT a random sample. 14/19 defective there says nothing reliable about the other ~480.
Prevalence across the full set is unknown → must be measured.

**Precision blind spot:** ssign emitted 1,933 (default) / 2,321 (t3ss) secreted proteins panel-wide;
only 39/51 are gold effectors → **1,894 / 2,270 emissions are unvalidated** (novel TP or FP, no
ground truth). Benchmark measures recall only, says nothing about precision.

**AGREED NEXT STEP (Teo, 2026-06-12):** deterministic citation-consistency sweep (CrossRef
title+abstract vs the row's gene/organism/SS-type; flag mismatches; NO agents in pass 1), then
targeted re-audit of flagged rows. **SCOPE: only the ssign-FOUND effectors for now** (the ~51
emitted_secreted gold effectors in `actual_per_effector.panel_genbank_t3ss.tsv`), to save time —
not the full 582. Items 3 (triage wrong-DOI-vs-mislabel by training impact) + 4 (precision estimate
for the 1.9k unlabeled emissions) deferred to discussion when we get there.

**PASS-1 DONE (2026-06-12)** — secretion-classifier-dataset task 6.1, `scripts/41_citation_consistency.py`
→ `data/dataset/citation_consistency_found.tsv`. **A join bug was caught and fixed mid-task** (the first
reported numbers were partly wrong): joining found→positives on raw UniProt collapsed the 17 found
effectors with `uniprot='-'` (a truthy dash) onto one arbitrary positives row (12 shared DOI nature11433)
→ garbage organism/DOI for those rows. Fixed by bridging found→`ceiling_per_effector` by (effector_locus,
gene) → gold `instance_id`, then keying positives by (instance_id, gene, ss). All 51 now resolve via
instance, 35 distinct DOIs (`join_method` column records how).

**Corrected result (51 found):** **21 CONSISTENT** (14 by gene name, 7 genus-only = weaker),
**8 FLAG_WRONG_TOPIC** (DOI in an unrelated field — celA/plaA→soil ecology, Map×2→Toxoplasma,
VirA→tooth development, BipB→colicins, aprA→CS1 pilin, and **EspZ→the T6SS-discovery PNAS paper**),
**3 FLAG_GENE_ABSENT** (prtB/prtA1/aprA cite a real Prt/Lip ABC-exporter T1SS paper that just doesn't
name the protease in its abstract — probably OK, confirm in pass-2), **3 DOI_UNRESOLVED** (YspE/CopN/
ChlaDub1 — classic ASM/JBC/Science DOIs failing the DOI.org handle check; likely legacy false-negatives
or typos), **16 INDETERMINATE** (CrossRef record exists but no abstract → can't refute deterministically).
Net: the deterministic method adjudicates ~30/51; 16 indeterminate + 14 flagged = 30 go to pass-2.
Validated vs the 19 prior manual verdicts: re-flagged celA/plaA/VirA independently AND caught **BipB's
wrong DOI the manual audit missed** (it was "sound" only because its agent was AUP-blocked; the DOI is a
colicin paper). Also: **aprA has two instances with two different DOIs** (one wrong CS1-pilin, one OK-ish
lipase) → the defect is per-row metadata, not per-protein. Caveat: ssign-found ~10% slice; full-582
prevalence still unmeasured.

**PASS-2 DONE (task 6.2)** — `pass2_results.tsv` + `pass2_results.md`; raw returns `pass2_raw/batch_*.json`;
verify `scripts/42_pass2_verify.py`. 30 non-CONSISTENT rows re-audited: 5 literature agents +
batch C (YspE/VirA/BipB/CopN/ChlaDub1) done in-session via PubMed (agent AUP-blocked on B. pseudomallei).
Every returned DOI deterministically re-verified (registered + on CrossRef + gene/genus present):
**26 rows now carry a verified-real on-topic primary DOI** (21 RESOLVED corrected, 3 CONFIRMED,
2 MISASSIGNED), **0 fabricated/unverified DOIs** — the agents agreed with the deterministic check on all
26, so they were reliable here.

**Dataset edits this produced (for the training set, NOT yet applied):**
- **2 MISASSIGNED ss_type:** BopA + BopE (B. thailandensis BTH_II rows) are **Bsa T3SS** effectors, not
  T6SS (name-collision; BopA=Cullinane 2008, BopE=Stevens 2003). Fix ss_type or drop. These are 2 of the
  benchmark's 6 "accidental cross-type" emissions.
- **4 NOT_FOUND (drop/down-tier):** ChlaDub1 (no T3SS-secretion evidence, DUB-only per Misaghi 2006);
  P.entomophila aprA (homology-only, no ortholog-specific T1SS paper); EFF00136 + EFF00150 (BTH DUF3274
  proteins not named in their cited paper).
- **1 duplicate:** Tle4 (idx15) = TplE_alias_Tle4 (idx24) = PA1510; BOTH also had wrong DOIs (protein-
  design / lung-regeneration papers) → correct = Jiang 2016 Cell Rep 10.1016/j.celrep.2016.07.012. De-dupe.
- EFF00142's pass-1 "unidentifiable" tag was WRONG: Russell 2012 names BTH_I2691 as a T6SS-1 substrate → CONFIRMED.
- EspZ's DOI was the T6SS-discovery PNAS paper → corrected to Kanack 2005 (10.1128/IAI.73.7.4327-4337.2005).

**ITEM 3 DONE (2026-06-12, task 6.3)** — applied to `positives_all.tsv` via `scripts/43_apply_citation_corrections.py`
(post-build overlay, **930→925 rows**). Triaged by training impact: **22 pure-provenance DOI fixes** (label
unchanged), **2 label-changing** (BopA/BopE T6SS→T3SS, re-labelled to type-level T3SS + DOI fixed), **5 removed**
(4 unsupported + TplE_alias_Tle4 duplicate → `positives_removed_citation.tsv`). Backup =
`positives_all.pre_citation_audit.tsv`; per-field log = `citation_corrections_log.tsv`. Verified surgical
(ss_type counts moved exactly T1−1/T3+1/T6−5; other instances of EspZ/Tle4 untouched), idempotent, and
byte-identical on re-run from the backup.

**Load-bearing caveats:**
- **Overlay, not in the build chain:** re-running `37_fold_t5ss.py` regenerates the uncorrected 930-row table,
  so `43` MUST be re-applied after any `37` rebuild. Documented in `data/dataset/README.md` "Citation-audit overlay".
- **Downstream re-run:** `40_pair_features.py` (and the eventual feature matrix) must be rebuilt — 2 ss_type
  changes + 5 dropped rows. pair_features.tsv is now stale.
- **Out of scope (real, deferred):** the audit only touched the ssign-FOUND rows. Other non-found instances may
  carry the SAME wrong DOIs — confirmed: a second EspZ instance (T3SS_24) still cites the T6SS-discovery PNAS
  paper. Fixing those is the full-corpus citation sweep (the original "prevalence unknown" item), still deferred.

**ITEM 4 DONE (2026-06-12, effector-recovery-benchmark §8)** — deterministic precision estimate, scripts
28→29→30→31, `data/phase2/PRECISION.md` + `figures/precision/01..03`. Teo chose deterministic tiers only
(DB-confirmed floor + obvious-FP), agent-sampled adjudication deferred. Denominator = the 1,572 proximity
substrate calls (default); 361 T5SS-self assessed separately (correct-by-construction; 68% effector /
0.3% housekeeping → T5SS detection sound). **Result: proximity precision is a wide band, ~3% provable
floor (DB-confirmed homology to SecReT4/6; T6SS the only well-covered type at 8.3%) → ~75% soft ceiling
(not-obviously-non-secreted), with ~19% clearly-FP housekeeping, ~6% machinery, ~11% annotation-effector,
and ~64% (hypothetical+other) unresolvable by DB or annotation.** Pairs with recall (8–10% emitted of
testable): the proximity rule is both permissive AND low-recall → the case for the classifier, which can
adjudicate the unresolvable middle. Tier-1 used pyhmmer phmmer (no external aligner; the repo venv is at
`../../.venv`, NOT `benchmark/.venv` — stdlib scripts silently fell back to system python3).

**Caveats baked into PRECISION.md:** no negative ground truth (ceiling overstates); annotation-based tier
inherits genome-annotation errors (NleB mis-annotated "IS3 transposase"); DB floor covers only T4/T6SS.

**STILL OPEN:** the full-corpus citation sweep (point 41→42→43 at the unaudited rows; defect confirmed
beyond the found slice), and optionally a stratified agent-adjudicated precision sample to place a point
estimate inside the 3–75% band. Both deferred to discussion with Teo.

### Earlier deferred repair items (from the 19-row audit; full table in `data/phase2/discordant_audit.md`):
- **Wrong/nonexistent sourcing DOI (≥9):** celA, plaA, VirA, ChlaDub1(404), CopN(404), Tle4, TplE, Tle1, Tae4_Stm. Re-source before trusting.
- **Genuine misassignment:** BopA, BopE labelled T6SS but are Bsa **T3SS** effectors (Bop name-collision) → reassign. ChlaDub1 T3SS call unsupported (inclusion-membrane DUB) → review.
- **Unidentifiable/unsupported rows:** EFF00142, EFF00150 (opaque DB placeholders), TseA_T6SS1 (cited paper never names it) → drop or repair.
- **Duplicate row:** Tle4 == TplE_alias_Tle4 (PA1510/Q9I3K2) → de-dupe.
- Same defect class likely across the wider corpus → run citation audit on ALL positives, not just these 19.
- BipB/BipC agents blocked by AUP filter on B. pseudomallei content; checked manually (sound Bsa T3SS translocon). If a future agent audit is needed, run those two by hand.

## effector-recovery-benchmark — FULL PANEL DONE (2026-06-12)

12 jobs (2988054-2988065) finished clean on CX3 gpu72: both tags 67/67 genomes, 0 empty/failed, 196,170 protein rows each. Synced to `data/phase2/runs/panel_genbank_{default,t3ss}/`; scored with scripts 24 + 25. Outputs `data/phase2/actual_{per_effector,vs_ceiling}.panel_genbank_{default,t3ss}.tsv`.

**Exact ssign-detected systems behind the found effectors** (`scripts/26_found_systems.py`, `found_systems.*.tsv`): default 39 found effectors via **32** distinct detected systems; t3ss 51 via **39**. Emission basis (why emitted): default own-type(legit) **33/39**, cross-type-only(accidental) **6/39**; t3ss own-type **45/51**, cross-type-only **6/51**. So the large majority of found effectors sit next to their OWN system type; only **6** are accidental cross-type emissions (the celA/plaA-via-T4aP, BopA/BopE-via-T5cSS, VirA/ChlaDub1-via-T5aSS set, same as the audit's cross-type rows). NB: required a family-normalization fix in script 26 (ssign labels detected systems by TXSScan subtype `T6SSi`, not the coarse `T6SS`); an earlier uncorrected run wrongly showed ~half accidental.

**Full-panel recall (emitted of testable, testable=499):**

| SS | testable | ceiling@7 | actual default | actual t3ss |
|----|---|---|---|---|
| T1SS | 25 | 80% | 60% | 60% |
| T2SS | 77 | 4% | 5% | 5% |
| T3SS | 227 | 33% | 1% | 7% |
| T4SS | 98 | 10% | 0% | 0% |
| T6SS | 72 | 26% | 24% | 24% |
| ALL | 499 | 25% | **8%** | **10%** |

Headline: even within +/-7 genes (ceiling 25%), as-shipped ssign emits 8% (T3SS-excluded) / 10% (T3SS-included). The ceiling-to-actual gap (25% -> 8-10%) is the core motivation for the secretion-classifier. T4SS actual=0% is the starkest gap (ceiling 10%, emits nothing). T3SS default=1% is by-design (excluded); the t3ss pass lifts it to 7%, still far under the 33% ceiling. 16-19 effectors emitted but ceiling-unreachable@7 = ssign's detected machinery disagrees with the literature answer key (mostly T6SSi self-adjacency, a few T5aSS/T4aP); benign, listed in script-25 output.

Denominator honesty: 14 no_run (effectors with no genome accession, all non-testable) + 75 not_in_input (only 6 testable; ssign's ORF set didn't carry those 6 loci = forced misses already counted against recall). Nothing inflated.

**T3SS detection CLARIFIED (2026-06-12), corrects the CLAUDE.md "MacSyFinder found 0 T3SS" framing for THIS panel.** `excluded_systems` (default `Flagellum,Tad,T3SS`) is a DOWNSTREAM FILTER, not a detection switch (`system_filtering.py:39`, `validate_macsyfinder_systems.py:148` separates included vs excluded AFTER parsing all MacSyFinder results). MacSyFinder/TXSScan models T3SS and runs it every time. In the t3ss-included panel run it **detected 30 real T3SS systems** (Pseudomonas/Bordetella/Chlamydia/Salmonella/Shigella, `excluded=False`); the default run shows 0 only because the filter drops them. So T3SS is OFF by default NOT because MacSyFinder can't find injectisomes (it can) but because DeepSecE over-predicts T3SS by misclassifying flagellar proteins (CLAUDE.md bug #4) — excluding it protects precision. Even with T3SS on, found only 3→15 because ~73% of T3SS effectors are genome-dispersed (unreachable @±3). The CLAUDE.md "0 across 74 genomes" note is about an older/different dev set; don't treat it as true for the benchmark panel.

**"Non-testable" (83) ≠ "no genome" — it's a mix** (corrects the figure-01 shorthand): effector_locus_not_found 26, own_instance_unknown 25 (the dropped net-new T6SS multi-instance effectors), no_genome 11+2 divergent, machinery_unanchored 10, no_instance_in_genome 8. All mean "couldn't fairly put it in front of ssign" (no assembly, ORF absent from the staged assembly, or no anchorable machinery to measure ±3 against), so excluded from the found/missed judgement.

**UNREACHABLE-MISSED ANALYSIS DONE (2026-06-12)** — `data/phase2/UNREACHABLE_ANALYSIS.md` + `figures/summary/04_distance_to_machinery.png` (script 33) + 5 literature agents. **Distances to machinery (testable): T1SS median 1 gene (20/23 within ±3 = operonic, Teo's intuition correct); T2SS 302, T3SS 45, T4SS 203, T6SS 232.** So T2/T3/T4SS effectors are GENUINELY genome-dispersed → high unreachable@3 is EXPECTED biology, not a ssign/benchmark failure (T2 substrates recruited post-translationally in periplasm; T3 effectors horizontally acquired on prophages/islands, shared-regulon-coordinated; T4 effectors recognized by portable C-terminal signal, location decoupled — Legionella/Coxiella ~300 scattered). The 5 T1SS misses are all genuine biology (Serralysin=LipBCD generalist exporter at separate locus; ApxIIA=in-trans via apxI-BD; FrpC=functional but scattered T1SS, no adjacent TolC/HlyB to detect; TRP47/32=T1SS but Hly transporter at separate locus). Core argument for the classifier: a fixed proximity window is the wrong abstraction for the dispersed systems.

**On the frpC/TolC question (Teo):** ssign detected NO T1SS in the Neisseria genome and there's no TolC/HlyB/HlyD adjacent to frpC — the apparatus is real but genome-scattered, so nothing to detect nearby. The answer-key "TolC at 1340 genes" is a spurious product-tier match.

**3 ACTIONABLE FOLLOW-UPS surfaced:**
1. **[DONE 2026-06-12]** TRP47/TRP32 (Ehrlichia, T1SS) were anchored on **VirB8 (a T4SS gene)**; +3. **[DONE]** frpC anchored on a spurious TolC 1340 genes away. `scripts/34_answer_key_corrections.py` reclassifies all 3 → `machinery_unanchored` / non-testable (correct machinery not anchorable; checked: ssign detects NO secretion system in either the Ehrlichia or Neisseria genome, and no Hly/TolC is adjacent — ssign could not have found them regardless). Found count unchanged (all not_emitted); testable 499→496; T1SS unreachable 5→2 (now just the genuine-biology Serralysin + apxIIA). Backups `*.pre_anchor_fix`. Figures 01/04 regenerated.
**FALSE-NEGATIVE DIAGNOSIS DONE (2026-06-12, task 8.7)** — `FALSE_NEGATIVES.md` + `figures/summary/05` (script 35). Of the 62 reachable@3-but-missed effectors, **50 (81%) are DETECTION failures**: ssign detected no secretion system, and the secreted-protein predictors run ONLY on the ±3 neighborhood of a detected system (`runner.py` dlp_whole_genome=False), so the effector was never evaluated (all tool signals blank). Per type: T1SS 5/5 detection, T3SS 37/40, T4SS 6/6, T6SS 2/11 (T6SS is the exception: 9/11 were processed-but-filter-rejected — the genuine prediction-side misses). **T1SS clincher:** all 5 are RTX toxins; HlyA in a COMPLETE genome (NZ_CP031766.1) has the full hlyC-A-B-D operon annotated (transporter 1-2 genes away) yet ssign called no T1SS. Likely cause: TXSScan T1SS model needs a co-localized TolC, but TolC is a distant housekeeping gene, so no complete system assembles. **Recall is bottlenecked by SS DETECTION, not the predictors.** Two takeaways: (1) ssign-side fix — T1SS detection should tolerate non-co-localized TolC / accept HlyB+HlyD operons; (2) a per-protein learned classifier (not gated on detection+proximity) would recover these obvious toxins = the recall-side classifier argument. **Open drill-down:** the 9 T6SS processed-but-rejected cases (why did cross-validation reject them?).

2. **[STILL OPEN] T6SS ceiling likely UNDERCOUNTED:** we anchor on the core tss cluster (TssM/ClpV), but many T6SS effectors sit beside an orphan vgrG/paar/hcp far from the core, with cognate immunity downstream. Add a nearest-vgrG/paar + immunity-pair anchor and recompute T6SS reachable@3 (likely rises; ssign detects vgrG so some "unreachable" T6SS may be reachable). The one system where apparent proximity-failure is partly a benchmark-anchor choice, not biology.

**Next (now unblocked):** 6.6/6.7 figures + docs; dataset group-4 4.1 feature join (the `actual_per_effector.*.tsv` tables ARE the per-protein tool signals). Benchmark-side bridge/SUBMIT/script changes still UNCOMMITTED (not yet approved).

## effector-recovery-benchmark — Phase 2 pilot results + FASTA bridge bug (2026-06-11)

4-genome pilot (NC_002516.2/003197.2/004337.2/004578.1) ran on CX3, 3 tags x 4 genomes, all 20/20 steps. Synced to `data/phase2/runs/`; scored with scripts 24 + 25. Outputs: `data/phase2/actual_{per_effector,vs_ceiling}.<tag>.tsv`.

- **In-panel pilot numbers (genbank_default):** of 118 testable effectors in the 4 genomes, ceiling-reachable@7 = 30 (25%), ssign emitted 4 = **3% of testable / 13% of proximity-reachable@7**. T3SS-enabled tag: 6 emitted = 5% / 20%. The gap is large and real (supports the classifier; doc 00 decision tree "gap >20% -> proceed"). NB script 25's printed table dilutes "actual" across all 499 effectors incl. 445 no_run (not-in-this-4-genome-pilot) -> shows 1%; the in-panel 3-5% is the honest pilot number. **Full panel run needed for real figures (6.6/6.7).**
- **FASTA-mode bridge — FIXED 2026-06-11 (coordinate bridge).** Was 0/137 matched: Bakta renames locus_tags (`NFOBNJ_00001`) and the seq fallback read a `sequence` column ssign's results_raw never emits. Fix: `bench_runout.RunOutput.find_by_coord` + `by_coord` index match on (contig_base, strand, 3'-stop) since Bakta calls the same ORF at the same RefSeq coordinate; script 24 looks up effector coords from the gene-order index. Result: **0→116 matched, and the emitted set is identical to genbank mode** (4 effectors, exact), confirming correctness. Remaining 21 unmatched are effectors with no locus_tag at all (unbridgeable). If ssign ever emits an aa `sequence` column, the seq bridge reactivates automatically (kept as a later tier).
- **3 emitted-but-ceiling-unreachable@7** (VirA/Q7BU69 T3SS nearby=T5aSS; Tle4 + TplE_alias_Tle4 T6SS nearby=T6SSi): ssign emitted these off a *different* nearby system than the literature answer-key assigns. This is the benchmark task-6b "ssign machinery != literature instance" signal; capture in the 6b write-up, not a bug.

## CX3 environment

- **signalp6 + deeplocpro PATH on fresh nodes** (task #22). Now auto-discovered via `_find_in_conda_envs` in `core/runner.py`. Binaries execute directly without `conda activate`, which works for pip-installed Python entry points because the shebang pins the env's Python. **Open risk:** if either tool's torch build needs system CUDA libs (rather than its own bundled CUDA wheels), `LD_LIBRARY_PATH` won't include `<env>/lib/` and the subprocess will crash with `libcudnn.so.X not found`. Mitigation if it ever fires: prepend `<env>/lib` to `LD_LIBRARY_PATH` in `run_signalp.py` / `run_deeplocpro.py` when the binary path lives under a conda env. Trigger to revisit: a CUDA-related ImportError from either wrapper.

## PLM-Effector performance

- **Ensemble model checkpoint cache** (task #16). `run_ensemble` now called 85× instead of 5×; re-reads ~150 MB of small files per call. Page cache makes real cost ~150 MB total, so lower priority than the raw number suggests. Revisit: if PLM-E step wallclock is still the long pole after `0445d94` cross-type caching is validated.
- **FP16/BF16 + `--batch-size 32`**. Could drop PLM-E ~74m → ~12-15m on whole-genome runs. Needs validation that FP16 doesn't shift predictions. Now low priority: PLM-E runs on the SS neighborhood (~128 proteins) by default, so absolute wallclock is small. Only revisit if `--plme-whole-genome` becomes a common workflow.

## Disk sizes (measured 2026-06-03 on CX3, scripts/audit_disk_sizes.py)

- **base 2 GB / extended 140 GB / full 1.3 TB**. BLAST nr is the long pole at 802 GB; users without nr-cross-genome BLASTp save 1.2 TB by staying on extended. UniRef30 dominates HH-suite at 261 GB. Several install.md per-tool estimates were ~50% off (EggNOG 25→47, IPS 24→35, PLM-E 18→26, HH-suite 55→340). Updated in docs/how-to/install.md.

## Torch.load safety

- **Migrate `run_deepsece.py` to `weights_only=True` + `add_safe_globals`.** DeepSecE's checkpoint is a state_dict (already loaded via `model.load_state_dict(...)`), not a whole-module pickle like PLM-Effector's, so it's a safe candidate for the stricter loader PyTorch 2.6 introduced. PLM-E itself can't migrate (whole-module saves; would need an upstream refactor). Revisit: any time a new ssign dep bumps torch and `weights_only=False` triggers a deprecation warning.

## Annotation parallel-group scheduling

- **Wave-scheduling or finish-and-redistribute** (deferred 2026-06-03). After fixing the 3-4x oversubscription bug (commit `82cece9`), each annotation tool now gets `effective_cpu_count / N` cores in the parallel group. Remaining inefficiency: when one tool (e.g. IPS at 36m) finishes before the others (EggNOG, pLM-BLAST), the cores it released sit idle because the surviving tools were started with a fixed thread count and can't add workers mid-run. Three options considered: (1) two-wave scheduling (fast tools then slow tools), (2) restart-survivors-with-more-threads (throws away expensive embed/diamond work), (3) leave as-is. Going with (3) for now. **Trigger to revisit:** once we have step_timings.csv from a few post-fix runs across different genome sizes — if there's a consistent 10-20%+ wallclock floor that's clearly tail-tool dominated, the wave-scheduling experiment becomes worth it.

## Benchmark: T1SS effector rescue (effector-recovery-benchmark, 2026-06-11)

The corpus left 19 of 28 validated T1SS effectors without a genome (old UniProt entries cite 1980s-90s EMBL gene clones, not assemblies), so they were untestable for the proximity ceiling. Per Teo's decision (representative-species placement; floor set to >=90% id after he asked to include the 93% prtA), `scripts/12_rescue_t1ss_ipg.py` (IPG, identical protein) + `13_rescue_t1ss_blast.py` (remote blastp, >=90% id / >=90% cov) + `14_finalize_t1ss_rescue.py` placed **16/19** into RefSeq genomes (9 IPG-identical, 7 BLAST; 13 exact-species, 3 genus-only). Output: `data/t1ss_rescue/t1ss_rescued.tsv`. (prtA Q07295 = 93.1% in same-species D. chrysanthemi SR64; the rest 98-100%.)

- **3 still unplaced, with reasons (genuinely unplaceable, not a decision):**
  - `prtG` (Q07162, Dickeya chrysanthemi): best genome match only **61%** id. Genuinely divergent.
  - `prtA` (P82115, Photorhabdus sp. Az29) + `lktA` (P55123, "Pasteurella haemolytica-like sp. 5943B"): "no significant similarity" because those exact strains have **no RefSeq genome**; the search was correctly species-restricted. *A genus-level fallback (flagged genus_only) could rescue the conserved leukotoxin; the Photorhabdus one is a partial UniProt entry and riskier.* Trigger to revisit: if Teo later wants >=18/19, add a genus-level BLAST fallback.
- **Adjacency verified (script 15 + 16, 2026-06-11):** instead of assuming ceiling=100%, we read the gene order around each placed effector for the T1SS transporter (ABC/HlyB-family + HlyD-family MFP; literature families off RefSeq annotation, ssign-independent). **14/16 CONFIRMED adjacent** (transporter 1-3 genes away -> reachable at N=3). **2 genuine exceptions:** apxIIA (apxIICA operon lacks a transporter; trans-secreted by ApxIB/D from the distant apxI operon) + serralysin (Lip/LipBCD transporter at a separate locus, confirmed across 13 complete Serratia genomes). These 2 are real "impossible" cases. `16_t1ss_replace_fragmented.py` fixed hlyA (was on a 3 kb single-CDS WGS contig that truncated hlyCABD -> re-placed into complete genome NZ_CP031766.1, operon intact). Tables: `data/t1ss_rescue/t1ss_ceiling.tsv`.
- **Phase 1 integration still TODO:** fold the 16 placements into `effector_gold_set.tsv` (add `placement_tier` ipg_identical/representative_strain + `species_match` columns), add the placement genomes to the Phase 2 panel, and record the T1SS ceiling as 14 reachable + 2 impossible. Placement genomes already cached in `refseq_cache/`. Trigger: start of Phase 1.
- 3 placements are sister-species (genus_only): `prtA1` P.luminescens->P.akhurstii, `prtB`/`prtC` D.chrysanthemi->D.dadantii. Sequence is 100% identical (IPG), so defensible, but flagged.

## Benchmark: Phase 1 ceiling (effector-recovery-benchmark, 2026-06-11)

Phase 1 complete: ceiling = % of testable verified effectors the +/-N proximity rule could reach. **T1 80% / T2 1-4% / T3 21-33% / T4 6-10% / T6 22-26%; 499 testable, 83 untestable.** Scripts 17-21 + `scripts/bench_index.py`; outputs in `data/phase1/`, figures `figures/01-04`, writeup `docs/phase1_ceiling.md`. Deferred / judgment-call items:

- **33 net-new SecReT6 effectors dropped to untestable** (Teo's Checkpoint-A-consistent call): 25 in multi-T6SS genomes + 8 in Serratia SMDB11 with no curated T6SS instance. We do NOT guess the nearest instance (circular). T6SS ceiling rests on the 72 rigorously-assignable effectors. **Trigger:** if Phase 2 makes T6SS coverage worth expanding, literature-assign each net-new effector to its specific T6SS (SecReT6 citations), then re-run 19.
- **pulA (T2SS_16) untestable.** Its instance was dropped at Checkpoint A (X12831 = 1.8 kb EMBL pul fragment, not a genome); not backfilled. The founding T2SS effector is therefore not in the ceiling. **Trigger:** if wanted, place pulA into a complete K. pneumoniae genome + curate the gsp machinery there (real Phase 0b work).
- **26 effector_locus_not_found** = corpus locus_tag scheme absent from the available RefSeq assembly (different assembly/annotation), no unique gene symbol to bridge. Genuinely unplaceable in cache. **Trigger:** coordinate-based placement via UniProt->genome if completeness matters.
- **Gene-symbol fallback (9 effectors, incl. canonical Yersinia Yops on pYV).** Located by *unique* /gene symbol when the tag scheme was missing; flagged `effector_match=gene_symbol` in `ceiling_per_effector.tsv`. Uniqueness-gated (no paralog ambiguity), but a weaker identifier than locus_tag. Documented in `docs/phase1_ceiling.md` strengths.
- **off-replicon machinery -> impossible.** Effectors whose own-instance machinery anchored only to a different replicon are counted impossible (structurally unreachable by a +/-N window), not untestable. Reason `machinery_off_replicon` in the per-effector table.

## ssign output: emit genomic coordinates (base pipeline)

- **`results.csv` and `results_raw.csv` carry no genomic coordinates** (only `locus_tag`, `sequence`, per-tool signals). The intermediate `gene_info` step already has `contig,start,end,strand` (`extract_proteins.py:384`) but `_build_master_csv` / `_build_raw_csv` (`core/runner.py`) drop them. **Add `contig,start,end,strand` to the raw CSV at least** (left-join from gene_info, which is already the raw base). **Why it matters:** any coordinate-based downstream use (operon context, the effector-recovery benchmark's Bakta->RefSeq bridge, the secretion-classifier model's positional features) currently has to recover coordinates by re-parsing the input or matching on sequence. The benchmark Phase 2 worked around it with locus_tag + protein-sequence identity matching, but emitting coordinates removes the workaround and is generally useful. **Trigger:** next base-pipeline pass; small change, gene_info is already the join base in `_build_raw_csv`. (Flagged 2026-06-11 during effector-recovery Phase 2.)

## Statistics

- **(Resolved 2026-06-02)** The broken permutation + biased Fisher path was replaced by the A+ rewrite: opt-in `--enrichment-stats` flag, null sample of N=200 random non-SS-neighborhood proteins per genome, scipy binomial test per real SS system + per broad type, BH FDR. Multi-genome runs also emit a pooled view. See `enrichment_testing.py` + `pool_enrichment_stats` in `core/runner.py`.

## Shared TSV/parsing helpers (#75a simplify follow-up, resolved 2026-06-05)

- **Tolerant int parsing — resolved.** `ssign_lib/parsing.py:parse_int_or_none(value, allow_range=False)` is the single source. `t5_passenger` imports it directly; `t5ss_handler._parse_sp_end` is now a one-liner wrapping it with `allow_range=True` so call sites still read as `_parse_sp_end(...)`.
- **"Load TSV → dict by key" — resolved.** `ssign_lib/tsv_io.py:load_tsv_by_key(path, key_columns, missing_ok=True)` is the single source. `cross_validate_predictions._load_tsv_by_locus` passes the tolerant fallback chain `("locus_tag", "protein_id", "seq_id")`; `enrichment_testing.load_predictions_keyed` keeps strict-locus_tag-only via `key_columns=("locus_tag",)`. `t5_passenger.load_t5_classifications` was left as-is — it does substantial per-field coercion on top of the TSV read, so the shared helper would only save the open-and-iterate boilerplate.

## secretion-classifier-dataset — deferred re-citation (group 1, 2026-06-11)

- **25 predicted positives kept with UNRESOLVED citations** (`evidence_tier=predicted`, `citation_status=UNRESOLVED` in `data/dataset/predicted_audited.tsv`). The audit kept them (protein/locus/instance label is sound; the broken DOI is metadata, not a label, same bar as the validated gold set which also never drops on citation alone). 8 distinct broken DOIs: lktA/mbxA/movA RTX leukotoxins (11), prtA/prtC proteases (4), Ehrlichia T4SS EBP/ECH (3, IAI.00513-13), Shigella OspD2/3 (2, science.1175302), Anaplasma T4SS APH (2, cmi.13405), Chlamydia CT813 (1), V. vulnificus rtxA (1, JB.187.10.3392), Dickeya pelZ (1, annurev-micro-102215). 36 of the original 61 were already recovered (Apx ×32 → Frey 1993 10.1099/00221287-139-8-1723; ApxIVA ×2 → Schaller 1999; V. cholerae rtxA ×2 → Boardman 2004) in `doi_recite.jsonl`. **Trigger:** before any public release that cites the predicted tier, run a per-paper re-citation for these 8 (one paper each, all well-characterized families); add verified DOIs to `doi_recite.jsonl` and re-run scripts 31→32. Held back because applying a shaky umbrella citation to benchmark provenance is worse than an honest flag. The agent re-citation route hit an API policy block on toxin literature, so this was done inline; the remaining 8 need careful manual lookups.

## secretion-classifier-dataset — predicted instance audit flags (group 2, 2026-06-11)

The 2.3 literature-audit agent (27 RESOLVED / 17 UNRESOLVED of 44 ambiguous) flagged 3 corpus data-quality issues to fix upstream (assignment stands at gene/family level; not blocking):
- **C. rodentium ICC168 NleG locus_tags mislabeled.** Corpus has ROD_31501/ROD_21621/ROD_25881/ROD_15971; NCBI annotates those as exported-protein/hisA/prophage/hypothetical. Real ICC168 NleG tags are ROD_16511 / ROD_48891 / ROD_40971. **Trigger:** correct these 4 locus_tags in `T3SS_verified.tsv` before feature-join (else bench_runout won't find them in run output).
- **VopX / VPA1374 name-locus mismatch.** UniProt annotates VPA1374 as "Uncharacterized"; "VopX" is a V. cholerae effector. Left UNRESOLVED. **Trigger:** verify the V. parahaemolyticus locus for VopX (or drop if a mislabel).
- **EC042_4675 (Tle3_Sci2) unresolvable.** Locus returns zero NCBI/UniProt hits and sits ~100 genes outside both E. coli 042 T6SS clusters; documented EAEC 042 Tle is Tle1 (Sci-1). Left UNRESOLVED. **Trigger:** re-source the correct locus from the bioRxiv ref (10.1101/2025.02.11.637775) or drop.

## ssign extended tier is GPU-gated (PLM-Effector cannot fall back to CPU via runner)

- **Finding (2026-06-11, while diagnosing CX3 pilot queue):** `run_plm_effector.py` defaults `--device cuda` and `plm_effector/predict_api.py:_resolve_device` **raises RuntimeError if CUDA is absent** (only `device='cpu'` explicitly avoids it). `core/runner.py:_step_plm_effector` (~line 2377) builds the PLM-E argv WITHOUT a `--device` flag, so it always uses the cuda default. Net: **extended/full tier crashes at the PLM-Effector step on any GPU-free node** (CPU-only HPC nodes, Amine's Mac, laptops). DeepSecE already auto-selects CPU (`run_deepsece.py:_select_device`), so PLM-E is the lone hard GPU gate.
- **Impact:** can't schedule benchmark Phase-2 runs on the abundant CPU queues (must wait for gpu72 GPU). Also contradicts the "runs anywhere / zero-maintenance" longevity pitch and breaks Mac/laptop extended runs.
- **Fix:** add a `plm_effector_device` config knob + env (mirror DeepSecE's `_select_device` auto-detect: cuda if available else cpu) and pass `--device` through in `_step_plm_effector`. Keep cuda the default when a GPU exists (CPU is ~1-2 min/protein per the code's own warning, but the +/-3 neighborhood keeps protein counts small). **Trigger:** next base-pipeline pass, or sooner if the CX3 GPU queue stays multi-hour and we want to run the panel on CPU nodes.

## secretion-classifier-dataset — T5SS sourcing gaps (group 3, 2026-06-11)

T5SS sourced by subtype (scripts 36/37): 29 proposed → **23 placed + verified** (DOI.org-resolved + locus in genome), folded into `positives_all.tsv` (930 total). Subtypes: T5a 7, T5b 11, T5c 3, T5d 1, T5e 1. self_secreted=true for the 12 autotransporters (T5a/c/d/e), false for the 11 T5b TpsA substrates.

Deferred:
- **Autotransporter sourcing was cut short** by a wifi/socket drop after 18 entries (agent was targeting ~25-40). The brief's anchor list has more uncovered T5a/T5c examples (SPATEs Pic/Sat/Vat, IgA proteases, more TAAs: Hia/UspA2/BadA/SadA). **Trigger:** if more T5SS depth is wanted, re-run one autotransporter sourcing agent (it appends to `t5ss_raw/`, then 36→37 re-fold). Low urgency — T5SS is a supplementary/partly-exploratory tier.
- **6 unplaceable** (`t5ss_unplaceable.tsv`): hap (HI_0248 absent from NC_000907.1 — likely wrong tag), tsh (APECO1_O1CoBM73 tag mismatch), hia/badA/uspA1/invA (agent gave no usable locus_tag / no genome accession). **Trigger:** re-source these specific loci if T5SS coverage matters; each needs one targeted NCBI lookup.
- **GCF→nucleotide resolution is NCBI-flaky.** `36.replicons_for` uses assembly→elink→nuccore but NCBI eutils dropped connections (EOF) repeatedly; the 5 Haemophilus T5b were rescued by hardcoding GCF→NC in `tps.json` (NC_000907.1/NC_002940.2/NC_017452.1). If re-sourcing, prefer instructing agents to emit nucleotide (NC_/NZ_) accessions, not GCF.

## secretion-classifier-dataset — model-handoff deferrals (group 5, 2026-06-11)

`positives_all.tsv` (930) is the label-complete deliverable; `data/dataset/README.md` documents schema + conventions. Items below belong to later changes, not this one:

- **validated:predicted weight ratio is unset.** `evidence_tier` is recorded but the loss weight is a model-training decision. **Trigger:** the model-training change (nnPU loss); start from validated=1.0 and sweep predicted ∈ {0.3-0.7}.
- **Type-level positives (298 rows, `type_level=yes`) have no instance, so pair-features are null.** Whether they feed the instance model or only a separate type-level head is undecided. **Trigger:** model architecture design; if instance-only, hold the 298 out or route them to the type-level head.
- **T5SS subtype depth is a first pass** (23 placed: T5a 7 / T5b 11 / T5c 3 / T5d 1 / T5e 1). Depth decision (source more vs leave as supplementary) waits on whether the self-secreted tier helps. **Trigger:** after the first model eval shows whether self_secreted rows move the needle; if yes, re-run autotransporter sourcing (see group-3 note).
- **T5a-neighbor DLP/DSE observation stays exploratory.** The "autotransporter neighbours look secreted" signal lives in benchmark task 6b as a side-study, never a training label. **Trigger:** none for this dataset; revisit only if 6b finds it predictive.
- **Feature side (group 4) is blocked on CX3 Phase-2 panel runs.** `training_dataset.tsv` (positives + DLP/DSE/SignalP/PLM-E + ESM ref + gene-distance + PU unlabeled set) joins `results_raw` via `bench_runout` once runs land. **Trigger:** Phase-2 pilots flip Q→R on CX3 gpu72; then rsync runs back and run group-4 scripts (4.1-4.5, not yet written). The 4.5 assembler's output schema is now **pinned**: it must emit the columns in `secretion-classifier/secretion_classifier/schema.py` (REQUIRED_COLUMNS), which the model's `TrainingDataset` loader validates against. 4.2 pair-features (`pair_features.tsv`) already supplies the pair/system columns.
- **Model-prep core is built** (2026-06-11, separate repo `reidmat/secretion-classifier`, commits f641d32..3b04d54): data-independent loader/losses/prior/model/splits/metrics/sweep, 33 tests green. Only the trainer + ESM extraction + `training_dataset.tsv` remain, all gated on the same CX3 runs. See memory [[project_secretion_classifier_model]].
- **Label sentinels are a bare-string cross-file contract.** `evidence_tier`, `instance_source`, `type_level`, `self_secreted` values are bare literals shared across scripts 36/37 and (soon) group 4 (e.g. `self_secreted == "true"` in 37 depends on 36 writing `str(...).lower()`). A typo mislabels rows silently. **Trigger:** when writing the group-4 feature scripts, hoist these into a small shared label-constants block (e.g. in `bench_io.py`) and reuse across 36/37/group-4 rather than re-typing literals.

## Full-table citation audit (secretion-classifier-dataset, 2026-06-15)

Two-pass provenance audit over ALL 925 positives (not just the ssign-found 51 from the earlier pass).
Scripts 44-47. Pass 1 = deterministic CrossRef gene/genus check (44); pass 2 = 20 batched agents read
each cited paper and judged SUPPORTED/REFUTED/INACCESSIBLE under an anti-hallucination contract
(`deepverify_input/CONTRACT.md`); 46 merges, 47 applies. User policy = **strict: drop every refuted row.**

Result: **positives_all 925 -> 458.** Removed 467 = 252 pass-1 (wrong-topic 92 / gene-absent 74 / dead-DOI 86)
+ 215 deep-verify refuted (wrong_organism 152 / no_effector_evidence 30 / wrong_protein 27 / wrong_system 6).
Kept 458: verified_paper 330 (with verbatim quote), unverifiable 121 (paywalled, no counter-evidence),
verified_external 3, fallback_consistent 4. New cols on positives_all: `citation_trust`, `citation_quote`.
Backup `positives_all.pre_deepverify.tsv`; removal log `deepverify_removed.tsv`; per-row verdicts
`deepverify_results_full.tsv`; pass-1 verdicts `citation_consistency_full.tsv`.

Headline: only ~36% (333/925) of the original "literature-sourced" rows had a citation that holds up when
the paper is actually read. Many refutes are fabricated-DOI cases (DOI resolves to an unrelated paper:
SptP->GroEL, IcsB->vitamin-D, a T6SS effector->an astronomy paper). Consistent with an LLM-built answer key.

### DOWNSTREAM CASCADE (STALE — must propagate before any dataset/benchmark claim ships)
positives_all shrank 925->458, which invalidates everything built on the old table:
- **Benchmark recall/precision/false-negative figures** (`data/phase2/figures/summary/01-05`,
  `precision/01-03`) and the per-effector tables (`ceiling_per_effector.tsv`, `actual_per_effector.*.tsv`)
  use the old positives. The recall DENOMINATOR shrinks; recall % almost certainly RISES (most dropped rows
  were never in the ssign-found set). **Trigger:** before presenting any recall number, re-run the ceiling +
  actual + figure scripts against the 458-row table. Decide first whether the benchmark should measure recall
  over (a) only citation-verified effectors, or (b) verified+unverifiable. Recommend (a) for the headline,
  (b) as a sensitivity check.
- **40_pair_features.tsv** (was re-run to 925 rows) needs re-run to 458.
- **secretion-classifier-dataset group-4 feature matrix** (gated on CX3 runs) must build on the 458-row table.
- **evidence_tier vs citation_trust**: these are now two separate axes. validated/predicted is the curator's
  claim; citation_trust is what the paper actually supports. 67 rows are evidence_tier=validated but were
  refuted (mostly wrong_organism) and got dropped — so the surviving "validated" set is itself cleaner now.
- **The earlier found-set audit (scripts 41/42/43, positives 930->925) is now subsumed** by this full-table
  pass. Don't double-count its removals; 47 operates on the post-43 table (925) as input.

## T3SS on by default (requested 2026-06-15)
Teo: make T3SS detection ON by default in ssign and test it, then remove the "T3SS off by default"
footnote from the recall figures. The benchmark already validates a T3SS-on path (panel_genbank_t3ss
tag); MacSyFinder detects ~30 T3SS in the panel and the false-positive risk (DeepSecE flagellar
misclassification, 1808 spurious T3SS calls on the dev set) is the historical reason it was excluded
by default. **Trigger:** ssign-side default flip = drop T3SS from `excluded_systems` default in
constants.py; re-run the panel; confirm T3SS precision (the DeepSecE-misclassification concern) is
acceptable before shipping the default change. Figure footnote removed from 52_system_recall.py now.

## 2026-06-15 poster figures (deferred simplify)
- Skipped: full 4-agent simplify pass on the poster edits (52_system_recall.py render() refactor; 03_plot_fragmentation.py fig_04_recovery_poster; build_paired_heatmaps.py poster gate). Edits are small font-scale multipliers + a poster output gate; figures verified by eye. Trigger to revisit: next simplify sweep over benchmark scripts.

## 2026-06-16 T6SS audit (task #74) deferred fixes
Reconciled file: validation_sweeps/benchmark/data/phase2/verification/reconciled_T6SS.tsv (37 proteins R190-R226; 0 fabricated, 0 wrong genomes). Per-row suggested_fix is in that table. Two items need a decision beyond a citation swap:
- R207 BopC / R208 BopE (B. thailandensis E264, T6SS-5 / T6SS_09): all 4 audit passes agree these are misattributed. UniProt accessions point to a TssG baseplate subunit (BTH_II0865) + a 107aa uncharacterized protein (BTH_II0874); cited DOI (chom.2024.04.012) is the TssM-esterase paper and never mentions them; BopC/BopE are canonically T3SS(Bsa) effectors. Decision needed: DROP both from the gold set, or reclassify the B.thai T6SS-5 effector to VgrG5. Trigger: post-audit fix pass.
- Serratia Ssp/Rhs rows R221-R226 (SMDB11_ locus tags): source papers (English 2012, Fritsch 2013) use SMA tags; agent3 could not cross-confirm the SMDB11->SMA mapping and flagged possible drift for Ssp5 (SMDB11_2278, may be RhsI1 immunity region) and Ssp6 (SMDB11_4259 vs 4673). Proteins/refs are real; only the locus_tag provenance is unverified. Trigger: post-audit, spot-check the 6 SMDB11 tags against the Db10 GenBank.
- R214 YPK_3548: effector status rests on a T6SS4-expression/regulation paper (ppat.1013356), not a secretion/translocation assay; needs a better primary ref + a non-deleted accession.

## 2026-06-16 T5SS audit (task #74) — FINAL batch, full 245 audit complete

T5SS batch (19 proteins R227-R245, self-secreting autotransporters/TPS, in_gold_set=no): 4 independent passes (claude + 3 blind agents). **Fully clean: 0 fabricated, 0 wrong genome, 0 wrong/deleted UniProt, 0 wrong_ref.** 14 verified, 5 weak_ref (citation is on-topic but not the ideal secretion primary): R232 fhaB (prodomain-mechanism paper), R233 flu/Ag43 (chaperone-QC paper), R234 hmw1A (cites HMW1B translocator not HMW1A), R235 hpmA (cross-complementation paper), R237 iga (IgG3-specificity paper; UniProt Q51163 is a 496aa fragment of the ~1500aa AT). All 15 supplied accessions + all genomes correct (incl. correctly plasmid-encoded espP/pO157 and yadA/pYV). See reconciled_T5SS.tsv.

### Consolidated deliverable (deferred fixes, all 6 batches)
`validation_sweeps/benchmark/data/phase2/verification/audit_fixes_consolidated.tsv` — 101 rows needing a provenance fix, built from the 6 reconciled_*.tsv by /tmp/build_audit_summary2.py.
- 245 proteins audited, **0 fabricated, 0 wrong genome-accessions** (the 2 T1SS "wrong_organism" are citation-species mismatches, genome fine).
- 144 clean/verified (incl. 12 T2SS secretome-verified).
- 101 need a fix: 39 weak_ref, 24 wrong_ref, 14 wrong/bad DOI, 13 deleted_uniprot, 7 wrong_uniprot, 2 wrong_organism, 2 misattributed (T6SS R207/R208 BopC/BopE = real loci but T3SS effectors, drop-or-reclassify). 4 rows carry a deleted_uniprot on top of a ref defect.
- Verified-DOI replacements already captured in the per-batch reconciled files; weak_ref suggested_fixes are mostly direction-only (marked "unverified") — Crossref-confirm before writing any replacement DOI into the answer key.

**Conclusion: the recall-figure biology is trustworthy (no hallucinated proteins, no wrong genomes). The answer-key's provenance fields (DOIs + a handful of UniProt accessions) need repair, and BopC/BopE need a drop-or-reclassify call.**

## 2026-06-16 Answer-key fix-finding (task #75-79): fresh, blind-agent + reconcile

Per user "do a fresh finding, don't rely on old suggested DOIs". 9 blind agents (2x ClassA-DOI, 2x ClassB1 T1/T2, 1x B2 T3, 2x B3 T4/T5/T6, 2x ClassC uniprot) + my own pass, all finding via PubMed search -> abstract -> Crossref/EuropePMC title gate (never from memory). Reconciled to `validation_sweeps/benchmark/data/phase2/verification/reconciled_fixes.tsv` (101 rows). Final DOI set re-verified in one batch: 0 unresolved.

Counts: 91 auto-apply (74 DOI replacements + 20 uniprot remaps + 4 uniprot->blank + 9 locus_tag fixes), 2 keep_current (R060 lip, R066 zmpA: no better primary exists), 8 FLAG (need Teo's call).

Fresh-pass catches that the OLD suggested fixes got wrong (vindicates the rerun):
- R177 VirE3: old hint DOI 10.1073/pnas.0500396102 is DEAD (404). Correct = Vergunst 2003 Plant Physiol 10.1104/pp.103.029223.
- R113 VopS: old hint accession Q87GE8 is VPA1368 (T3SS inner-rod, 85aa), WRONG. Correct = Q87P32 (reviewed VopS, 387aa).
- R187 VceC: row locus is BAB2_0123 -> Marchesini 2011 (locus-exact); old hint de Jong 2008's VceC is BAB1_1058 (different locus).
- R006 aprA -> Liehl 2006 (P. entomophila, correct organism, not the P. aeruginosa homolog).
- R014/R015 lktA -> Davies 2002 (covers B. trehalosi + M. glucosida, fixes the wrong-organism citation).

NEW defect class uncovered: 9 wrong locus_tags in the answer key (not just citations):
- Ralstonia swaps: R104 PopC RSc0608->RSp0875, R106 AvrA RSp1130->RSc0608, R107 RipJ RSp0871->RSc2132 (RSp0871 was HrpD apparatus).
- R061 chiY YE3650->YE3576; R160 MavQ lpg2813->lpg2552; R161 MitF lpg2818->lpg1976; R168 SdjA lpg2155->lpg2508; R174 CinF CBU0041->CBU0513; R185 BtpB BAB1_0782->BAB1_0756.
- Also locus-uncertain (DOI fixed, locus flagged): R162 PieE lpg1924?->lpg1969, R163 PieF lpg1959?->lpg1972.

4 uniprots genuinely purged by UniProt (proteome removed, no same-locus replacement) -> set blank, keep refseq+locus: R062 engY(GenBank CAL12863), R063 ye3650(CAL13678), R064 cbpE(ABR85197), R213 YPK_0952(NCBI Gene 6091262).

8 FLAG / judgment calls pending Teo:
- R207 BopC, R208 BopE: misclassified T6SS dup; gold set already has correct T3SS/Bsa rows 203/204. RECOMMEND DROP.
- R214 YPK_3548: demonstrated substrate is YPK_3549=YezP (=R215); likely DUPLICATE or wrong locus. DROP-or-relocus.
- R179/R186 TcpB/Btp1: no VirB-translocation assay exists; keep-weak (Salcedo 2008, best-available) or drop.
- R159 MavF: lpg2391=SdbC vs MavF=lpg2351 name/locus collision; decide which protein is intended.
- R162/R163 PieE/PieF: DOI applied; CONFIRM locus (canonical lpg1969/lpg1972).

APPLY not yet run. Targets = data/gold_build/effector_gold_set.tsv (primary_ref/uniprot/locus_tag for T1-4,6) + data/dataset/positives_all.tsv (T5SS + superset). Follow script-43 pattern: backup + identity-matched edits + reversible log + quarantine for drops. Then regenerate recall_figure_proteins.tsv (script 54) + re-run citation_consistency check.

## 2026-06-16 APPLIED (script 55) + figure-table regen (script 54)

Ran scripts/55_apply_audit_corrections.py. Decisions (Teo, this session): drop R207/R208 (BopC/BopE T6SS dups), R214 (YPK_3548 dup of R215), R179/R186 (TcpB, no VirB-translocation assay); R159 MavF->SdbC (lpg2391 IS sdbC, a Dot/Icm substrate; uniprot Q5ZSX5; DOI Huang 2011); keep_current R060/R066.

Applied to data/gold_build/effector_gold_set.tsv (582->577) + data/dataset/positives_all.tsv (T5SS DOIs): 71 primary_ref + 25 uniprot + 9 locus_tag + 1 gene edits = 113 changelog entries, 5 quarantined.
Reversibility: *.pre_audit_fix.tsv backups, effector_gold_set.removed_audit.tsv (quarantine), audit_fix_changelog.tsv.

Script 54 reworked: the recall-figure table backbone is the ssign actual_per_effector tables (pre-audit identity), so it now (a) sources provenance from the corrected gold row via gold_provenance(), (b) applies the MavF->SdbC rename via a changelog-derived gene_alias, (c) excludes quarantined effectors via dropped_keys from removed_audit. Regenerated recall_figure_proteins.tsv: 245 -> 240 proteins (T4SS 42->40, T6SS 37->34). Verified: all 5 drops gone, corrected accessions/loci/DOIs surface (VopS Q87P32, SdbC, CinF Q83E23, VirE3 Vergunst2003, RipJ RSc2132, MitF lpg1976).

### STILL TODO (follow-ups)
- **Plotted figure 52 (52_system_recall.py) NOT yet regenerated.** It reads actual_per_effector + instance classification directly, independent of gold drops. The 5 dropped effectors slightly change instance-level recall; script 52 needs the same dropped_keys exclusion for the plotted figure to match the corrected answer key. (T4SS/T6SS bars affected.)
- positives_all.tsv only had its T5SS rows touched here; the T1-4,6 citation fixes were NOT mirrored into positives_all (it's the secretion-classifier-dataset training table with its own audit, script 43). If training labels should share the corrected DOIs/accessions, run a parallel sync.
- R162 PieE / R163 PieF: DOI fixed; locus left as-is (key lpg1924/lpg1959 vs canonical lpg1969/lpg1972) - confirm if it matters.
- simplify review on scripts 54/55: self-checked (55 follows script-43 pattern + bench_io; 54 edits reuse read_tsv, add one helper). Full 4-agent simplify not run (context).

## 2026-06-16 RE-AUDIT of corrected answer key (task #81) — acceptance test, STRICT bar

Method: re-audit all 240 corrected recall-figure proteins, my manual pass + 4 blind agents per batch, reconcile union. Worklists in data/phase2/verification/reaudit/ (batch_*.tsv from the regenerated recall_figure_proteins.tsv; agent{1-4}_<TYPE>.tsv; claude_<TYPE>.tsv; reconciled_<TYPE>.tsv).

USER DECISION (evidence bar): **STRICT + DROP** rows that lack a direct same-species secretion assay naming the protein (wrong-species citation, or sequence/structure/review/homology-only with no secretion assay). Keep rows with a same-species secretion-mutant/secretome assay (operon-level within the correct species/system counts as kept).

HEADLINE so far (T1SS+T2SS = 75 rows, 5 passes each): provenance 100% clean — 0 fabricated, 0 wrong accession, 0 wrong genome, 0 dead/wrong DOI. The corrected key passes the hallucination axes. The only issues are evidence-strength (the strict-drop targets).

### Running STRICT-DROP list
- T1008 cya (B. bronchiseptica): cites Glaser1988 = B. pertussis (wrong species; no bronchiseptica CyaA secretion paper exists).
- T1012 lktA (B. trehalosi): Davies2002 = lkt-operon sequence/evolution, no secretion assay.
- T1013 lktA (M. glucosida): same Davies2002, no secretion assay.
- T2042 plcB (P. aeruginosa PA0026): cited to Filloux2011 REVIEW; primary Barker2004 shows Sec secretion NOT T2SS -> possible misclassification.
- T2054 zmpA (B. cenocepacia): mic.0.26243-0 = homology-only ("likely GSP-secreted"), no T2SS-secreton-mutant assay.

### Fix-not-drop (T2SS)
- T2023 lipA (X. euvesicatoria): real Xps-T2SS substrate (JB.00322-15), but accession A0ABS8LNA2 points to a different assembly (LN463_12775) than stated XCV0536/NC_007508.1. Needs accession re-map, not a drop.

### Borderline KEPT (same-species secretome/operon evidence, strict-pass)
T1SS indirect-but-kept: T1001 Serralysin, T1003 aprA(entomophila), T1006 apxIIIA, T1015/T1016/T1017 prtA, T1018 prtB. T2SS indirect-but-kept: paeY, pelA, pelD, pelL, paAP, lasA, lapA (same-species Out/Xcp/Hxc secretome or foundational secreton paper).

### REMAINING re-audit batches: T3SS(72), T4SS(40), T5SS(19), T6SS(34) — agents not yet launched.
Then: consolidated strict-drop pass (extend script 55 pattern), regenerate figure table, final clean report.

### Re-audit (task #81) — T3SS batch done 2026-06-16
- **Agent channel BLOCKED**: all 4 blind agents (and a reframed probe) tripped the usage-policy content filter ~8 tool-calls in, on the T3SS pathogen-virulence content (plague/Yersinia/Shigella/Burkholderia). Same expected for T4/T5/T6 (intracellular pathogens). Substituted a rigorous main-session pass: every DOI resolved via crossref+europepmc, every UniProt accession resolved live, BopA row checked via PubMed. Redundancy (4 independent agents) NOT achievable for these batches by that route — user must know method deviated.
- **Provenance 100% clean** (34/34 DOIs real+on-topic T3SS papers; 49/49 UniProt accessions correct incl. held locus fixes RipJ→RSc2132, AvrA→RSc0608, PopC→RSp0875, VopS→Q87P32).
- **T3SS verdict: 54 keep / 3 drop / 15 fixable.** File: reaudit/reconciled_T3SS.tsv.
  - DROP (3): T3005 BopA (cited DOI is the TssM T6SS-DUB paper Szczesna 2024, not BopA/not T3SS; instance tag "T6SS-5" on a T3SS row — incoherent), T3013 (EXACT dup of T3011 EspA/ROD_30191), T3025 (EXACT dup of T3023 EspZ/ROD_30281).
  - FIX-not-drop (15): A/E nle/esp ortholog rows in EPEC or C.rodentium cited to a sister-species paper (mostly O157 Tobe 2006 pnas.0604891103, or the EPEC Map/NleA paper). Real verified effectors; same-species primaries exist. Rows: T3011,T3016,T3018,T3031,T3033,T3035,T3036,T3041,T3042,T3044,T3045,T3047,T3048,T3050,T3051.
- **DECISIONS PENDING (before T4/T5/T6):** (a) the 15 fixable — swap to same-species refs / accept as panel-convention / drop? (b) agent-block: accept main-session pass for remaining pathogen batches, or pause?
- Running drop list now: T1008,T1012,T1013,T2042,T2054,T3005,T3013,T3025. Fix-not-drop: T2023 lipA accession + the 15 T3SS rows above.

### Re-audit (task #81) — ALL 6 batches done 2026-06-16 (main-session pass)
HEADLINE: **provenance 100% clean across all 240 rows** — every DOI resolves to a real on-topic paper (crossref/europepmc), every UniProt accession correct (live), all prior locus/accession fixes held (RipJ→RSc2132, AvrA→RSc0608, PopC→RSp0875, VopS→Q87P32, CinF→Q83E23, MavF→SdbC, VceC→Marchesini, YezP→YPK_3549). 0 fabricated, 0 wrong-genome. All remaining issues are evidence-STRENGTH (strict bar) + a few duplicates + citation refinements.

Per-batch (reconciled_<TYPE>.tsv): provenance clean everywhere.
- T1SS 21: drop 3 (T1008,T1012,T1013). 
- T2SS 54: drop 2 (T2042,T2054); fix accession T2023 lipA.
- T3SS 72: drop 3 (T3005 BopA mis-sourced, T3013+T3025 exact dups); FIX 15 cross-species A/E orthologs (same-species refs).
- T4SS 40: drop 0; FIX 5 (T4016 Lem8 locus lpg0945->lpg1290?, T4025 PieE lpg1924->lpg1969?, T4026 PieF lpg1959->lpg1972?, T4039 VirE3 DOI is a VirE2 paper, T4007 BtpB verify translocation assay).
- T5SS 19: drop 0; FIX 1 (T5011 iga cross-species gonorrhoeae->meningitidis).
- T6SS 34: drop 1 (T6021 Tle4 = dup of T6022/PA1510); FIX T6027 TseH (bioRxiv preprint->published Altindis 2015), normalize T6002 PMID:16432199->10.1073/pnas.0510322103; VERIFY T6020 Tle3 ref + T6033 YPK_0952; DEFINITIONAL FLAG T6001/T6002 Hcp = structural tube protein but secreted (hallmark), generic gene names EFF00001/EFF00006 — user call: keep-as-secreted vs reconsider.

TOTALS: 9 drops (T1008,T1012,T1013,T2042,T2054,T3005,T3013,T3025,T6021); ~24 citation/locus/accession fixes; 3 verifies; 1 definitional flag (Hcp x2).
USER DECISIONS (this session): main-session pass OK for pathogen batches (agents blocked); FIX the 15 A/E orthologs with same-species refs (Citrobacter rows with genuinely no same-species secretion assay fall back to strict drop, flag them).
NEXT: fix-finding pass (verify+source each fix), then consolidated apply (extend script 55 quarantine+changelog pattern) for the 9 drops + fixes, regenerate recall_figure_proteins.tsv (script 54), then plotted figure 52 (task #80).

### Decision 2026-06-16: DROP Hcp rows T6001 + T6002 (user call)
Hcp is a structural T6SS tube protein (secreted as hallmark but not a cargo effector); generic gene names. User: drop both. Drop list now 11: T1008,T1012,T1013,T2042,T2054,T3005,T3013,T3025,T6021,T6001,T6002. Still open: BtpB (T4007, verifying), Citrobacter A/E orthologs (drop-vs-keep).

### BtpB resolved 2026-06-16: KEEP (not drop)
Salcedo 2013 (fcimb.2013.00028) abstract: BtpB "is a novel Brucella effector that is translocated into host cells." Same-species translocation claim present -> clears the bar TcpB/Btp1 failed (those had no translocation assay). T4007 BtpB = keep with current ref.
Citrobacter A/E orthologs: NO new decision needed — prior decisions ("fix same-species refs" + "drop thin") determine it: fix-if-same-species-assay-exists, else strict-drop. Execute as part of fix-finding.
Drop list stays 11. Correction phase = 11 drops + fix-finding (EPEC orthologs ~5 have refs; Citrobacter ~9 search-then-fix-or-drop; T4 loci x3; VirE3; T5011 iga; T6027 TseH; T2023 lipA; normalize T6002 PMID->DOI).

### Fix-finding started 2026-06-16 — checkpoint after 1 result + 1 near-miss
VERIFIED so far:
- T3047 NleE (EPEC) -> Nadler 2010 PLoS Pathog "NleE inhibits NF-kB activation" 10.1371/journal.ppat.1000743 (abstract: NleE INJECTED by EPEC into host cell = same-species translocation). GOOD FIX.
WARNING / do-not-use:
- PMID 15496394 is Mundy 2004 J Med Microbiol (espI epidemiology survey, 10.1099/jmm.0.45684-0), NOT the NleA/EspI translocation paper. The real NleA identification+translocation = Gruenheid 2004 Mol Microbiol (find correct PMID). Caught by metadata check — reaffirms: verify every PMID->paper, never trust the search-rank guess.
PENDING fix-finding (task #82): EPEC NleA(T3035), NleC(T3041), NleD(T3044), NleF(T3050); O157 Map(T3031); Citrobacter x9 (T3011 EspA,T3016 EspH,T3018 EspJ,T3033 Map,T3036 NleA,T3042 NleC,T3045 NleD,T3048 NleE,T3051 NleF) = search same-species, fix-or-strict-drop. Then task #83 (T4 loci/VirE3/iga/TseH/lipA/T6002-normalize) + task #84 (apply 11 drops + fixes, regen 54, then 52).
RECOMMEND: /compact before continuing fix-finding (precision-critical per-paper work; context heavy from full audit). All state durable here + reconciled_*.tsv.

### Fix-finding T3SS A/E orthologs COMPLETE 2026-06-16 (all 15 verified via PubMed get_article_metadata)
Every DOI below confirmed = real paper, abstract shows the claimed same-species secretion/translocation/mutant. None drop; all keep with same-species (or same-system for EspA) ref.

EPEC E2348/69 (5):
- T3035 NleA  -> Thanabalasuriar 2009 Cell Microbiol 10.1111/j.1462-5822.2009.01376.x ("NleA...type III translocated...during EPEC infection")
- T3041 NleC  -> Pearson 2011 Mol Microbiol 10.1111/j.1365-2958.2011.07568.x ("Delivery of NleC by the T3SS of EPEC")
- T3044 NleD  -> Creuzburg 2017 Infect Immun 10.1128/IAI.00620-16 ("NleD...translocated into host enterocytes"; EPEC-infected cells)
- T3047 NleE  -> Nadler 2010 PLoS Pathog 10.1371/journal.ppat.1000743 (NleE injected by EPEC) [prior]
- T3050 NleF  -> Pallett 2014 Infect Immun 10.1128/IAI.02131-14 ("T3SS-dependent translocation of NleF" by EPEC)

O157:H7 Sakai (1):
- T3031 Map   -> Tobe 2006 PNAS 10.1073/pnas.0604891103 (O157 Sakai repertoire; 28 confirmed by translocation assay incl. Map/IpgB family) = same paper as other O157 rows

Citrobacter rodentium (9):
- T3011 EspA  -> Deng 2015 J Bacteriol 10.1128/JB.02401-14 (EspA/B/D translocator T3 secretion, A/E T3SS incl. C.rodentium). NOTE: secretion experiments in EPEC; C.rodentium named as same system. SYSTEM-LEVEL same-species (flag for Teo). EspA = translocon filament, in-scope (cf. EPEC EspA row T3012 kept).
- T3016 EspH  -> Mundy 2004 Infect Immun 10.1128/IAI.72.4.2288-2302.2004 (C.rodentium espH mutant, in vivo)
- T3018 EspJ  -> Dahan 2005 Infect Immun 10.1128/IAI.73.2.679-686.2005 (EspJ translocated TTSS effector; C.rodentium mouse infection dynamics)
- T3033 Map   -> Mundy 2004 10.1128/IAI.72.4.2288-2302.2004 (C.rodentium map mutant; significant colonization defect)
- T3036 NleA  -> Mundy 2004 10.1128/IAI.72.4.2288-2302.2004 (identifies EspI=NleA, T3SS-dependent secretion IN C.rodentium = direct same-species secretion assay)
- T3042 NleC  -> Sham 2011 Infect Immun 10.1128/IAI.05033-11 (delta-nleC C.rodentium mice -> worsened colitis)
- T3045 NleD  -> Kelly 2006 Infect Immun 10.1128/IAI.74.4.2328-2337.2006 (nleD deletion constructed+tested in C.rodentium; null colonization). WEAK same-species genetic (flag).
- T3048 NleE  -> Kelly 2006 10.1128/IAI.74.4.2328-2337.2006 (nleE deletion in C.rodentium; NleE shown translocated by A/E LEE-T3SS). same-species genetic.

FLAGS for Teo (keep, but evidence is system/genetic-level not a strain-specific functional secretion assay): T3011 EspA-Cr (Deng2015 experiments are EPEC), T3045 NleD-Cr + T3048 NleE-Cr (Kelly2006 null deletions). All defensible under "operon-level same-species counts as kept"; surfacing for transparency.

Drop list UNCHANGED at 11: T1008,T1012,T1013,T2042,T2054,T3005,T3013,T3025,T6021,T6001,T6002.
Machine-readable fix table: reaudit/fixes_verified.tsv (built next). NEXT task #83: T4 loci (Lem8/PieE/PieF), VirE3 DOI, T5011 iga, T6027 TseH, T2023 lipA accession, normalize T6002 PMID->DOI, verify T6020 Tle3 + T6033 YPK_0952.

### Fix-finding T4/T5/T6 + misc COMPLETE 2026-06-16 (task #83) — all verified live
Net result: only 3 applicable fixes; several of my own re-audit flags were FALSE ALARMS (caught by verifying the actual abstracts/loci, not the search-rank guess). This is the value of the pass.

APPLY (3):
- T4016 Lem8: locus lpg0945 -> lpg1290 (+ accession -> Q5ZVZ8). Confirmed: eLife 2022 PMID 35175192 ("Lem8 (Lpg1290), 528aa Cys protease, Phldb2") + Huang 2011 (cited ref) + UniProt Q5ZVZ8 lpg1290 528aa YopT-peptidase. Row's lpg0945 = legL1 (wrong). Ref 01531.x kept.
- T6027 TseH: DOI 10.1101/868539 (bioRxiv preprint) -> 10.1128/mBio.00075-15 (Altindis 2015 mBio; "VCA0285 (TseH)...secreted by T6SS" V.cholerae). [it's mBio not PNAS]
- T2023 lipA: accession A0ABS8LNA2 -> Q3BXQ3 (reviewed LIPA_XANE5, XCV0729, 337aa secreted lipase) + locus XCV0536 -> XCV0729. A0ABS8LNA2 = LN463_12775 420aa, NOT lipA (wrong protein). Ref JB.00322-15 (Sole 2015) CONFIRMED correct: abstract names "a lipase...virulence factor" as Xps-T2S substrate.

VERIFIED-NO-CHANGE (my re-audit re-flags were wrong; current refs already correct):
- T4039 VirE3: pp.103.029223 (Vergunst 2003) DOES assay VirE3 translocation by VirB/D4 (CRAFT): "C-terminal 50 aa of VirE2 and VirE3 sufficient to mediate Cre translocation". KEEP. (My reconciled-T4SS note "it's a VirE2 paper" was WRONG.)
- T6020 Tle3: fmicb.2019.01218 (Berni 2019) DOES characterize Tle3 secretion ("secretion mechanism of Tle3...VgrG2b spike...H2-T6SS"). KEEP.
- T6033 YPK_0952: spectrum.04278-23 (Yang 2024) names "YPK_0952...effector of T6SS-3...secreted by T6SS-3...DNase". KEEP.

FLAG-NO-CHANGE (don't apply unverified locus changes that conflict with genome annotation):
- T4025 PieE, T4026 PieF: literature gives loci (mSphere 2024 "PieF (Lpg1972)"; PieE~lpg1969) but UniProt annotates lpg1972 = 125aa "Dot protein" and lpg1969 = LegC3 635aa -- CONFLICT with the genome the benchmark fetches from. Keeping current row loci (PieE=lpg1924, PieF=lpg1959); flag for Teo to reconcile against the specific assembly. Effectors are real; DOIs already fixed prior pass.
- T5011 iga (N.meningitidis MC58): keep Pohlner 1987 (325458a0) founding IgA-protease autotransporter mechanism. Cross-species (gonorrhoeae) but acceptable for a T5aSS self-secretor under the lenient T5SS bar (mechanism papers OK). van Ulsen 2001 alt is only a genome ORF survey (weaker). Flag.

### CONSOLIDATED correction set for task #84 (apply)
DROPS (11): T1008,T1012,T1013,T2042,T2054,T3005,T3013,T3025,T6021,T6001,T6002.
FIXES (18 total): 15 T3SS A/E ortholog DOIs (see prior section) + T4016 Lem8 locus/acc + T6027 TseH DOI + T2023 lipA acc/locus.
Build reaudit/fixes_verified.tsv next, then extend script 55 pattern -> apply -> regen script 54 -> regen figure 52 (task #80).

### APPLIED re-audit corrections (task #84) + figure regen + plotted-figure regen (task #80) — 2026-06-16
Scripts: 56_apply_reaudit_corrections.py (gold), edits to 54 (figure table) + 52 (plotted instance-recall).

CAUGHT A MIS-DROP before it shipped: first run of script 56 used script-55's locate() which pre-filters by
ss_type. The batch file tags T3005 BopA as T3SS, but the incoherent BopA gold row is ss_type=T6SS (locus
BTH_II0876); meanwhile a LEGIT B.pseudomallei T3SS BopA (Q63K42/BPSS1524) also exists. locate() matched and
dropped the LEGIT one. Reverted (git-clean: gold was uncommitted working-tree, removed_audit untracked;
restored from pre_reaudit2 backup), rewrote script 56 matching to GENOME identity (gene+locus/uniprot, ss_type
NOT used) + organism-mismatch guard + grouped multi-field fixes. Re-ran: correct BopA (B.thailandensis T6SS)
dropped, B.pseudomallei kept. Audited all 27 applications: every T3SS DOI landed on the right species row.

Gold: 577 -> 568 (9 quarantined). Field edits: 16 primary_ref + 2 locus_tag + 2 uniprot = 20.
Reversibility: effector_gold_set.pre_reaudit2.tsv (backup), effector_gold_set.removed_reaudit2.tsv (9 drops),
audit_fix_changelog_reaudit2.tsv (29 entries). removed_audit.tsv left at round-1 (5 rows); script 54/52 read BOTH.

Script 54 changes: (a) dropped_keys -> dropped_id keyed by genome identity (gene,uniprot)+(gene,locus), NOT
(ss_type,gene) -- required because cya (B.bronchiseptica drop vs B.pertussis keep) and lktA (drop 2 of 4 hosts)
share a gene; coarse key would over-drop. Reads removed_audit + removed_reaudit2. (b) instance-dedup: collapse
same (ss,gene,locus,organism) keeping 'found' over 'unreach' -> removes the 2 EspA/EspZ-Cr T3SS_28 dup rows
(logged, not silent).
recall_figure_proteins.tsv: 240 -> 230 (8 drops-in-figure + 2 dedup; BopA-thai wasn't in figure so its gold
drop removes 0 figure rows). Per-type: T1SS 18, T2SS 52, T3SS 70, T4SS 40, T5SS 19, T6SS 31.

Script 52 (plotted instance-recall 06_recall_systems.png): added same dropped_id filter (dedup N/A at instance
level). Instance-recall before->after audit drops (excl T5SS): found 36/46 -> 33/43 reachable. Deltas: T1SS
testable 21->18 (cya-Bb + lktA x2 were singleton-effector instances; recall RATE unchanged 100%), T2SS 12->11,
T4SS 10->9 (round-1 TcpB/Btp1). T3SS/T6SS unchanged (dropped effectors shared instances or were instance-less).
Final plotted: ssign found 48/60 reachable systems (incl T5SS 15/17).

TASK #81 ACCEPTANCE TEST RESULT: provenance 100% clean across all 240 (0 fabricated/wrong-genome/dead-DOI);
strict-bar corrections = 11 drops + 18 verified fixes, all applied + propagated. Re-audit COMPLETE.
Open flags for Teo (keep, low-priority): T3011 EspA-Cr (Deng2015 system-level), T3045/T3048 NleD/E-Cr (Kelly2006
genetic), T4025/T4026 PieE/PieF loci (lit vs UniProt annotation conflict), T5011 iga (cross-species mechanism).

## 2026-06-19 — Permutation refinement of the enrichment test (task #70, commit a8bca53)

Validated ssign's binomial enrichment test against a clustering-aware circular-permutation
spatial null on the 67-genome fleet. Scripts: validation_sweeps/benchmark/analysis/fleet_67/
03a_regen_neighborhoods.py (macsyfinder regen of per-genome +/-3 neighborhood masks, cached in
neighborhoods/) + 03_permutation_refined.py. Figures 06/07/08.

FINDING: the production binomial is anti-conservative. Of 87 systems it calls significant
(q<0.05), the masked spatial null confirms only 19 (22%). Cause: the binomial assumes each
neighborhood gene is independently secreted at the genome background rate, but secreted genes
cluster (operons/islands), so dense neighborhoods are less surprising than the binomial thinks.
Masking other systems' neighborhoods out of the null RAISED the spatial count 11->19 (vs script
02's unmasked floor), so masking is the fair test, not a power loss. n=1000-sampled permutation
recovers 95% of exact calls (24 vs 19, slightly liberal) — Monte-Carlo placement is fine, it's
the independence assumption that's the problem. Per type: binomial over-call worst for T5SS
(19->2; tiny autotransporter neighborhoods) and T6SS; T3SS-DLP (9/15) and T1SS are robust under
both. DSE over-called less than DLP.

DEFERRED DECISION (affects #69 "enrichment on by default"): the binomial over-states confidence.
Options when whole-genome predictions exist: (a) keep binomial + caveat, (b) switch enrichment
significance to the spatial-permutation null, (c) report both. Permutation needs whole-genome
ordered predictions (not available in default neighborhood-mode) and has a ~1e-3 resolution
floor when Monte-Carlo'd. Revisit when deciding #69. Trigger: before flipping enrichment to
default-on.

## 2026-06-19 — Circular-shift enrichment with fold values + correct T5SS (task #70 cont.)

Teo: the sig/not-sig counts don't show HOW enriched systems are, and asked which T5SS proteins
were measured (hitchhikers vs autotransporters). Found the Xanthobacter-era reference
(/home/teo/Desktop/Billerbeck - SS Identification/: make_enrichment_figures.py figure4/figure5,
t5ss_aware_analysis.py) and reproduced it on the fleet: 04_circular_shift_enrichment.py.

METHOD: circular-shift null (genome-structure-preserving) — rotate each genome's predictor-
positive positions by a random offset, count how many land within +/-3 of an SS component, sum
across genomes, 10k shifts. fold = observed / null-mean. Exact per-genome all-rotations count via
FFT cross-correlation (unit-tested vs brute force, c[0]==observed). Figures 09 (genome-wide null
histogram), 10 (per-type fold), 11 (autotransporter self-detection).

RESULTS (genome-structure-preserving, so lower than binomial but honest):
- genome-wide: DLP 400 vs null 121 = 3.3x, DSE 362 vs null 95 = 3.8x, both p<1e-4.
- per-type DLP all *** except T2SS(1.4x n.s.): T1SS 4.3x, T3SS 3.8x, T5bSS 3.4x, T6SS 3.6x, T4SS 3.0x*.
- per-type DSE: T6SS 6.3x*** (strongest; effectors are DSE-typed), T1SS 4.4x***; T4SS 0.2x n.s.
  (DSE misses translocated T4SS substrates), T5bSS 1.6x / T2SS 1.4x n.s. T3SS-DSE excluded by design.
- enrichment is REAL and strong everywhere except T2SS; binomial inflated the MAGNITUDE (6-12x) not
  the existence.

T5SS FIX (the binomial/03 lumped all T5 -> meaningless for autotransporters):
- T5aSS/T5cSS are autotransporters (protein = machinery AND substrate). Window is mostly empty
  (binomial: only 36/206 T5a, 9/53 T5c systems had any neighbour positive). Correct test =
  SELF-DETECTION: is the component itself OM-or-extracellular (DLP) / secreted-type (DSE)?
  -> T5aSS 35% DLP / 11% DSE ; T5cSS 5% / 5%. DeepLocPro under-calls autotransporters at >=0.8,
  esp. trimeric T5cSS adhesins (YadA-like). Real DLP sensitivity gap.
- T5bSS stays a window type (TpsA hitchhikers near TpsB): DLP 3.4x*** real enrichment.

DEFERRED CLEANUP (Teo: separate task once stats analysis is finalized per-run):
- BUG in fleet results_raw: dlp_max_probability always == dlp_extracellular_prob even when
  dlp_max_localization == Outer Membrane (e.g. NMB1969 max_loc=OM, max_prob=0.294=extra, OM=0.55).
  Cosmetic: is_dlp_positive uses dlp_extracellular_prob directly, so calls are unaffected. Fix in
  the DLP-output writer when cleaning up the enrichment/stats code.
- The per-run enrichment stat to SHIP (#69) should be this circular-shift fold + null, not the
  binomial, with T5a/T5c on self-detection and T5b/others on windows.

## 2026-06-19 — annotation-accuracy sheet (#71) + fleet PAO1 duplicate

ANNOTATION SHEET (recall-gated): known experimental annotation vs ssign's predicted
annotation, manual-eyeball TSV (no LLM/keyword scoring, per Teo).
- script: validation_sweeps/benchmark/analysis/annotation_accuracy_sheet.py
- out:    validation_sweeps/benchmark/analysis/annotation_accuracy_sheet.tsv
- universe = 51 emitted_secreted corpus effectors (only proteins ssign annotates).
- ground truth from positives_all: known_family 11/51, known_quote 42/51 (quote is the
  main signal; family sparse). Genome via ssign_locus -> results_raw scan.
- ssign yield: interpro 39/51 (76%), eggnog 43/51 (84%), pLM-BLAST/ECOD 44/51 (86%),
  Pfam 38/51 (75%); >=1 tool 44/51 (86%); 7 got nothing from any tool (genuine misses,
  not join failures; several are uniprot-less corpus rows with thin ground truth too).
- spot-check: serralysin effectors -> interpro "Metallopeptidase catalytic domain" +
  eggnog "peptidase M10A matrixin"; RsaA -> "RTX Ca-binding / COG2931 RTX toxins". Tools
  land on the right family. Teo to eyeball the full sheet for correctness.

FLEET PAO1 DUPLICATE (deferred check for the enrichment fleet, NOT this sheet):
- /tmp/ssign_fleet_67 contains PAO1 TWICE: AE004091 (INSDC) and NC_002516.2 (RefSeq),
  same genome (~5571 proteins each), identical annotations. So the "67-genome" pooled
  circular-shift enrichment slightly double-counts PAO1. Trigger: when finalizing the
  per-run stat (#69) / re-running pooled enrichment, dedup to 66 distinct genomes.

## 2026-06-19 — circular-shift enrichment shipped per-run (#69)

SHIPPED openspec `enrichment-circular-shift-per-run`: the opt-in `--enrichment-stats`
test is now a per-SS-type circular-shift permutation (fold + permutation p + BH q) +
an auto per-type null-distribution figure. Binomial retired as the significance source.
- enrichment_testing.py rewritten (rotation_counts FFT, window/self masks, per-type
  run_enrichment, new OUT_FIELDS: ss_type/tool/mode/observed/n_mask/null_mean/fold/
  p_perm/qvalue/significant/n_rotations); is_dlp_self_positive added; binom_pvalue +
  broad_type retained only for the validation scripts.
- runner.py: __post_init__ forces dlp/dse_whole_genome when enrichment_stats on (+note);
  _step_sample_null_proteins removed from DAG AND deleted; _step_enrichment rewired
  (no --null-ids, adds --nulls-out npz); new _step_enrichment_figure; pool_enrichment_stats
  rewritten to Monte-Carlo pool the per-genome rotation nulls (ENRICH_POOL_REPS=10000).
- new scripts/run_enrichment_figure.py; constants ENRICH_* added.
- tests: test_enrichment_testing.py + test_runner pooling + integration rewritten; all
  1373 unit tests pass; integration passes.

DEFERRED CLEANUP (now-unwired null-background machinery) — DONE 2026-06-22:
- Retired `sample_null_proteins.py` + its test, the `n_null_proteins`/`null_seed` config
  fields + `--n-null-proteins`/`--null-seed` CLI args, the `include_null_concat` param +
  branch in `_resolve_step_input_fasta`, and the dead `dlp_dse_input` pooling block in
  multi_runner `_pool_segment_b_inputs`. The rotation null is the exact background, so none
  of it was reachable. 1362 unit tests pass; CLAUDE.md key-params line updated.
- VERIFY on a real multi-genome run with --enrichment-stats: that forcing whole-genome
  DLP/DSE composes correctly with multi_runner's neighborhood-pooling/segment logic.
  (IN PROGRESS — the 4-genome CX3 enrichment job is the validation run.)
- PAO1 duplicate in /tmp/ssign_fleet_67 (AE004091 == NC_002516.2): dedup to 66 genomes
  if the fleet is re-run for pooled-enrichment validation.

## 2026-06-22 — annotation sheet: UniProt known-annotation column (#71)

annotation_accuracy_sheet.py now fetches UniProt protein_name/families/function
(cached in annotation_uniprot_cache.json, offline after first run) as the richer
"known" column. Coverage: uniprot_name 30/51, function 14/51 (many effectors are
unreviewed TrEMBL w/o a function comment), + existing corpus quote 42/51 -> 50/51
rows now have >=1 known descriptor (was family 11/51 only). 17/51 effectors have
NO uniprot accession in the panel (recorded "-") so can't be looked up; their known
side stays gene+quote. Sheet already surfaces real misses (e.g. Q7N8R3/Q84F70 known
serralysins: ssign eggnog mis-calls "Serine 3-dehydrogenase"; interpro correct).

## 2026-06-22 — pooled cross-genome enrichment figure (#69 follow-up)

Multi-genome runs now also get a COMBINED figure (the pooled view Teo originally
liked from the Xanthobacter fleet). pool_enrichment_stats gained nulls_output=,
dumping the pooled Monte-Carlo null arrays; Home.py renders
figures/pooled_enrichment_null_distributions.png from the pooled TSV+npz via
run_enrichment_figure.py. Per-genome figures unchanged. Tests: unit (npz keys/shape)
+ integration (pooled figure renders). Part of openspec enrichment-circular-shift-per-run.

## 2026-06-22 — CX3 validation of enrichment runs (handoff to Teo)

Branch enrichment-circular-shift-per-run pushed (4 commits). CX3 test is Teo-driven:
I can't SSH (Imperial login refuses publickey; needs password + MS 2FA). Plan: git
pull the branch on ~/blastp_t5a/ssign, qsub run_batched_multi.pbs with
SSIGN_EXTRA_ARGS="--enrichment-stats" for (a) 1 genome -> per-genome figure, (b) 4
genomes -> per-genome + pooled figure + pooled_enrichment_stats.tsv. Verify:
<run>/<sid>/figures/<sid>/<sid>_enrichment_null_distributions.png (single),
<run>/pooled_enrichment_null_distributions.png + pooled_enrichment_stats.tsv (multi).
KEY THING TO WATCH: whether forcing whole-genome DLP/DSE composes correctly with
MultiGenomeRunner's segment-B prediction pooling (the one unverified interaction).

## 2026-07-15 — annotation-subsystem-cleanup: golden-fixture regen owed

PLM-Effector removed entirely (openspec annotation-subsystem-cleanup). The e2e
golden `tests/fixtures/golden/t5ass_minimal/t5ass_minimal_predictions.tsv` had its
three `plm_effector_*` columns surgically dropped (deterministic column deletion).
BUT that fixture is ALSO stale for a pre-existing, unrelated reason: it lacks the
`cytoplasmic_membrane_prob` DeepLocPro column that current production emits (golden
jumps cytoplasmic_prob -> dse_ss_type; production has cyt_mem between them). So the
opt-in golden test (`-m integration`, needs local licensed DeepLocPro) will still
byte-mismatch on that file until a FULL DeepLocPro-driven regen is run per
`tests/fixtures/golden/REGENERATE.md`. Trigger: next time on a host with SignalP/DLP
installed (SSIGN_DEEPLOCPRO_PATH set), run the regen and commit the refreshed refs.
Also confirm the golden `_EXPECTED_FIGURES` list (still `01`-`07` old names) matches
figures-v2 output while there. Can't do here (no licensed DTU tools on this laptop).

## 2026-07-15 — scripts/analyse_k12_runs.py still parses PLM-E log lines (deliberate)

`scripts/analyse_k12_runs.py` (+ `tests/unit/test_analyse_k12_runs.py`) still contain
PLM-E log-parsing (34 refs). LEFT IN ON PURPOSE: it's a retrospective benchmark
log-analyzer, and historical K-12 validation logs contain PLM-E lines. Ripping it out
would lose the ability to re-analyze past runs. Decision for Teo: keep the historical
parser, or strip it too for a fully-clean tree? Not a pipeline consumer, so out of the
annotation-subsystem-cleanup scope.

## 2026-07-15 — PLM-E residual: scripts/fetch_weights.sh — RESOLVED (deleted)

RESOLUTION (verified + Teo's call): deleted the script. Checked what each of its
6 downloads served: 5 were PLM-E-only (full ProtT5 `prot_t5_xl_uniref50`, ProtBert,
PLM-E trained_models, AND the ESM1b/ESM2 HF downloads — those went to a
`plm_effector/transformers_pretrained/` path that NO kept tool reads: DeepSecE/
DeepLocPro load ESM via `esm.pretrained` from the torch-hub cache, a different
location/format). Only the DeepSecE-checkpoint fetch served a kept tool, and that
overlaps with runtime auto-download. So the script was ~vestigial → deleted; refs in
dependency_manifest.py / scripts/README.md / install_test_runbook.md updated. The
image never used it (bakes weights directly in %post). ORIGINAL NOTE BELOW (obsolete):

annotation-subsystem-cleanup removed PLM-E from src/tests/core-docs, and this
session also cleaned scripts/fetch_databases.sh (removed the ~19 GB dead
fetch_plm_effector_weights + URL + preflight + tier comments), troubleshooting.md,
scripts/README.md, the run_batched_multi.pbs dead export, and an audit comment.
STILL TODO: scripts/fetch_weights.sh. It's entangled — its default-path PLM
fetches (ProtT5/ESM1b/ESM2 via fetch_hf_model) live UNDER a `$TARGET/plm_effector/
transformers_pretrained/` directory, alongside PLM-E-only bits
(fetch_plm_effector_trained_models = the 1.7 GB sourcecode.zip, fetch_protbert,
PLM_EFFECTOR_SOURCECODE_URL, the SSIGN_PLM_EFFECTOR_WEIGHTS log). Decision for Teo:
is fetch_weights.sh meant to pre-stage the DEFAULT-path weights (DeepSecE ESM1b /
DeepLocPro ESM2 / pLM-BLAST ProtT5) for native offline use, or was it PLM-E-only?
If default-path: drop only trained_models + protbert + the plm_effector URL/log,
and rename the `plm_effector/` layout to something tool-neutral (check nothing reads
that path). If PLM-E-only: the whole script goes (the container bakes those weights
itself; native users' tools auto-download). NOTE: scripts/README.md was already
edited to describe fetch_weights as fetching "DeepSecE, ProtT5, ESM" (the
default-path reading), so resolve the script to match or re-edit the README.

## 2026-07-15: v7 image built + validated; branches merged to main (DONE)

**v7 container** (`~/ssign_v7.sif`, 19 GB) built with HH-suite baked (#35) +
build-weight cache (#36). Offline smoke passed: base-tier `doctor` all green;
full-tier binaries 7/8 (bakta/eggnog/hhsearch/hhblits/blastp 2.17/DeepLocPro
resolve from the image; InterProScan the expected 1 FAIL, host-mounted by
design). Model weights 1/1 (DeepSecE checkpoint only, ZERO PLM-E). DBs 0/8 FAIL
= expected (mount at runtime). Only host-provided: SignalP (DTU) + IPS + DBs.
NOT yet run: the CX3 parity proof at extended/full with only DB+DTU mounts
(Teo-driven; expect 23/23). reproducible-install openspec still 17/22 (its
5.10/rebuild + CX3-parity tasks); update after the CX3 run.

**Branch merge (task #30 DONE).** Fast-forwarded `main` to the dev branch
(`enrichment-circular-shift-per-run`) at 6e5c380, pushed main + dev to origin,
deleted `pre-claude-purge` (local-only backup, 0 unique commits) + the remote
`fix/ci-lint-mypy-test` (superseded: its one still-valuable piece, the complete
pyhmmer mypy override with `follow_imports_for_stubs`, was ported into dev).
Recovery SHAs if ever needed: old main da96c1f, pre-claude-purge 7469c90, fix
374a7da. Before merging, made dev CI-green: fixed 10 ruff + 22 mypy errors
(commit 6e5c380), all behavior-preserving, 1398 unit tests pass. NOTE: main ==
dev now; the feature-branch name is stale, so future work can just go on main
(dev branch kept per Teo, but it's a cleanup candidate).

## 2026-07-15: PLM-E doc/metadata sweep + base-tier size drift (DONE)

Finished the PLM-E removal outside src/ (annotation-subsystem-cleanup only did
code+core-docs). Swept it as a live tool from README, pipeline_overview,
troubleshooting, configure, licensing, env_vars, run_on_hpc, install_test_runbook,
models/README, plus the release metadata CITATION.cff + codemeta.json (Zenodo/GitHub
read these) and dead pyproject.toml config (keyword, coverage-omit block, mypy
exclude, per-file lint ignore, all pointed at deleted files). Added a CHANGELOG
"Removed" entry for PLM-E + fetch_weights.sh. Verified: `git grep` shows PLM-E only
in openspec/, NOTES, CHANGELOG/design_decisions (removal notes), the historical
analyse_k12_runs parser+test, test_scaled_timeout (deliberate "unknown tool"
example), and .gbff sequences (amino-acid false positives).

Same commit fixed the base-tier size drift (docs disagreed: 2/4/17/22 GB). Root
cause: some docs counted model weights in "tier size", some counted only DBs, and
none said which. **Standardized decomposition (measured from the v7 build rootfs +
ssign.def bake manifest):**
- base **DBs ~4 GB** (taxdump 1.5 + Bakta light ~2.5); extended ~100 GB; full ~500 GB.
- **weights ~14 GB** auto-downloaded/baked, shared across tiers (ESM2 2.5 +
  DeepSecE ckpt 2.5 + DeepSecE's ESM-1b 7.3 + DeepLocPro ~1.4; ProtT5 ~2.4 is
  extended-only via pLM-BLAST). Old "~18 GB weights" was inflated by PLM-E's
  ProtBert 1.6 + PLM-E weights 1.7 (now gone).
- HH-suite is **full-only** (design_decisions wrongly split Pfam+PDB70 to extended);
  Bakta light belongs to **base** not extended (install.md had it under extended);
  Bakta full ~84 GB not ~30 (install.md prose was stale).

FOLLOW-UPS (small, non-blocking):
- **openspec plme-prediction spec** (`openspec/specs/plme-prediction/spec.md`) still
  describes PLM-E as live. It's retired by the spec-sync when annotation-subsystem-
  cleanup is archived (its proposal.md says "plme-prediction: removed"). Trigger:
  `/opsx:archive annotation-subsystem-cleanup`. Also runtime-effort-model spec still
  lists PLM-Effector in the predictor set (stale vs the code, which dropped it).
- **pool_utils seq_id handling** (`ssign_lib/pool_utils.py` custom-id-column split,
  tested by test_pool_utils.test_split_custom_id_column): the only documented consumer
  was PLM-E (emitted seq_id not locus_tag). Check whether any current tool still emits
  seq_id; if none, the custom-id path + its test are dead code. Trigger: repo TLC audit (#32).
