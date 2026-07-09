# Fresh-install test runbook

Pre-v1.0.0 release QA. The goal is to **install ssign strictly from the published docs, as a brand-new user would**, and record every place the docs fail. Two parts:

- **Part A** — manual fresh-install test, one clean-room install per tier, following the real docs verbatim.
- **Part B** — Docker build + in-container install test, following `containers/README.md` verbatim.

Run Part A first. Part B is gated on Part A **and** the full-tier CX3 smoke test both passing.

## The one rule

**Follow the docs exactly as written. Do not "know better".** If a documented command fails or you have to deviate to make it work, **stop and log it** (table below), then apply the minimal deviation to continue. Every deviation is a doc bug to fix before release, that is the entire point of this exercise. Do **not** consult this runbook for install commands; it deliberately contains none. It only tells you *which* doc to follow, how to *verify*, and how to *wipe* between tiers.

## Doc entry points a new user actually hits

| Read in this order | What it covers |
|---|---|
| `README.md` § Quickstart | first install command + GUI launch |
| `README.md` § Install tiers / § Verify | tier table + `ssign doctor` |
| `docs/tutorials/first_run.md` | first end-to-end run (functional test) |
| `docs/how-to/install.md` | per-tool install for extended/full (Bakta, EggNOG, InterProScan, HH-suite, BLAST+, SignalP, DeepLocPro, pLM-BLAST, PLM-Effector) |
| `containers/README.md` | Docker / Singularity (Part B) |

## Prerequisite: are the docs actually reachable as a new user?

The README Quickstart installs from GitHub (`pip install git+https://github.com/billerbeck-lab/ssign.git@<tag>`). For a faithful test the **repo must be public on GitHub at that tag**. If it is not published yet, either (a) push it first, or (b) run the test against your local checkout *standing in for* GitHub — but then note that "clone/install from GitHub" itself is untested. Log which mode you used.

## Findings log (fill this in as you go)

| Tier | Doc + section | Command as written | What happened | Deviation needed | Fix |
|---|---|---|---|---|---|
| | | | | | |

Known issues to expect (already suspected, confirm + log them):
- README Quickstart says `pip install git+https://…@v0.9.0-prerefactor`, but the tier table's install column says `pip install ssign` — the PyPI name 404s until v1.0.0. **Inconsistent install command in the README.**
- README tier table lists **base ≈ 22 GB**; `install.md` lists **base ≈ 2 GB**. **Size mismatch between the two docs.**
- `ssign doctor` on an HPC login node flags `hhsearch`/`hhblits` missing unless the `hhsuite` conda env is on PATH (the PBS wrappers add it; a hand-run doctor does not).

## Two wipe levels

| Level | Removes | Keeps | When |
|---|---|---|---|
| **Light** (default, between tiers) | the ssign venv + the clone dir | fetched DBs, model weights, DTU/conda tool envs, `~/.ssign` cache | between each tier — re-tests `pip install …[tier]` + dep resolution fresh |
| **Deep** (full clean-room) | the above **plus** `~/.ssign/{models,db_root}`, `~/.macsyfinder`, and optionally the `signalp6`/`deeplocpro`/`hhsuite`/`bakta-deps` conda envs | fetched reference DBs on `$EPHEMERAL` (630 GB, validated; re-pulling is wasteful) | once, to re-test `fetch_weights.sh` + the DTU install sections end-to-end |

Between-tier wipes default to **light**. Re-running the doc's `fetch_databases.sh --tier N` is cheap when the DBs are present (skip-guards short-circuit) and doubles as a test of the skip logic.

---

# Part A — per-tier fresh install, from the docs

Sequence: **tier 1 (base) HPC → tier 1 (base) laptop → confirm both, wipe → tier 2 (extended) HPC, confirm, wipe → tier 3 (full) HPC, confirm (keep).**

For **each** tier below the procedure is identical; only the doc sections and the verify tier change:

1. **Clean environment.** New shell, fresh clone dir, fresh venv, exactly as the doc's Quickstart tells you to create it.
2. **Install** by copy-pasting the doc commands for that tier — README Quickstart/tiers for base; add the relevant `install.md` sections for the extra tools at extended/full. No commands from this file.
3. **Verify**: `ssign doctor --tier <base|extended|full>` — must exit green for the tier you installed. (Base has no required databases; its only tier-specific entries are the *optional* local SignalP/DeepLocPro, so base is green straight after the pip install.)
4. **Functional test**:
   - base → follow `docs/tutorials/first_run.md` to completion (its no-DTU-tools variant uses the DTU webserver, so **no GPU or databases needed**; ~30 min, internet required). `tests/run_tutorial.sh` automates exactly this variant if you'd rather run it unattended, but for the *doc* test, follow the tutorial prose.
   - extended / full → a real run with the tools enabled. On CX3 that needs a GPU node; submit via `scripts/cx3/submit_batched_overnight.sh --tier <extended|full>` on 2 genomes and confirm N/N steps succeeded (retrieve per the cx3-retrieve flow).
5. **Record** doctor result + functional pass/fail + wall-clock in the findings log.
6. **Wipe** (light) before the next tier.

## Per-tier specifics

- **A1 — tier 1 (base), HPC (CX3).** `module load Java/11.0.25` first (the tutorial's MacSyFinder step needs it). Everything else from README Quickstart + `first_run.md`.
- **A2 — tier 1 (base), laptop (ThinkPad).** Same, cross-platform check on AMD / no-NVIDIA-GPU. The tutorial's webserver fallback means no GPU is needed. Also confirm the bare `ssign` command launches the Streamlit GUI (the Easy-Mode entry point the README points at).
- **A3 — confirm both, light-wipe both.**
- **A4 — tier 2 (extended), HPC.** Follow `install.md` for each extended tool + `fetch_databases.sh --tier extended`. The DTU tools (SignalP, DeepLocPro) and EggNOG/pLM-BLAST are separate installs per their `install.md` sections; a light wipe kept the existing conda envs, so a strict test either reuses them (faster) or does a **deep** wipe of just those envs to re-test their install docs. Functional run via the PBS submitter. Confirm, then wipe.
- **A5 — tier 3 (full), HPC (final, keep).** Add `install.md` § HH-suite + § BLAST+ and `fetch_databases.sh --tier full`. `ssign doctor --tier full` should report 8/8 databases. **The functional test for this tier is the full-tier CX3 smoke test already in flight** (2 genomes, HEAD with the EggNOG/Bakta/HH-suite fixes) — cross-reference its result rather than launching a duplicate. Keep this install; the real full-tier panel runs from it.

## Part A pass criteria

- Each tier: `ssign doctor --tier <tier>` green **and** a functional run completes with expected outputs, **using only the documented commands**.
- Every deviation you had to make is logged and turned into a fix to `README.md` / `docs/how-to/install.md` / `docs/tutorials/first_run.md`.

---

# Part B — Docker build + install test (runbook only)

Gated on Part A passing **and** the full-tier CX3 smoke test passing. Follow `containers/README.md` verbatim (same one-rule: log any deviation). The `containers/Dockerfile` is **not** modified except to lock the CUDA base SHA, which is a documented build-time step.

1. **Lock the CUDA base SHA** — the Dockerfile ships a placeholder digest (line 35, all zeros). Resolve + lock it per the FRAGILE comment at the top of the Dockerfile (`docker pull` the base, `docker inspect --format='{{index .RepoDigests 0}}'`, paste the digest into `CUDA_BASE`). This is the only Dockerfile edit.
2. **Build** per `containers/README.md` § Build.
3. **In-container smoke test (no GPU, no DBs needed):** `--help` resolves the entrypoint; `ssign doctor --tier base` runs inside the image (it will report DBs/weights missing — expected, they bind-mount at run time); `pip show ssign` (reports `0.9.0` until the version bump). Log anything that isn't as `containers/README.md` describes.
4. **Full pipeline run in-container** — needs `--gpus all` + bind-mounted DBs/weights, per `containers/README.md` § Run. **Deferred to a GPU host**; the ThinkPad can build + smoke-test but not run the GPU stages.
5. **Singularity (HPC)** — optional, per `containers/README.md` § Run (Singularity).

## Part B pass criteria

- Image builds reproducibly with the SHA-locked base.
- `--help` + `doctor` run in-container.
- (On a GPU host) a full pipeline run against a bind-mounted tutorial genome completes.
- Every deviation from `containers/README.md` is logged + fixed.
