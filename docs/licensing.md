# Licensing — ssign Dependencies & Redistribution

This page documents the redistribution status of every external tool, model,
and database ssign uses. The audits below answer the question
**"Can the ssign Docker image bundle this, or must the user download it
themselves?"** for each component.

ssign itself is **GPL-3.0-or-later** (see `LICENSE`). Every dependency listed
here is compatible with GPL-3.0 redistribution at the level of the ssign
package; the per-component questions below are about whether **ssign's
Docker image** can include the binary/weight/database file inside it, vs.
require the user to fetch it after pulling.

| Component | License | In ssign Docker image? | If not bundled, fetch via |
|---|---|---|---|
| **DeepSecE checkpoint** | (vendored, MIT) | ✅ Bundled (mirrored to Zenodo) | — |
| **PLM-Effector weights** | MIT | ✅ Bundled (mirrored to Zenodo) | — |
| **ProtT5 weights** (Rostlab/prot_t5_xl_uniref50) | AFL-3.0 | ✅ Bundled (mirrored to Zenodo) | — |
| **SignalP 6.0** (binary + weights) | DTU academic | ⏳ Pending DTU response | DTU portal (manual) |
| **DeepLocPro** (binary + weights) | DTU academic | ⏳ Pending DTU response | DTU portal (manual) |
| **EggNOG database** | unspecified | ⏳ Pending EMBL response | `download_eggnog_data.py` |
| **InterProScan tarball + member DBs** | Apache 2.0 (core) + mixed (members) | ❌ Not bundled | `scripts/fetch_databases.sh` |
| **BLAST NR / Swiss-Prot** | NCBI public | ❌ Not bundled (size) | `scripts/fetch_databases.sh` |
| **Bakta DB** | CC-BY 4.0 | ❌ Not bundled (size) | `bakta_db download` |
| **HH-suite Pfam / PDB70 / UniRef30** | mixed (Söding lab + Tübingen) | ❌ Not bundled (size) | `scripts/fetch_databases.sh` |
| **TXSScan models** | CECILL | Bundled with MacSyFinder install | — |

---

## ProtT5 weights — ✅ MIRROR + BUNDLE OK

Used by PLM-Effector and pLM-BLAST. Released under the **Academic Free License
3.0** (AFL-3.0), an OSI-approved permissive license.

> "The ProtTrans pretrained models are released under the terms of the
> Academic Free License v3.0 License."
> — [ProtTrans README](https://github.com/agemagician/ProtTrans)

AFL-3.0 explicitly grants public redistribution. No non-commercial or
field-of-use restriction. Compatible with ssign's GPL-3.0. Ship with the
AFL-3.0 license text and a citation pointer to Elnaggar et al. 2022 (TPAMI).
ColabFold and bio_embeddings already redistribute these weights in their
Docker images — established precedent.

**Action:** mirror to ssign's Zenodo deposit + bake into the v1.0.0 Docker
image.

## InterProScan — ❌ BUNDLE NOT OK, user fetches

Core engine is **Apache 2.0** (redistribution-friendly). The bundled member
databases have mixed licenses, and the show-stoppers are:

- **ProSite (PROSITE Profiles + Patterns)** — SIB custom license: non-commercial,
  no derivatives. Forbids redistribution outside the InterPro bundle.
- **SMART** — EMBL academic-only; "redistribution path is via the InterPro
  consortium" per the SMART FAQ.
- **SUPERFAMILY** — no formal license (flagged ambiguous by reusabledata.org).

**Industry signal:** EBI's own
[`interpro/interproscan` Docker image](https://hub.docker.com/r/interpro/interproscan)
does *not* bundle the data — it expects the user to mount it. Bakta avoids IPS
entirely. No mainstream academic Docker image redistributes the IPS bundle.

**Action:** ssign Docker image has IPS *binary* (Apache 2.0) but no member-DB
data files. `scripts/fetch_databases.sh --tier extended` pulls the tarball
from EBI FTP at first install. User experience is one extra command, run
once.

Sources: [InterPro license](https://interpro-documentation.readthedocs.io/en/latest/license.html) ·
[HowToDownload](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html) ·
[`interpro/interproscan` Docker](https://hub.docker.com/r/interpro/interproscan) ·
[SMART FAQ](https://smart.embl.de/help/FAQ.shtml) ·
[SUPERFAMILY reuse](https://reusabledata.org/supfam.html)

## EggNOG — code is AGPL-3.0, database is ⚠️ AMBIGUOUS

These are two separate things and only one is the licensing question we care
about for the Docker image:

|     | What it is | License | Redistribution OK? |
| --- | --- | --- | --- |
| **eggnog-mapper code** | ~5,000-line Python tool that runs DIAMOND + parses output | AGPL-3.0 ([`LICENSE.txt`](https://github.com/eggnogdb/eggnog-mapper/blob/master/LICENSE.txt)) | ✅ Yes — standard AGPL terms (include LICENSE + source pointer) |
| **EggNOG database** | ~50 GB of precomputed ortholog data the tool queries against | **Unspecified** | ⚠️ Ambiguous — see below |

This split (open-source code, separately-licensed data) is the common
bioinformatics pattern. Bakta does it (Apache code, Bakta DB), Prokka does it
(GPL code, mixed-licence bundled DBs), InterProScan does it (Apache core,
mixed-licence member DBs). Bundling the *code* is straightforward; bundling
the *data files* is the question we have to ask EMBL.

For the **database**: we checked the
[EggNOG website](http://eggnog5.embl.de/), the entire download tree at
`/download/`, the `eggnog-mapper` repo (which only carries the AGPL covering
the *code*), the EggNOG papers (Huerta-Cepas 2019 for v5.0,
Hernández-Plaza 2023 for v6.0), and the wiki. **No license clause is stated
anywhere for the data files** — only the AGPL-3.0 on the code repo and
third-party MIT notices for the website's bundled JS.

Under default copyright law (EU/Germany, EMBL-Heidelberg), absence of a
license on the data files means **all rights reserved** — silence is not
permission. Bakta, Prokka, and nf-core/funcscan all use
`download_eggnog_data.py` as the install path; no public Docker image
redistributes the ~50 GB database itself.

**Action:** v1.0.0 Docker image bundles the eggnog-mapper *code* (AGPL is
straightforward to comply with) but fetches the *database* via
`download_eggnog_data.py` at first install — same UX as IPS. In parallel,
the Billerbeck Lab has emailed `eggnog@embl.de` to request explicit
redistribution permission for the database files; if granted, a future
release can bundle the DB to make the Docker image fully self-contained.

## SignalP 6.0 + DeepLocPro — ⏳ Awaiting DTU response

Both released by DTU under an academic-only license. The license issued to
the Billerbeck Lab covers internal academic use. Whether ssign can include
the **binary + model weights** inside its public Docker image is currently
under discussion with DTU.

Two outcomes possible:

1. **DTU permits redistribution** — bundle inside Docker image, document
   academic-use restriction, users get a smooth `docker pull` experience.
2. **DTU requires per-user acquisition** — Docker image expects a bind-mount
   pointing at the user's own DTU-licensed install at runtime. Documented in
   `docs/how-to/install.md`.

Until DTU responds, the default **remote API mode** (via BioLib) requires no
license and is the current shipping path.

## DTU/SignalP — Phobius and TMHMM omitted from ssign

ssign's IPS configuration explicitly excludes Phobius, SignalP, and TMHMM
from InterProScan analyses (see `DEFAULT_IPS_APPLICATIONS` in
`run_interproscan.py`) — these member analyses each require their own
license and the IPS distribution itself ships pre-stripped of them. Not
relevant unless a user opts back in via `--applications`.

---

## Outstanding correspondence

External actions still in flight before v1.0.0 release (tracked under the
license-audit task in the publication roadmap):

- **DTU email** (Sonja Billerbeck → DTU) — SignalP 6.0 + DeepLocPro
  redistribution in our Docker image. Draft in roadmap.
- **EMBL-EggNOG email** (Sonja Billerbeck → eggnog@embl.de) — EggNOG database
  redistribution. Draft in roadmap.

Once both replies are in, this page gets updated and either the Docker image
gets fatter (bundling more) or the fetch-script story gets fleshed out.

---

## How license decisions feed into the install experience

Ssign's three-tier install (`pip install ssign[base|extended|full]`) maps onto
the licensing decisions above:

- The Docker image carries everything the user is *legally* allowed to receive
  from us in one pull.
- `scripts/fetch_databases.sh` (Phase 4b) handles the rest in a single
  `--tier` invocation that reads the tier name and downloads each missing DB
  from its canonical source.

Net user experience: `docker pull` → `ssign install-databases --tier extended`
→ `ssign run input.gbff`. Three commands, two of them once-only.
