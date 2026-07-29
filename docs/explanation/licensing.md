# Licensing: ssign dependencies and redistribution

This page documents the redistribution status of every external tool, model, and
database ssign uses.

ssign's own code is **GPL-3.0-or-later** (see `LICENSE`). The **container image**
is a separate matter: it bakes in DeepLocPro (CC BY-NC-SA 4.0, a non-commercial
licence) and EggNOG-mapper (AGPL-3.0, strong copyleft). So **the image as a
whole is for non-commercial / academic research use**, and anyone who
redistributes it must honour the AGPL source-offer (satisfied by the public
ssign + eggnog-mapper repositories). The rows below say whether the image can
include each binary / weight / database, or whether the user fetches it later.

| Component | License | In the container image? | If not, obtain via |
|---|---|---|---|
| **ssign** (this package) | GPL-3.0-or-later | Yes (it is the image) | n/a |
| **DeepLocPro** (tool + weights) | CC BY-NC-SA 4.0 (non-commercial) | Baked | n/a |
| **eggnog-mapper** (tool / code) | AGPL-3.0 | Baked (isolated env; makes the image AGPL-copyleft) | n/a |
| **DeepSecE** checkpoint (fine-tuned ESM-1b) | MIT | Baked | n/a |
| **ProtT5** encoder (Rostlab/prot_t5_xl_half_uniref50-enc) | AFL-3.0 | Baked | n/a |
| **ESM-1b / ESM2** backbones | MIT | Baked | n/a |
| **Bakta / HH-suite / MacSyFinder** (tools) | GPL-family | Baked | n/a |
| **TXSScan models** | CeCILL | Bundled with MacSyFinder | n/a |
| **SignalP 6.0** (binary + weights) | DTU academic | Not baked (cannot redistribute; DTU reply 2026-05-07) | `ssign-setup-dtu` (DTU portal), or `--signalp-mode remote` |
| **InterProScan** (engine + member DBs) | Apache-2.0 core + mixed members | Not baked (member DBs forbid redistribution) | `scripts/fetch_databases.sh` (public EBI, no licence gate) |
| **EggNOG database** (~47 GB) | unspecified | Not baked (size + licence) | `scripts/fetch_databases.sh` |
| **Bakta DB** | CC-BY 4.0 | Not baked (size) | `scripts/fetch_databases.sh` |
| **BLAST Swiss-Prot / NR** | NCBI public | Not baked (size) | `scripts/fetch_databases.sh` |
| **HH-suite DBs** (Pfam / PDB70 / UniRef30) | mixed (Söding lab + Tübingen) | Not baked (size) | `scripts/fetch_databases.sh` |

---

## DeepLocPro: baked, open source

DeepLocPro is open source
([github.com/Jaimomar99/deeplocpro](https://github.com/Jaimomar99/deeplocpro)),
released under **CC BY-NC-SA 4.0**, so ssign bakes both the tool and its weights
(plus the ESM2 backbone it needs) into the image at a pinned commit.

## EggNOG: code is AGPL-3.0 (baked), database is ambiguous (not baked)

These are two separate things:

|     | What it is | License | In the image? |
| --- | --- | --- | --- |
| **eggnog-mapper code** | ~5,000-line Python tool that runs DIAMOND and parses output | AGPL-3.0 ([`LICENSE.txt`](https://github.com/eggnogdb/eggnog-mapper/blob/master/LICENSE.txt)) | Baked (isolated env) |
| **EggNOG database** | ~47 GB of precomputed ortholog data the tool queries | **Unspecified** | Host-fetched |

The v1.0.0 image **bakes the eggnog-mapper code** (AGPL-3.0) into its own
micromamba environment. That isolation sidesteps a version conflict: eggnog-mapper
hard-pins `biopython==1.76` while Bakta needs `biopython>=1.78`, so the two live
in separate conda environments. Baking AGPL code is what gives the image its
AGPL source-offer (satisfied by the public repository).

The **database** is a different question. We checked the
[EggNOG website](http://eggnog5.embl.de/), the `/download/` tree, the
`eggnog-mapper` repo (whose AGPL covers only the code), the EggNOG papers
(Huerta-Cepas 2019, Hernández-Plaza 2023), and the wiki: **no license clause is
stated anywhere for the data files.** Under default copyright law (EU), silence means all rights reserved, so the ~47 GB database is
**not** baked. Users fetch it with `scripts/fetch_databases.sh` (which wraps
`download_eggnog_data.py`), the same install path Bakta, Prokka, and
nf-core/funcscan use. EggNOG annotation is on by default at the extended and full
tiers (off only at `--tier base`). The Billerbeck Lab has asked `eggnog@embl.de`
for permission to redistribute the database files, no response was received.

## ProtT5 weights: bundle OK

Used by pLM-BLAST. Released under the **Academic Free License 3.0** (AFL-3.0), an
OSI-approved permissive license
([ProtTrans README](https://github.com/agemagician/ProtTrans)) that grants public
redistribution with no non-commercial or field-of-use restriction and is
compatible with ssign's GPL-3.0. ColabFold and bio_embeddings already
redistribute these weights (established precedent); ssign archives them inside
the Zenodo-deposited image rather than as a separate deposit.

## InterProScan: not baked, user provides it

The core engine is **Apache 2.0** (redistribution-friendly), but the bundled
member databases have mixed licenses, and the show-stoppers are:

- **PROSITE** (Profiles + Patterns): SIB custom license, non-commercial, no
  derivatives; forbids redistribution outside the InterPro bundle.
- **SMART**: EMBL academic-only; "redistribution path is via the InterPro
  consortium" per the SMART FAQ.
- **SUPERFAMILY**: no formal license (flagged ambiguous by reusabledata.org).

**Industry signal:** EBI's own
[`interpro/interproscan` image](https://hub.docker.com/r/interpro/interproscan)
does *not* bundle the data either; no mainstream academic image redistributes
the IPS bundle.

**Action:** ssign does **not** bake InterProScan. Both its engine and member
databases are host-provided; the image ships only a Java runtime to run a
host-mounted IPS install. `scripts/fetch_databases.sh --tier extended` pulls the
IPS distribution from EBI FTP, one extra command run once.

Sources: [InterPro license](https://interpro-documentation.readthedocs.io/en/latest/license.html) ·
[HowToDownload](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html) ·
[`interpro/interproscan` image](https://hub.docker.com/r/interpro/interproscan) ·
[SMART FAQ](https://smart.embl.de/help/FAQ.shtml) ·
[SUPERFAMILY reuse](https://reusabledata.org/supfam.html)

## SignalP 6.0: cannot bundle (DTU reply 2026-05-07)

DTU confirmed by email on 2026-05-07 that the SignalP 6.0 license does not
permit redistribution, so ssign cannot bundle the binary or weights. It is the
**only predictor not baked into the image**. Install
SignalP locally from the [DTU HealthTech portal](https://services.healthtech.dtu.dk/)
(the `ssign-setup-dtu` helper walks through it) and wire it in with
`--signalp-path <dir>`. Without a DTU licence, the remote DTU webserver is an
opt-in fallback (`--signalp-mode remote`) needing no licence on the user's part.
By default ssign runs SignalP locally and never contacts DTU. See
`docs/how-to/install.md`.

---

## How this maps to install

The image carries every tool and weight ssign may *legally* redistribute; `scripts/fetch_databases.sh` pulls the databases too
large or not licence-clear to bundle, each from its canonical source.
