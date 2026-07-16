# Security policy

## Supported versions

ssign is in active development toward its first publication release (v1.0.0).
Until then, only the latest tagged release and the `main` branch receive
security fixes.

| Version | Supported |
| ------- | --------- |
| `main`  | ✅        |
| `v1.0.x` | ✅       |
| older   | ❌        |

## Reporting a vulnerability

**Please do not file public GitHub issues for security vulnerabilities.**

Report security issues privately using GitHub's
[private vulnerability reporting](https://github.com/billerbeck-lab/ssign/security/advisories/new),
or by email to:

- **M. Teo Reid** (maintainer): `t.reid25@imperial.ac.uk`
- **Dr. Sonja Billerbeck** (PI): `s.billerbeck@imperial.ac.uk`

Include:

- A description of the issue and the affected component (script, function,
  external dependency).
- Steps to reproduce, ideally with a minimal input.
- Your assessment of impact (e.g. arbitrary code execution, data leak,
  denial of service).

We aim to acknowledge reports within **14 days**, and to triage and either
fix or document mitigation within **60 days** for confirmed
medium-or-higher severity issues. Reporters are credited in the release
notes unless you ask otherwise.

## Scope

ssign is a research-software pipeline, not a hosted service. The primary
attack surfaces are:

- **Local input handling**: parsing user-supplied GenBank, FASTA, GFF3 files.
- **Subprocess invocations**: ssign shells out to external bioinformatics
  tools (BLAST, MacSyFinder, Bakta, etc.); injection via genome filenames
  or sequence headers is in-scope.
- **Network fetches**: `scripts/fetch_databases.sh`, the model-weight
  downloads each predictor makes on first run (torch-hub / Hugging Face hub),
  and the optional remote modes (DTU webserver for SignalP / DeepLocPro).

The hosted web service (planned, post-v1.0.0) will have its own threat
model, documented when it is built.

Out of scope: vulnerabilities in upstream tools (BLAST, MacSyFinder, etc.) —
please report those to their respective projects.
