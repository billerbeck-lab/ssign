# containers/

Build artifacts for the **ssign** Singularity/Apptainer image. To **run** ssign,
start at [`docs/how-to/install.md`](../docs/how-to/install.md)
(install, fetch databases, run). This directory is for maintainers who **build and
publish** the image.

Contents:
- **`ssign.def`** the image definition (digest-pinned CUDA base + frozen `uv.lock`).
- **`build_sif.sh`** build + offline smoke-test helper.
- **`requirements-base.lock.txt`** pinned deps for the macOS/base pip path.
- The launcher scripts users run (`ssign-run`, `ssign-setup-dtu`) live in `scripts/`
  and are baked into the image at `/opt/ssign/scripts/`.

The image bundles DeepLocPro (CC BY-NC-SA) and EggNOG-mapper (AGPL), so it is for
**non-commercial** research use.

## Build the image

Build on a machine with **open internet** (laptop, workstation, or CI). An HPC
compute node usually **can't** build it: clusters allowlist package hosts and block
the general Ubuntu servers the base image's `apt` needs. Build once, copy the `.sif`
up.

```bash
# needs apptainer + a real-disk build tmp (NOT a small tmpfs). ~19 GB image, ~30-60 min.
containers/build_sif.sh          # build + offline smoke test
# ...or by hand:
export APPTAINER_TMPDIR=$HOME/.ssign_build/tmp APPTAINER_CACHEDIR=$HOME/.ssign_build/cache TMPDIR=$APPTAINER_TMPDIR
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
apptainer build --fakeroot ssign.sif containers/ssign.def
```

- The build downloads ~10 GB of model weights (ESM2, ESM-1b (the DeepSecE backbone),
  the DeepSecE checkpoint from a slow academic mirror, ProtT5). Keep the connection up: a container build has no
  resume, so a dropped download restarts that file from 0 and can abort the build.
  The DeepSecE mirror is the slow part; everyone who pulls the finished image skips it.
- `--fakeroot` builds unprivileged (standard `/etc/subuid` mapping; no system root).
- The base is digest-pinned; re-pin per the `ssign.def` comment only if a base CVE
  forces it.

## Publish the image (release)

After this, others get the image with one `apptainer pull` (GHCR) or a download
(Zenodo), then run per
[`install.md`](../docs/how-to/install.md).

**Gate:** publish ONLY after a tier-2 (extended) container run passes the full
extended-tier step count on a real genome. A Zenodo DOI is permanent; publishing a
broken image is not undoable.

You need:
- The validated `.sif` (referred to below as `$SIF`).
- A **Zenodo** account + a personal access token with `deposit:write` +
  `deposit:actions`: <https://zenodo.org/account/settings/applications/tokens/new/>.
- **(Optional, GHCR)** a GitHub personal access token (classic) with `write:packages`.

```bash
SIF=/path/to/validated.sif      # the image that passed the full extended-tier run
VERSION=1.0.0
```

### 1. Zenodo (primary: archival, permanent DOI, no pull-time auth)

Zenodo gives the citable DOI the paper points at, allows 50 GB per file, and needs no
login to download. The web uploader chokes on a large file; use the API (resumable).
With `ZENODO_TOKEN` set:

```bash
Z=https://zenodo.org/api/deposit/depositions
dep=$(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" \
      -X POST "$Z" -d '{}')
bucket=$(printf '%s' "$dep" | python3 -c 'import sys,json;print(json.load(sys.stdin)["links"]["bucket"])')
curl -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file "$SIF" "$bucket/ssign_$VERSION.sif"
```

Then in the Zenodo web UI for that deposition fill title (`ssign v$VERSION
(Singularity image)`), authors, description, and a licence note (the image bundles
DeepLocPro CC BY-NC-SA 4.0 and EggNOG-mapper AGPL, so non-commercial research use;
ssign's own code is GPL-3.0-or-later). Publish, then copy the DOI. A future image
version becomes a new version under the same **concept DOI** (cite the concept DOI so
`latest` always resolves).

### 2. GHCR (optional: convenience `apptainer pull`)

Nicer UX, but a large single blob can exceed a registry's per-layer limit; if the push
errors on size, skip it and ship via Zenodo.

```bash
echo "$GH_PAT" | apptainer registry login -u <gh-user> --password-stdin oras://ghcr.io
apptainer push "$SIF" oras://ghcr.io/billerbeck-lab/ssign:$VERSION
# make the package public: GitHub org -> Packages -> ssign -> Package settings -> visibility
apptainer pull oras://ghcr.io/billerbeck-lab/ssign:$VERSION   # verify a clean, unauthenticated pull
```

### 3. Wire the docs + tag the release

- Point the `apptainer pull` line in `README.md` and
  `docs/how-to/install.md` at the source that actually worked (GHCR if the
  push took, otherwise `download from Zenodo: <DOI>`).
- Add the DOI to the README badge + `CITATION.cff`.
- Tag it: `git tag v$VERSION && git push origin v$VERSION`.

### Notes

- The `.sif` is self-contained except SignalP 6 (DTU licence) and the reference
  databases; those stay user-provided (`fetch-databases` handles the DBs).
- Re-publishing a fixed image: new GHCR tag (`:1.0.1`) + a new Zenodo version on the
  same concept DOI. Never overwrite a published Zenodo file.
