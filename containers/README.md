# containers/

Build artifacts for the **ssign** Singularity/Apptainer image. To **run** ssign,
start at [`docs/how-to/install.md`](../docs/how-to/install.md)
(install, fetch databases, run). This directory is for maintainers who **build and
publish** the image.

Contents:
- **`ssign.def`** the image definition (digest-pinned CUDA base + frozen `uv.lock`).
- **`build_sif.sh`** build + test helper.
- **`requirements-base.lock.txt`** pinned deps for the macOS/base from-source install.
- The launcher scripts users run (`ssign-run`, `ssign-setup-dtu`) live in `scripts/`
  and are baked into the image at `/opt/ssign/scripts/`.

The image bundles DeepLocPro (CC BY-NC-SA) and EggNOG-mapper (AGPL), so it is for
**non-commercial** research use.

## Build the image

Build on a machine with **open internet** (laptop, workstation).

```bash
# needs apptainer + a real-disk build tmp (not a small tmpfs). ~19 GB image, ~30-60 min.
containers/build_sif.sh          # build + test
# ...or by hand:
export APPTAINER_TMPDIR=$HOME/.ssign_build/tmp APPTAINER_CACHEDIR=$HOME/.ssign_build/cache TMPDIR=$APPTAINER_TMPDIR
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
apptainer build --fakeroot ssign.sif containers/ssign.def
```

- The build downloads ~15 GB of model weights (ESM2, ESM-1b (the DeepSecE backbone),
  the DeepSecE checkpoint from a slow mirror, ProtT5). Keep the connection up: a container build has no
  resume, so a dropped download restarts that file from 0 and can abort the build.
- `--fakeroot` builds unprivileged (standard `/etc/subuid` mapping; no system root).

## Publish the image (release)

**Gate:** publish ONLY after a tier-2 (extended) container run passes the full
extended-tier step count on a real genome.

You need:
- The validated `.sif` (referred to below as `$SIF`).
- A **Zenodo** account + a personal access token with `deposit:write` +
  `deposit:actions`: <https://zenodo.org/account/settings/applications/tokens/new/>.

```bash
SIF=/path/to/validated.sif      # the image that passed the full extended-tier run
VERSION=1.0.0
```

### 1. Zenodo

Zenodo gives the citable DOI the paper points at, allows 50 GB per file, and needs no
login to download. The web uploader sometimes chokes on a large file; use the API (resumable).
With `ZENODO_TOKEN` set:

```bash
Z=https://zenodo.org/api/deposit/depositions
dep=$(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" \
      -X POST "$Z" -d '{}')
bucket=$(printf '%s' "$dep" | python3 -c 'import sys,json;print(json.load(sys.stdin)["links"]["bucket"])')
curl -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file "$SIF" "$bucket/ssign_$VERSION.sif"
```

Then in the Zenodo web UI for that deposition fill title (`ssign v$VERSION
(Singularity image)`), authors, description, and a licence note. Publish, then copy the DOI.

### 2. Wire the docs + tag the release

- Add the DOI to the README badge + `CITATION.cff`.
- Tag it: `git tag v$VERSION && git push origin v$VERSION`.
