# Publishing the ssign image (release)

Turns the `super easy install` on: after this, others get the image with one
`apptainer pull` (GHCR) or a download (Zenodo), then run per `containers/README.md`.

**Gate:** do this ONLY after a tier-2 (extended) container run passes **23/23** on
a real genome. A Zenodo DOI is permanent; publishing a broken image is not undoable.

## What you need

- The validated `.sif` (e.g. `~/ssign_v5.sif`, ~13 GB).
- **Zenodo** account + a personal access token with `deposit:write` +
  `deposit:actions` scopes: <https://zenodo.org/account/settings/applications/tokens/new/>.
- **(Optional, GHCR)** a GitHub personal access token (classic) with `write:packages`.

## 1. Zenodo (primary: archival, permanent DOI, no pull-time auth)

Zenodo is the canonical host. It gives the citable DOI the paper points at, allows
50 GB per file, and needs no login to download.

The web uploader chokes on 13 GB; use the API (resumable). With `ZENODO_TOKEN` set:

```bash
Z=https://zenodo.org/api/deposit/depositions
# create an empty deposition
dep=$(curl -s -H "Authorization: Bearer $ZENODO_TOKEN" -H "Content-Type: application/json" \
      -X POST "$Z" -d '{}')
bucket=$(printf '%s' "$dep" | python3 -c 'import sys,json;print(json.load(sys.stdin)["links"]["bucket"])')
echo "$dep" | python3 -c 'import sys,json;print("deposition id:", json.load(sys.stdin)["id"])'
# upload the image (PUT straight to the bucket; resumable)
curl -H "Authorization: Bearer $ZENODO_TOKEN" --upload-file ~/ssign_v5.sif "$bucket/ssign_v1.0.0.sif"
```

Then in the Zenodo web UI for that deposition: fill title (`ssign v1.0.0
(Singularity image)`), authors, description, and a licence note, the image
**bundles DeepLocPro (CC BY-NC-SA 4.0) and EggNOG-mapper (AGPL), so it is for
non-commercial research use**; ssign's own code is GPL-3.0-or-later. Publish, then
copy the DOI. A future image version becomes a new version under the same concept
DOI (cite the concept DOI so `latest` always resolves).

## 2. GHCR (optional: convenience `apptainer pull`)

Nicer UX, but a 13 GB single blob can exceed a registry's per-layer limit; if the
push errors on size, skip it and ship via Zenodo.

```bash
echo "$GH_PAT" | apptainer registry login -u <your-gh-user> --password-stdin oras://ghcr.io
apptainer push ~/ssign_v5.sif oras://ghcr.io/billerbeck-lab/ssign:1.0.0
# make the package public: GitHub org -> Packages -> ssign -> Package settings -> Change visibility
apptainer pull oras://ghcr.io/billerbeck-lab/ssign:1.0.0    # verify a clean, unauthenticated pull
```

## 3. Wire the docs + tag the release

- Point the `apptainer pull` line in `README.md` and `containers/README.md` at the
  source that actually worked (GHCR if the push took, otherwise `download from
  Zenodo: <DOI>`).
- Add the DOI to the README badge + `CITATION.cff`.
- Tag it: `git tag v1.0.0 && git push origin v1.0.0`.

## Notes

- The `.sif` is self-contained except SignalP 6 (DTU licence) and the reference
  databases; those stay user-provided, `fetch-databases` handles the DBs. Nothing
  in this step changes that.
- Re-publishing a fixed image: new GHCR tag (`:1.0.1`) + a new Zenodo version on
  the same concept DOI. Never overwrite a published Zenodo file.
