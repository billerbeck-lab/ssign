# `data/`: reference databases and raw inputs

This directory is a **placeholder**. The reference databases used by ssign are
too large (~1.5 TB combined) to distribute through Git, so they are hosted
externally and fetched at install time.

## What goes here at runtime

After a full install, this directory (or a path you pass to ssign via config)
will contain:

| Database                                 | Tier            | Size    | Purpose                                               |
| ---------------------------------------- | --------------- | ------- | ----------------------------------------------------- |
| NCBI taxdump (`nodes.dmp` + `names.dmp`) | base            | ~1.5 GB | Local taxonomy resolution (replaces NCBI E-utilities) |
| Bakta (light)                            | base            | ~2 GB   | Gene annotation (light variant)                       |
| Bakta (full)                             | full            | ~84 GB  | Whole-genome gene annotation                          |
| EggNOG-mapper                            | extended / full | ~50 GB  | Orthologous groups, COGs, KEGG                        |
| HH-suite (Pfam + PDB70)                  | extended        | ~26 GB  | Pfam (~3 GB) + PDB70 (~23 GB) profile HMM databases   |
| HH-suite (+ UniRef30)                    | full            | +~25 GB | Adds UniRef30 for deep remote-homology MSA generation |
| InterProScan                             | extended / full | ~24 GB  | Domain / motif annotation                             |
| pLM-BLAST (ECOD70)                       | extended / full | ~24 GB  | Embedding-based remote homology                       |
| BLAST NR                                 | full            | ~390 GB | Primary homology search                               |

## Getting the databases

ssign ships a tier-aware fetch script:

```bash
bash scripts/fetch_databases.sh --tier base       # ~3 GB
bash scripts/fetch_databases.sh --tier extended   # ~150 GB
bash scripts/fetch_databases.sh --tier full       # ~630 GB
```

The fetcher pulls each database from its canonical academic mirror (Bakta
GitHub release, Tübingen MPI for HH-suite + ECOD70, NCBI for taxdump and
BLAST NR, EBI for InterProScan, EMBL for EggNOG). At v1.0.0 release, a
mirror copy of the full set is also deposited on **Zenodo** under a pinned
DOI for long-term reproducibility; the fetch script gains a
`--source zenodo` flag at that time.

If a deposit exceeds Zenodo's per-record size limit, it is split across
multiple linked records, or hosted on Imperial's Research Data Store with
a persistent URL and checksum manifest.

## Integrity checking (post-v1.0.0)

Once the Zenodo mirror lands, every file will have an SHA-256 checksum
recorded alongside it. To verify after download:

```bash
cd data && sha256sum -c checksums.sha256
```

For now, integrity relies on HTTPS + tar's own corruption detection.

## What must **not** go here

- User-uploaded genomes or proteins for a specific run (those live in the
  user's chosen output directory).
- Any cleaned, preprocessed, or derived data (that belongs in `results/`).
- Model weights (those live in `models/`).
