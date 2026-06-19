#!/usr/bin/env python3
"""Regenerate each fleet genome's SS-component neighborhoods (the ±3 gene mask).

The circular-permutation null in 02 compares a system's secreted density to
random windows anywhere in the genome — including windows that land on OTHER
systems' substrate-rich neighborhoods. That makes the null too dense, so the
test is a conservative floor. To fix it we need to know which genes belong to
ANY system's neighborhood, then build the null only from non-neighborhood
windows. The fleet's pulled results don't carry that membership (nearby_ss_types
holds only the final substrate candidates, and the run's neighborhood files were
in purged job scratch), so we recompute it.

This reproduces the cheap detection chain locally (no GPU) for all 67 genomes,
matching the fleet config (TXSScan all, excluded Flagellum,Tad -> T3SS kept,
wholeness 0.8, window 3), and caches the union of all neighborhood locus_tags
per genome. ~20 s/genome of macsyfinder, run once.

    .venv/bin/python 03a_regen_neighborhoods.py            # all 67 (cached)
    .venv/bin/python 03a_regen_neighborhoods.py AE002098   # one genome (test)
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SCRIPTS = os.path.join(REPO, "src", "ssign_app", "scripts")
PY = os.path.join(REPO, ".venv", "bin", "python")
VENV_BIN = os.path.join(REPO, ".venv", "bin")
GBDIR = os.path.join(REPO, "validation_sweeps", "benchmark", "inputs_gb")
FLEET = "/tmp/ssign_fleet_67"
CACHE = os.path.join(HERE, "neighborhoods")  # one {genome}.json per genome
WORKROOT = os.path.join(HERE, "_regen_work")

WHOLENESS = "0.8"
EXCLUDED = "Flagellum,Tad"  # fleet ran full_t3ss -> T3SS kept
WINDOW = "3"


def run(script: str, args: list[str]) -> None:
    cmd = [PY, os.path.join(SCRIPTS, script), *args]
    r = subprocess.run(cmd, cwd=SCRIPTS, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"{script} failed (rc={r.returncode}):\n{r.stderr[-1500:]}")


def run_macsyfinder(proteins: str, out_dir: str) -> None:
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    # macsyfinder + the hmmsearch shim both live in venv/bin, off the inherited PATH.
    env = {**os.environ, "PATH": VENV_BIN + os.pathsep + os.environ.get("PATH", "")}
    cmd = [
        os.path.join(VENV_BIN, "macsyfinder"),
        "--sequence-db",
        proteins,
        "--db-type",
        "ordered_replicon",
        "--models",
        "TXSScan",
        "all",
        "--out-dir",
        out_dir,
        "-w",
        "8",
        "--mute",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
    if r.returncode != 0:
        raise SystemExit(f"macsyfinder failed (rc={r.returncode}):\n{r.stderr[-1500:]}")


def gbff_for(genome: str) -> str:
    # exact stem or stem + version suffix only; a loose startswith would let
    # AE002098 collide with a hypothetical AE0020980.gbff.
    for f in os.listdir(GBDIR):
        if f.endswith((".gbff", ".gb", ".gbk")) and (f == genome or f.startswith(genome + ".")):
            return os.path.join(GBDIR, f)
    raise SystemExit(f"no gbff for {genome}")


def regen_one(genome: str) -> dict:
    """Run detection -> return {neigh_ids:[...], n_systems, n_components}."""
    work = os.path.join(WORKROOT, genome)
    os.makedirs(work, exist_ok=True)
    proteins = os.path.join(work, "proteins.faa")
    gene_info = os.path.join(work, "gene_info.tsv")
    gene_order = os.path.join(work, "gene_order.tsv")
    msf_out = os.path.join(work, "macsyfinder_out")
    ss_components = os.path.join(work, "ss_components.tsv")
    valid_systems = os.path.join(work, "valid_systems.tsv")
    neigh_fasta = os.path.join(work, "neighborhood.faa")
    neigh_ids = os.path.join(work, "neighborhood_ids.tsv")

    run(
        "extract_proteins.py",
        ["--input", gbff_for(genome), "--sample", genome, "--out-proteins", proteins, "--out-gene-info", gene_info],
    )
    run("extract_gene_order.py", ["--gene-info", gene_info, "--output", gene_order])
    run_macsyfinder(proteins, msf_out)
    run(
        "validate_macsyfinder_systems.py",
        [
            "--msf-dir",
            msf_out,
            "--gene-info",
            gene_info,
            "--sample",
            genome,
            "--wholeness-threshold",
            WHOLENESS,
            "--excluded-systems",
            EXCLUDED,
            "--out-components",
            ss_components,
            "--out-systems",
            valid_systems,
        ],
    )
    run(
        "extract_neighborhood.py",
        [
            "--gene-order",
            gene_order,
            "--ss-components",
            ss_components,
            "--proteins",
            proteins,
            "--window",
            WINDOW,
            "--output",
            neigh_fasta,
            "--output-ids",
            neigh_ids,
        ],
    )

    ids = sorted({line.strip() for line in open(neigh_ids) if line.strip()})
    n_systems = (
        sum(1 for _ in csv.DictReader(open(valid_systems), delimiter="\t")) if os.path.exists(valid_systems) else 0
    )
    n_components = sum(1 for _ in csv.DictReader(open(ss_components), delimiter="\t"))
    return {"genome": genome, "neigh_ids": ids, "n_systems": n_systems, "n_components": n_components}


def main():
    os.makedirs(CACHE, exist_ok=True)
    only = sys.argv[1:] or sorted(os.listdir(FLEET))
    done = skipped = 0
    for g in only:
        cache_path = os.path.join(CACHE, f"{g}.json")
        if os.path.exists(cache_path):
            skipped += 1
            continue
        res = regen_one(g)
        with open(cache_path, "w") as fh:
            json.dump(res, fh)
        done += 1
        print(
            f"  {g}: {res['n_systems']} systems, {res['n_components']} components, "
            f"{len(res['neigh_ids'])} neighborhood loci"
        )
    print(f"\nregenerated {done}, cached-skip {skipped}; cache -> {CACHE}")
    if os.path.isdir(WORKROOT):
        shutil.rmtree(WORKROOT)  # scratch; the json cache is the durable output


if __name__ == "__main__":
    main()
