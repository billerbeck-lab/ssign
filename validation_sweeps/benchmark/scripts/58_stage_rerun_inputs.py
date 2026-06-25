#!/usr/bin/env python3
"""Assemble one ssign input per distinct genome for the tier-2 rerun, from rerun_panel_manifest.tsv.

For each genome in the panel we produce exactly one input `.gbff` (whole-assembly), so MacSyFinder and
the whole-genome DLP/DSE enrichment see the complete genome:

  - single staged unit       -> the whole assembly: prefer inputs_gb_fullasm/<unit>.gbff (the parent
                                assembly script 50 staged when the corpus accession is only a
                                plasmid/WGS contig), else reuse inputs_gb/<unit>.gbff
  - split assembly (a group   -> if the units are the SAME sequence (equal total length = a RefSeq/INSDC
    with >1 unit)               duplicate, e.g. PAO1 AE004091==NC_002516.2) keep ONE (the unit with more
                                 answer-key effectors); if they are DIFFERENT replicons (e.g. A. fabrum
                                 C58 linear chromosome + Ti plasmid) concatenate them into one input.
  - cache-not-staged T5SS     -> concatenate the genome's cached replicon .gb files into a new gbff.

Drops the 15 system-less units and the duplicate copies. Emits the submit list (one input filename per
distinct genome) + reports which files are NEW (must be copied to CX3; the rest are already on CX3 from
the prior fleet).

Inputs : data/phase2/rerun_panel_manifest.tsv, data/phase2/_t5_acc_cache.json, inputs_gb/*.gbff,
         data/refseq_cache/*.gb
Outputs: inputs_gb/<merged-or-staged>.gbff (new whole-assembly inputs), data/phase2/rerun_inputs.txt
Run    : <repo>/.venv/bin/python scripts/58_stage_rerun_inputs.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_index import accession_base as norm  # noqa: E402
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
MANIFEST = P2 / "rerun_panel_manifest.tsv"
ACC_CACHE = P2 / "_t5_acc_cache.json"
INPUTS = BENCH / "inputs_gb"
FULLASM = BENCH / "inputs_gb_fullasm"  # script 50: parent assembly for fragment-only corpus accessions
CACHE = BENCH / "data" / "refseq_cache"
OUT_LIST = P2 / "rerun_inputs.txt"


def stage_whole_assembly(genome: str, new_files: list[str]) -> str:
    """Filename for a single staged accession, preferring the script-50 full assembly when present.

    When the corpus refseq_genome is only a plasmid/WGS contig (e.g. hlyA on a 175-CDS plasmid),
    script 50 stages the PARENT whole assembly in inputs_gb_fullasm/<acc>.gbff. The rerun must run
    that, not the fragment, or MacSyFinder never sees the chromosomal machinery and the effector
    cannot emit. Copy it over the fragment so the single inputs_gb dir is always whole-assembly,
    and flag it for re-transfer to CX3."""
    name = f"{genome}.gbff"
    fa = FULLASM / name
    dst = INPUTS / name
    if fa.exists() and (not dst.exists() or dst.stat().st_size != fa.stat().st_size):
        dst.write_bytes(fa.read_bytes())
        new_files.append(name)
    return name


def gb_total_len(path: Path) -> int:
    """Sum of the `LOCUS ... <n> bp` lengths in a GenBank file (sequence-identity proxy)."""
    total = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                parts = line.split()
                for i, tok in enumerate(parts):
                    if tok == "bp" and i:
                        total += int(parts[i - 1])
    return total


def write_concat(dst: Path, srcs: list[Path]) -> None:
    with open(dst, "wb") as out:
        for s in srcs:
            out.write(s.read_bytes())


def main() -> int:
    rows = read_tsv(MANIFEST)
    run = [r for r in rows if r["action"] == "run"]
    add = [r for r in rows if r["action"] == "stage-from-cache"]

    inputs: list[str] = []  # input filename per distinct genome
    new_files: list[str] = []  # newly created -> need transfer to CX3

    # staged genomes: one input per genome_group
    by_group: dict[str, list[dict]] = defaultdict(list)
    for r in run:
        by_group[r["genome_group"]].append(r)
    for group, units in sorted(by_group.items()):
        if len(units) == 1:
            inputs.append(stage_whole_assembly(units[0]["genome"], new_files))
            continue
        # split assembly: cluster units by total sequence length (equal length = same sequence)
        by_len: dict[int, list[dict]] = defaultdict(list)
        for u in units:
            by_len[gb_total_len(INPUTS / f"{u['genome']}.gbff")].append(u)
        reps = [max(us, key=lambda u: int(u["n_secreted_proteins"])) for us in by_len.values()]
        if len(reps) == 1:  # pure duplicate (e.g. PAO1 RefSeq==INSDC) -> keep the richer copy
            inputs.append(stage_whole_assembly(reps[0]["genome"], new_files))
        else:  # genuine multi-replicon (e.g. C58) -> concatenate into one whole-assembly input
            reps.sort(key=lambda u: u["genome"])
            name = "_".join(u["genome"] for u in reps) + ".gbff"
            write_concat(INPUTS / name, [INPUTS / f"{u['genome']}.gbff" for u in reps])
            inputs.append(name)
            new_files.append(name)

    # cache-not-staged T5SS genomes: concat their cached replicon files into a new input
    acc_cache = json.loads(ACC_CACHE.read_text())
    cache_by_base = {norm(p.stem): p for p in CACHE.glob("*.gb")}
    for r in add:
        g = r["genome"]
        srcs = [cache_by_base[b] for b in acc_cache.get(g, [norm(g)]) if b in cache_by_base]
        if not srcs:
            print(f"  WARN no cached replicon for {g}; skipped")
            continue
        name = "_".join(sorted(p.stem for p in srcs)) + ".gbff"
        write_concat(INPUTS / name, sorted(srcs))
        inputs.append(name)
        new_files.append(name)

    inputs = sorted(set(inputs))
    OUT_LIST.write_text("\n".join(inputs) + "\n")

    print(f"distinct genomes in rerun : {len(inputs)}")
    print(f"  reused (already on CX3)  : {len(inputs) - len(new_files)}")
    print(f"  NEW (copy to CX3)        : {len(new_files)}")
    print(f"\nwrote submit list -> {OUT_LIST.relative_to(BENCH)} ({len(inputs)} inputs)")
    print("\nNEW input files to transfer to CX3:")
    for n in sorted(new_files):
        print(f"  inputs_gb/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
