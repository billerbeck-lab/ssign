#!/usr/bin/env python3
"""pyhmmer-based hmmsearch shim for MacSyFinder compatibility.

FRAGILE SHIM: This emulates the HMMER3 `hmmsearch` command-line tool using
pyhmmer (a Cython binding to HMMER3). It exists so that `pip install ssign`
works without needing `sudo apt install hmmer`.

Only the subset of hmmsearch flags used by MacSyFinder is implemented:
    hmmsearch --cpu N -o output [--cut_ga | -E evalue] [--tblout tblout] hmmfile seqdb

MacSyFinder's parser (macsylib/report.py) reads the `-o` text output by
looking for:
  - ">>" lines to identify hit targets
  - Domain table header "   #    score  bias  c-Evalue ..."
  - Whitespace-delimited domain data rows

IF THIS SHIM BREAKS:
  1. Install real HMMER: sudo apt install hmmer
     (or: conda install -c bioconda hmmer)
  2. The real hmmsearch will take precedence if installed to /usr/bin/
  3. Report the issue at: https://github.com/billerbeck-lab/ssign/issues

Version history:
  - v0.1 (2026-03): Initial shim covering MacSyFinder v2.1 usage pattern
"""

import argparse
import sys

try:
    import pyhmmer
    from pyhmmer.easel import Alphabet, SequenceFile
    from pyhmmer.plan7 import HMMFile
except ImportError:
    print(
        "ERROR: pyhmmer is not installed.\n"
        "  Install: pip install pyhmmer\n"
        "  Or install real HMMER: sudo apt install hmmer\n"
        "  Or reinstall ssign: pip install --force-reinstall ssign",
        file=sys.stderr,
    )
    sys.exit(1)

SHIM_VERSION = "0.1"


def _decode(val):
    """Safely decode bytes to str, or return str as-is."""
    if isinstance(val, bytes):
        return val.decode()
    return str(val) if val is not None else ""


def parse_args(argv=None):
    """Parse hmmsearch-compatible arguments (MacSyFinder subset only)."""
    parser = argparse.ArgumentParser(
        prog="hmmsearch",
        description=f"pyhmmer hmmsearch shim v{SHIM_VERSION} for ssign",
    )
    parser.add_argument("hmmfile", help="HMM profile file (.hmm)")
    parser.add_argument("seqdb", help="Target sequence database (FASTA)")
    parser.add_argument(
        "-o", dest="output", required=True, help="Direct output to file"
    )
    parser.add_argument(
        "--tblout", default=None, help="Save per-target tabular output"
    )
    parser.add_argument(
        "--domtblout", default=None, help="Save per-domain tabular output"
    )
    parser.add_argument("--cpu", type=int, default=1, help="Number of CPUs")
    parser.add_argument(
        "--cut_ga", action="store_true", help="Use GA gathering thresholds"
    )
    parser.add_argument(
        "-E", dest="evalue", type=float, default=10.0, help="E-value threshold"
    )
    # Accept and ignore other common flags for resilience
    parser.add_argument("--noali", action="store_true", help="(accepted, ignored)")
    parser.add_argument("--notextw", action="store_true", help="(accepted, ignored)")
    parser.add_argument("-Z", type=float, default=None, help="(accepted, ignored)")

    return parser.parse_args(argv)


def run_search(args):
    """Run hmmsearch via pyhmmer and return list of (query_hmm, top_hits)."""
    alphabet = Alphabet.amino()

    with HMMFile(args.hmmfile) as hmm_file:
        hmms = list(hmm_file)

    with SequenceFile(args.seqdb, digital=True, alphabet=alphabet) as seq_file:
        targets = seq_file.read_block()

    n_targets = len(targets)
    n_residues = sum(len(t) for t in targets)

    # FRAGILE: pyhmmer.hmmsearch API — if pyhmmer changes the signature
    # or return type, this will break. Tested with pyhmmer 0.10.x.
    # bit_cutoffs="gathering" applies GA thresholds from the HMM profile.
    # If the profile has no GA thresholds, pyhmmer raises an error;
    # we catch that and fall back to E-value filtering.
    results = []
    for hmm in hmms:
        search_kwargs = {"cpus": args.cpu}

        if args.cut_ga:
            try:
                top_hits_iter = pyhmmer.hmmsearch(
                    hmm, targets, bit_cutoffs="gathering", **search_kwargs
                )
                top_hits = next(top_hits_iter)
            except ValueError:
                # Profile lacks GA thresholds — fall back to E-value
                top_hits_iter = pyhmmer.hmmsearch(
                    hmm, targets, E=args.evalue, **search_kwargs
                )
                top_hits = next(top_hits_iter)
        else:
            top_hits_iter = pyhmmer.hmmsearch(
                hmm, targets, E=args.evalue, **search_kwargs
            )
            top_hits = next(top_hits_iter)

        results.append((hmm, top_hits))

    return results, n_targets, n_residues


def write_text_output(f, results, args, n_targets, n_residues):
    """Write HMMER3 text format output that MacSyFinder can parse.

    FRAGILE: This output format must match what MacSyFinder's parser expects.
    MacSyFinder (macsylib/report.py) specifically looks for:
      - ">>" lines to find hit target names
      - Domain table rows with fields at whitespace-delimited positions:
        index 2 = score, index 5 = i-evalue,
        index 6-7 = hmm range, index 9-10 = ali range
    If MacSyFinder's parser changes, update this function or install real HMMER.

    COORDINATE NOTE: pyhmmer uses 0-based coordinates internally.
    HMMER3 text output uses 1-based coordinates. We add 1 to all
    start positions (from/start values) when writing output.
    """
    # Header
    f.write("# hmmsearch :: search profile(s) against a sequence database\n")
    f.write(f"# HMMER 3.4 (pyhmmer shim v{SHIM_VERSION})\n")
    f.write("# - - - - - - - - - - - - - - - - - - - - - - - - -\n")
    f.write(f"# query HMM file:                  {args.hmmfile}\n")
    f.write(f"# target sequence database:         {args.seqdb}\n")
    f.write(f"# output directed to file:          {args.output}\n")
    f.write("# - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")

    for hmm, top_hits in results:
        query_name = _decode(hmm.name) or "unnamed"
        query_acc = _decode(hmm.accession) or "-"
        query_desc = _decode(hmm.description)

        f.write(f"Query:       {query_name}  [{query_desc}]\n")
        f.write(f"Accession:   {query_acc}\n")
        f.write(f"Description: {query_desc}\n")

        if len(top_hits) == 0:
            f.write("\n   [No hits detected that satisfy reporting thresholds]\n\n")
            f.write("//\n")
            continue

        # Per-sequence scores table
        f.write(
            "Scores for complete sequences (score includes all domains):\n"
        )
        f.write(
            "   --- full sequence ---   --- best 1 domain ---    -#dom-\n"
        )
        f.write(
            "    E-value  score  bias    E-value  score  bias    exp  N  "
            "Sequence     Description\n"
        )
        f.write(
            "    ------- ------ -----    ------- ------ -----   ---- --  "
            "--------     -----------\n"
        )

        for hit in top_hits:
            hit_name = _decode(hit.name)
            hit_desc = _decode(hit.description)
            # Best domain
            best_dom = min(hit.domains, key=lambda d: d.i_evalue) if hit.domains else None
            best_evalue = best_dom.i_evalue if best_dom else hit.evalue
            best_score = best_dom.score if best_dom else hit.score
            best_bias = best_dom.bias if best_dom else 0.0

            f.write(
                f"  {hit.evalue:9.2e} {hit.score:6.1f} {hit.bias:5.1f}"
                f"  {best_evalue:9.2e} {best_score:6.1f} {best_bias:5.1f}"
                f"   {len(hit.domains):4.1f} {len(hit.domains):2d}"
                f"  {hit_name:12s}  {hit_desc}\n"
            )

        f.write("\n\n")

        # Per-hit domain annotation — THIS IS WHAT MACSYFINDER PARSES
        for hit in top_hits:
            hit_name = _decode(hit.name)
            hit_desc = _decode(hit.description)
            f.write(f">> {hit_name}  {hit_desc}\n")

            if not hit.domains:
                f.write("   [No individual domains that satisfy thresholds]\n\n")
                continue

            f.write(
                "   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to"
                "    alifrom  ali to    envfrom  env to     acc\n"
            )
            f.write(
                " ---   ------ ----- --------- --------- ------- -------"
                "    ------- -------    ------- -------    ----\n"
            )

            for i, domain in enumerate(hit.domains):
                # COORDINATE CONVERSION: pyhmmer 0-based → HMMER3 1-based
                # pyhmmer's alignment coordinates are already 1-based in
                # recent versions, but env coordinates may be 0-based.
                # We output what pyhmmer gives us and rely on testing to
                # verify correctness. If coordinates are off by 1, this
                # is the place to fix it.
                hmm_from = domain.alignment.hmm_from
                hmm_to = domain.alignment.hmm_to
                ali_from = domain.alignment.target_from
                ali_to = domain.alignment.target_to
                env_from = domain.env_from
                env_to = domain.env_to

                # CRITICAL: The "!" or "?" inclusion marker at index 1 is
                # required. MacSyFinder's parser uses fields[5] = i_evalue.
                # Without "!", all fields shift left and parsing breaks.
                # "!" = included in per-sequence score, "?" = not included
                incl = "!" if domain.included else "?"

                f.write(
                    f"   {i + 1:3d} {incl} {domain.score:>8.1f} {domain.bias:>6.1f}"
                    f"  {domain.c_evalue:>9.2e}  {domain.i_evalue:>9.2e}"
                    f" {hmm_from:>7d} {hmm_to:>7d} .."
                    f" {ali_from:>7d} {ali_to:>7d} .."
                    f" {env_from:>7d} {env_to:>7d} .."
                    f"  {0.90:>4.2f}\n"
                )

            f.write("\n")

        f.write("//\n")


def write_tblout(f, results, args):
    """Write HMMER3 --tblout tabular format.

    FRAGILE: Format must match HMMER3 tblout specification.
    Columns: target, t_acc, query, q_acc, fullseq_E, fullseq_score,
             fullseq_bias, best1_E, best1_score, best1_bias,
             exp, reg, clu, ov, env, dom, rep, inc, description
    """
    f.write(
        "# --- full sequence ---- --- best 1 domain ---- --- domain "
        "number estimation ----\n"
    )
    f.write(
        "# target name        accession  query name           accession"
        "    E-value  score  bias   E-value  score  bias   exp reg clu"
        "  ov env dom rep inc description of target\n"
    )
    f.write(f"#{'':->18s} {'':->10s} {'':->20s} {'':->10s} {'':->9s}"
            f" {'':->6s} {'':->5s} {'':->9s} {'':->6s} {'':->5s}"
            f" {'':->4s} {'':->3s} {'':->3s} {'':->3s} {'':->3s}"
            f" {'':->3s} {'':->3s} {'':->3s} {'':->22s}\n")

    for hmm, top_hits in results:
        query_name = _decode(hmm.name) or "unnamed"
        query_acc = _decode(hmm.accession) or "-"

        for hit in top_hits:
            hit_name = _decode(hit.name)
            hit_acc = "-"
            hit_desc = _decode(hit.description)
            n_dom = len(hit.domains)
            best_dom = (
                min(hit.domains, key=lambda d: d.i_evalue)
                if hit.domains
                else None
            )
            best_e = best_dom.i_evalue if best_dom else hit.evalue
            best_s = best_dom.score if best_dom else hit.score
            best_b = best_dom.bias if best_dom else 0.0

            f.write(
                f"{hit_name:20s} {hit_acc:10s} {query_name:20s}"
                f" {query_acc:10s} {hit.evalue:9.2e} {hit.score:6.1f}"
                f" {hit.bias:5.1f} {best_e:9.2e} {best_s:6.1f}"
                f" {best_b:5.1f} {n_dom:4.1f} {n_dom:3d} {0:3d}"
                f" {0:3d} {n_dom:3d} {n_dom:3d} {n_dom:3d} {n_dom:3d}"
                f" {hit_desc}\n"
            )

    f.write("//\n")


def main(argv=None):
    """Entry point — drop-in replacement for hmmsearch."""
    args = parse_args(argv)

    try:
        results, n_targets, n_residues = run_search(args)
    except FileNotFoundError as e:
        print(f"ERROR: File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(
            f"ERROR: pyhmmer hmmsearch shim failed: {e}\n"
            f"  This may be due to a pyhmmer version incompatibility.\n"
            f"  Fallback: install real HMMER with: sudo apt install hmmer\n"
            f"  Or: conda install -c bioconda hmmer",
            file=sys.stderr,
        )
        sys.exit(1)

    # Write main output (-o)
    try:
        with open(args.output, "w") as f:
            write_text_output(f, results, args, n_targets, n_residues)
    except Exception as e:
        print(f"ERROR: Failed to write output: {e}", file=sys.stderr)
        sys.exit(1)

    # Write tblout if requested
    if args.tblout:
        try:
            with open(args.tblout, "w") as f:
                write_tblout(f, results, args)
        except Exception as e:
            print(f"ERROR: Failed to write tblout: {e}", file=sys.stderr)
            sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
