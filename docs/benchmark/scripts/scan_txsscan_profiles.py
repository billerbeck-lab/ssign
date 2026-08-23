#!/usr/bin/env python3
"""Scan the 85 benchmark proteins with the TXSScan HMM profiles.

MacSyFinder has no training set, so the leakage question for it is different:
is the *secreted protein itself* modelled as a machinery component? For most
systems the answer is no (TXSScan models machinery). The exceptions matter:
autotransporters carry the translocator domain on the same polypeptide as the
passenger, and T6SS VgrG/Hcp are both structural components and secreted.

This answers it empirically rather than by argument, reproducing MacSyFinder's
own hit criteria as implemented in macsylib (`report.py:_parse_hmm_body`,
defaults in `config.py`), not the looser paraphrase in the 2016 paper:

  * `--cut_ga` is ON by default, so a profile carrying a GA cutoff is filtered
    by GA, not by an E-value. 278 of the 280 TXSScan profiles carry one; the
    remaining two fall back to `e_value_search` (0.1).
  * a domain is kept when `i_evalue <= 0.001` AND
    `(ali_to - ali_from + 1) / profile_length >= 0.5`. Both comparators are
    inclusive; using strict `<` and `>` drops boundary hits.
  * i-evalues scale with the size of the searched database. MacSyFinder searches
    a whole proteome, so scanning only these 85 sequences would make every
    i-evalue ~60x optimistic. Z is therefore fixed to a representative proteome
    size rather than left at the query count. GA is the binding constraint for
    all but two profiles, so the exact value barely matters.
"""

import glob
import os
import sys

import pyhmmer
from _paths import WORK

PROFILE_DIR = os.path.expanduser("~/.macsyfinder/models/TXSScan/profiles")
METADATA = os.path.expanduser("~/.macsyfinder/models/TXSScan/metadata.yml")
FASTA = os.path.join(WORK, "benchmark_proteins.fasta")
OUT = os.path.join(WORK, "txsscan_profile_hits.tsv")

I_EVALUE_MAX = 0.001  # macsylib config.i_evalue_sel
MIN_PROFILE_COV = 0.5  # macsylib config.coverage_profile
E_VALUE_SEARCH = 0.1  # macsylib config.e_value_search, used only without GA
PROTEOME_Z = 5000  # representative diderm proteome; see module docstring


def txsscan_version():
    if not os.path.exists(METADATA):
        return "unknown"
    for line in open(METADATA):
        if line.strip().startswith("vers:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


alphabet = pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(FASTA, digital=True, alphabet=alphabet) as sf:
    seqs = sf.read_block()

hmm_files = sorted(glob.glob(os.path.join(PROFILE_DIR, "*.hmm")))
if not hmm_files:
    sys.exit(f"no HMM profiles under {PROFILE_DIR}. Install them with:\n    macsydata install --user TXSScan==1.1.4")

# Seven pairs of TXSScan profiles share an internal HMM NAME (e.g. both
# T6SSii_iglB and T6SSiii_tssB are named "fam1.tab"), so a name-keyed lookup
# would mislabel hits. Rename each HMM to its file stem, which is the profile
# identifier the model XMLs actually reference.
hmms = []
for hf in hmm_files:
    with pyhmmer.plan7.HMMFile(hf) as f:
        for hmm in f:
            hmm.name = os.path.basename(hf)[:-4].encode()
            hmms.append(hmm)

with_ga = [h for h in hmms if h.cutoffs.gathering_available()]
without_ga = [h for h in hmms if not h.cutoffs.gathering_available()]
print(f"query proteins: {len(seqs)}", file=sys.stderr)
print(
    f"TXSScan v{txsscan_version()} profiles: {len(hmms)} ({len(with_ga)} with GA, {len(without_ga)} without)",
    file=sys.stderr,
)

rows = []
for batch, opts in ((with_ga, {"bit_cutoffs": "gathering"}), (without_ga, {"E": E_VALUE_SEARCH})):
    if not batch:
        continue
    for hits in pyhmmer.hmmsearch(batch, seqs, Z=PROTEOME_Z, **opts):
        M = hits.query.M
        pname = hits.query.name
        pname = pname if isinstance(pname, str) else pname.decode()
        for hit in hits:
            # strongest domain, not the first in alignment order
            dom = min(hit.domains, key=lambda d: d.i_evalue)
            ie = dom.i_evalue
            cov = (dom.alignment.hmm_to - dom.alignment.hmm_from + 1) / M
            if ie <= I_EVALUE_MAX and cov >= MIN_PROFILE_COV:
                name = hit.name if isinstance(hit.name, str) else hit.name.decode()
                parts = name.split("|")
                rows.append(
                    {
                        "instance_id": parts[0],
                        "uniprot": parts[1] if len(parts) > 1 else "",
                        "gene": parts[2] if len(parts) > 2 else "",
                        "ss_type": parts[3] if len(parts) > 3 else "",
                        "txsscan_profile": pname,
                        "i_evalue": f"{ie:.2e}",
                        "profile_coverage": round(cov, 3),
                        "bitscore": round(dom.score, 1),
                    }
                )

rows.sort(key=lambda r: (r["instance_id"], float(r["i_evalue"])))

cols = ["instance_id", "uniprot", "gene", "ss_type", "txsscan_profile", "i_evalue", "profile_coverage", "bitscore"]
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

prots = sorted({r["instance_id"] for r in rows})
print(f"\nproteins matched by >=1 TXSScan profile: {len(prots)}/{len(seqs)}", file=sys.stderr)
for r in rows:
    print(
        f"  {r['instance_id']:<10} {r['gene'][:20]:<20} {r['ss_type']:<6} -> "
        f"{r['txsscan_profile']:<28} iE={r['i_evalue']:<10} cov={r['profile_coverage']}",
        file=sys.stderr,
    )
print(f"\nwrote {OUT}", file=sys.stderr)
