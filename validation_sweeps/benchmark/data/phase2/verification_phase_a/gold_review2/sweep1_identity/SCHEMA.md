# Sweep 1 — IDENTITY review (gold list second pass)

You are auditing the **identity fields** of a curated gold list of bacterial secreted-protein effectors.
Each row is one secretion-system instance and its single representative substrate. A deterministic
precompute already cleared every row whose UniProt accession matched its genome coordinates; the rows you
are given are ONLY the ones it could not clear:

- `no_acc` — the row has **no UniProt accession** (`uniprot_in` is `-`/blank). Your job: find the UniProt
  accession for **this exact protein in this exact genome/strain**, or determine that none exists.
- `len_mismatch` — the cited accession's length disagrees with the gene's CDS span at the given
  coordinates (`len_ratio` far from 1.0). Your job: decide whether the cited accession is actually the
  right protein for this gene at these coordinates, or whether a different accession fits better.

You are GIVEN, per row: genome accession + contig + start/stop/strand, the RefSeq-style `effector_locus`,
organism, gene name, the expected protein length in aa (`expected_aa`, = CDS span/3 − 1), and the
precomputed `len_ratio`/`signal_flags`. Use them. Do not recompute coordinates; adjudicate from them.

## How to find the genome-matched accession (do this, don't guess)

Query the UniProt REST API (via WebFetch or curl). Reliable patterns:

- By locus tag: `https://rest.uniprot.org/uniprotkb/search?query=%22<effector_locus>%22&fields=accession,reviewed,length,gene_names,organism_name,xref_refseq&format=tsv`
- By gene + organism: `https://rest.uniprot.org/uniprotkb/search?query=gene:<gene>+AND+organism_name:%22<organism words>%22&fields=accession,reviewed,length,gene_names,organism_name&format=tsv`
- Confirm a candidate: `https://rest.uniprot.org/uniprotkb/<acc>.txt` (read OS=, OX=, GN=, DR RefSeq, SQ length)

A candidate is a **match** only if: (a) its organism/strain is the row's organism (or a documented synonym
of that exact strain), AND (b) its length is within ~10% of `expected_aa` (small signal-peptide/processing
differences are fine; a 2× or ½× difference is NOT the same protein). Prefer a **reviewed (Swiss-Prot)**
entry when both reviewed and TrEMBL exist for the same protein.

## Output — ONE TSV row per input row, columns in this exact order

`row_id  gene  uniprot_in  ss_type  identity_verdict  proposed_uniprot  proposed_reviewed  organism_match  proposed_len_aa  len_consistent  evidence_source  confidence  notes`

| column | values | rule |
|---|---|---|
| `row_id` | echoed | MUST equal the input verbatim |
| `gene` | echoed | verbatim |
| `uniprot_in` | echoed | verbatim (`-`/blank as-is) |
| `ss_type` | echoed | verbatim |
| `identity_verdict` | `ok` \| `replace` \| `add` \| `none_exists` \| `unclear` | see below |
| `proposed_uniprot` | accession or blank | the accession to use; blank for `ok`/`none_exists`/`unclear` |
| `proposed_reviewed` | `yes` \| `no` \| blank | is `proposed_uniprot` a Swiss-Prot entry |
| `organism_match` | `yes` \| `strain_diff` \| `no` \| `unclear` | does the proposed (or cited) accession's organism match the row |
| `proposed_len_aa` | integer or blank | length of the proposed (or confirmed) accession |
| `len_consistent` | `yes` \| `no` \| `unclear` | is that length within ~10% of `expected_aa` |
| `evidence_source` | UniProt acc / RefSeq / blank | what you actually retrieved |
| `confidence` | `high` \| `medium` \| `low` | |
| `notes` | free text | name every accession + query you checked |

### Verdicts
- **ok**: the cited `uniprot_in` is correct for this gene+organism at these coordinates (use this to clear a
  `len_mismatch` you judge benign, e.g. a propeptide difference). `proposed_uniprot` blank.
- **replace**: a different accession fits the gene+organism+coordinates better than `uniprot_in`. Put it in
  `proposed_uniprot` with its length/organism/reviewed status.
- **add**: the row had no accession and you found a genome-matched one → `proposed_uniprot`.
- **none_exists**: no UniProt entry maps to this exact protein in this strain (blank accession is correct).
  Say what you searched and why nothing qualifies.
- **unclear**: you could not resolve it (no access, ambiguous). Do NOT guess an accession.

## Anti-hallucination rules (non-negotiable)
1. Echo `row_id`/`gene`/`uniprot_in`/`ss_type` exactly. Never alter them.
2. NEVER invent a UniProt accession. Every `proposed_uniprot` must be one you actually retrieved from
   UniProt and whose length + organism you report. No qualifying entry found → `none_exists` or `unclear`.
3. Report the proposed entry's real length and reviewed status from UniProt, not an estimate.
4. A different-strain accession of the same gene is `strain_diff`, not a match — only propose it if no exact
   strain entry exists AND you flag it; prefer `none_exists` over a wrong-strain accession.

Return ONLY the TSV (header + one row per input row), nothing else, then write it to the path given in your
task prompt.
