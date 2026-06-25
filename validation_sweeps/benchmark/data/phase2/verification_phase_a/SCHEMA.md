# Phase A re-verification: output schema + agent prompt (benchmark-final-validation 1.2)

Fresh, complete re-verification of all 337 answer-key effectors (`phasea_ledger.tsv`). Method reuses
the tasks-74-84 blind-agent design, hardened against hallucination: a verifier may only **echo** the
identifiers it was given and must **quote real provenance** for any verdict. It may never invent a DOI,
PMID, accession, or quote. If it cannot find support, it says so (`unclear`), it does not guess.

## Per-effector output (one row per ledger row)

Tab-separated, columns in this order:

| column | values | rule |
|---|---|---|
| `row_id` | echoed PA### | MUST equal the input verbatim (drift = discard the row) |
| `gene` | echoed | MUST equal the input verbatim |
| `uniprot` | echoed | MUST equal the input verbatim (`-`/blank echoed as-is) |
| `ss_type` | echoed | MUST equal the input verbatim |
| `verdict` | `confirmed` \| `fix` \| `drop` | see below |
| `defect` | controlled vocab or blank | required when verdict != confirmed |
| `secretion_claim_supported` | `yes` \| `no` \| `unclear` | does a REAL resolvable primary ref support "this protein is a secreted substrate/effector of this SS type"? |
| `doi_resolves` | `yes` \| `no` \| `paywalled` | the ledger `primary_ref` DOI/PMID |
| `uniprot_ok` | `yes` \| `no` \| `na` | accession live + matches gene+organism (`na` if blank/`-`) |
| `provenance_quote` | verbatim ≤300 chars | REQUIRED for `confirmed`/`fix`; a real sentence from the source, not paraphrased |
| `provenance_source` | PMID / PMCID / DOI / UniProt acc | REQUIRED whenever a quote is given; the exact id the quote is from |
| `new_primary_ref` | DOI/PMID or blank | only for `fix`; a REAL ref the verifier located |
| `new_uniprot` | accession or blank | only for `fix` |
| `new_locus_tag` | locus or blank | only for `fix` |
| `confidence` | `high` \| `medium` \| `low` | |
| `notes` | free text | the reasoning; name every id checked |

### Verdicts
- **confirmed**: a resolvable primary reference supports the secretion claim AND the identity fields
  (uniprot/locus/organism/ss_type) are correct. `provenance_quote` + `provenance_source` required.
- **fix**: the effector is real and secreted, but a field is wrong (dead/off-topic DOI, obsolete
  accession, wrong locus, misclassified subtype). Give the corrected value in a `new_*` column + the
  provenance for the correction. Do NOT drop a real effector just because one field is stale.
- **drop**: the claim does not hold: no primary ref supports secretion, the protein is not a secreted
  substrate/effector of this SS type, it is a duplicate of another row, or it is misclassified into a
  system it does not belong to. `defect` + `notes` required; quote the disconfirming evidence if any.

### `defect` controlled vocabulary
`dead_doi`, `off_topic_ref` (resolves but is about something else), `wrong_uniprot`, `wrong_locus`,
`wrong_genome`, `misclassified_ss` (real effector, wrong SS type), `not_secreted` (not a secreted
substrate/effector), `duplicate`, `weak_evidence` (only homology/prediction, no experimental secretion
assay), `no_primary` (no usable primary reference exists).

## Anti-hallucination rules (non-negotiable)
1. Echo `row_id`/`gene`/`uniprot`/`ss_type` exactly. Never alter them to "correct" a perceived typo.
2. Every `confirmed`/`fix` needs a `provenance_quote` that is a **verbatim** span from a **real**
   source named in `provenance_source`. No source located → cannot be `confirmed`; use
   `secretion_claim_supported=unclear` and explain.
3. Never fabricate a DOI, PMID, PMCID, accession, or author/year. A proposed `new_primary_ref` must be
   a reference you actually retrieved (give its PMID/DOI).
4. Paywalled/unreadable full text is fine: verify from the abstract + UniProt + NCBI; mark
   `doi_resolves=paywalled` and lower `confidence`. Absence of access is not evidence of a defect.
5. Pathogen/toxin effectors are in scope as ordinary biology. Verify them from the literature like any
   other protein; do not refuse or editorialize.

## Per-effector agent prompt (template; one batch of N ledger rows per agent)

> You are auditing a curated answer key of bacterial secreted-protein effectors. For each row below,
> verify that the cited primary reference supports the claim that this protein is a **secreted
> substrate/effector of the stated secretion system**, and that its identifiers are correct. Use the
> PubMed tools (search/metadata/full-text/id-convert) and UniProt; consult NCBI for loci/genomes.
>
> Output ONE TSV row per input row with EXACTLY these columns: row_id, gene, uniprot, ss_type, verdict,
> defect, secretion_claim_supported, doi_resolves, uniprot_ok, provenance_quote, provenance_source,
> new_primary_ref, new_uniprot, new_locus_tag, confidence, notes. Echo row_id/gene/uniprot/ss_type
> verbatim. Follow the anti-hallucination rules: real verbatim quotes with real source ids only; never
> invent a DOI/PMID/accession; if you cannot find support, say `unclear`, do not guess. Return only the
> TSV (header + rows), nothing else.
>
> Rows:
> <tab-separated ledger rows: row_id, ss_type, subtype, gene, uniprot, locus_tag, organism, primary_ref, citation_quote>
