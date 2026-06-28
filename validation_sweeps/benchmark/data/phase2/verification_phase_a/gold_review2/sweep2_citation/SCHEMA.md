# Sweep 2 — CITATION review (gold list second pass)

You are auditing the **primary reference** of a curated gold list of bacterial secreted-protein effectors.
Each row pairs one secretion-system substrate with a `primary_ref` (DOI or PMID) and a `citation_quote`
that is supposed to be a verbatim sentence from that reference establishing the protein is **secreted**.

For each row, answer three questions, then assign a verdict:

1. **Organism match** — is the reference about the row's organism? A different *strain* of the same species
   (or a closely related species for a conserved effector family) is acceptable but must be flagged
   (`strain_diff`): e.g. a LEE effector cited from EPEC E2348/69 for an EHEC Sakai row is fine and common.
   A reference about an unrelated organism is `no`.
2. **Documents secretion (STRICT)** — does the reference actually show this protein is a **secreted /
   translocated substrate of this secretion-system type**? Acceptable evidence: detection in culture
   supernatant, translocation into host cells (CyaA/BlaM/fluorescent reporter), secretion-machinery-
   dependent export, or a direct in-vitro secretion assay. NOT acceptable on its own: "is a virulence
   factor", "is required for colonization", "is an adhesin/toxin", sequence homology, or a bare
   bioinformatic prediction. If the only support is homology/prediction, that is `documents_secretion=no`
   with verdict `weak`.
3. **Quote verbatim** — does `citation_quote` appear (verbatim, allowing trivial whitespace/quote-mark
   differences) in the reference's title/abstract/accessible full text?

## How to check (use the tools)

Use the PubMed MCP tools (search_articles, get_article_metadata, get_full_text_article, convert_article_ids,
lookup_article_by_citation) and Crossref/Europe PMC via WebFetch:
- Resolve a DOI: `https://api.crossref.org/works/<doi>` (title, authors, year, container).
- Abstract / full text: Europe PMC `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<doi>&resultType=core&format=json`, or PubMed get_full_text_article for PMC-available papers.
- If only the abstract is accessible, verify from it + UniProt; mark `quote_verbatim=cant_access` for
  full-text-only quotes and lower confidence. Lack of access is NOT a defect.

## Output — ONE TSV row per input row, columns in this exact order

`row_id  gene  uniprot  ss_type  citation_verdict  organism_match  documents_secretion  quote_verbatim  doi_resolves  proposed_primary_ref  proposed_quote  provenance_source  confidence  notes`

| column | values | rule |
|---|---|---|
| `row_id` | echoed | verbatim |
| `gene` | echoed | verbatim |
| `uniprot` | echoed | verbatim |
| `ss_type` | echoed | verbatim |
| `citation_verdict` | `ok` \| `fix_ref` \| `fix_quote` \| `weak` \| `unclear` | see below |
| `organism_match` | `yes` \| `strain_diff` \| `no` \| `unclear` | |
| `documents_secretion` | `yes` \| `no` \| `unclear` | STRICT rubric above |
| `quote_verbatim` | `yes` \| `paraphrased` \| `not_found` \| `cant_access` | |
| `doi_resolves` | `yes` \| `no` \| `paywalled` | the row's `primary_ref` |
| `proposed_primary_ref` | DOI/PMID or blank | only for `fix_ref`; a REAL ref you retrieved |
| `proposed_quote` | verbatim ≤300 chars or blank | only for `fix_quote`/`fix_ref`; a real secretion sentence |
| `provenance_source` | PMID / DOI / PMCID | the source the quote/decision is from |
| `confidence` | `high` \| `medium` \| `low` | |
| `notes` | free text | name every id checked; explain organism nuance |

### Verdicts
- **ok**: ref resolves, organism matches (`yes` or flagged `strain_diff` for a conserved effector),
  documents secretion (strict), and the quote is verbatim.
- **fix_ref**: the cited ref is wrong, off-topic, dead, or does not document secretion, AND you found a
  better real reference → `proposed_primary_ref` + a `proposed_quote` from it.
- **fix_quote**: the ref is correct but `citation_quote` is not verbatim (or is not actually a secretion
  sentence) → give a verbatim `proposed_quote` from the same ref.
- **weak**: the protein's secretion rests only on homology/prediction in this ref; no experimental secretion
  evidence. Flag for human review; do not auto-replace.
- **unclear**: you could not access/resolve enough to judge. Do not guess.

## Anti-hallucination rules (non-negotiable)
1. Echo `row_id`/`gene`/`uniprot`/`ss_type` exactly.
2. Every `proposed_quote` must be a **verbatim** span from a **real** source named in `provenance_source`.
   No source located → cannot be `ok`/`fix_*`; use `unclear`.
3. NEVER fabricate a DOI, PMID, accession, author, or year. A `proposed_primary_ref` must be one you
   actually retrieved (give its DOI/PMID).
4. Paywalled is fine: verify from abstract + UniProt; mark `doi_resolves=paywalled`, lower confidence.
   Absence of access is not evidence of a defect.
5. Pathogen/toxin effectors are ordinary biology; verify from the literature, do not refuse or editorialize.

Return ONLY the TSV (header + one row per input row), nothing else, then write it to the path given in your
task prompt.
