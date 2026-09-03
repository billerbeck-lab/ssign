# Figure reference lists

Evidence bases for the secretion-system overview figure used in the ssign
manuscript and in Teo Reid's MRes thesis. They are kept here so the figure's
claims can be checked source by source without the reference list bloating the
manuscript bibliography.

| File | Covers | Entries |
| --- | --- | --- |
| [`REFERENCES_MACHINES.md`](REFERENCES_MACHINES.md) | The machine figure: one entry per drawn component, with the PDB accession, method and resolution behind it, plus the in-situ maps and copy-number evidence relied on but not drawn. | 35 DOIs |
| [`REFERENCES.md`](REFERENCES.md) | The substrate figure: one section per secretion signal, tied to the specific claim it supports. | 52 DOIs |

No reference appears in both files.

Every structural citation in `REFERENCES_MACHINES.md` was retrieved from the RCSB
PDB entry record for the exact accession the figure draws, so authors and DOIs
cannot drift from the deposition.

## Related

Expected subcellular localizations of secretion-system machinery components, with
a primary citation per component, are a separate list:
[`src/ssign_app/scripts/ssign_lib/data/component_localizations.tsv`](../../src/ssign_app/scripts/ssign_lib/data/component_localizations.tsv)
(92 components, 53 unique primary references). The experimentally supported
system/substrate pairs used for benchmarking carry their own `primary_ref` and
`citation_quote` columns in
[`docs/benchmark/ssign_benchmarking_list.csv`](../benchmark/ssign_benchmarking_list.csv).
