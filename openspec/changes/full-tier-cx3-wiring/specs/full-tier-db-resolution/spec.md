## ADDED Requirements

### Requirement: Automatic BLAST NR database resolution

The pipeline SHALL resolve the BLAST NR database for BLASTp automatically when the operator has not set `blastp_db` explicitly, using the same environment-variable and marker-root resolution path as every other optional database. Because BLAST+ `-db` takes a database **name prefix** while the resolver locates the containing **directory** (via the `nr.pdb` sentinel), the pipeline SHALL append the `nr` basename to the resolved directory to form the `-db` value. An explicit `--blastp-db` value SHALL always take precedence over automatic resolution.

#### Scenario: NR resolved from environment variable

- **WHEN** `blastp_db` is unset and `SSIGN_BLAST_NR` points at a directory containing `nr.pdb`
- **THEN** the runner sets `blastp_db` to `<that directory>/nr` and BLASTp runs without an explicit `--blastp-db` flag

#### Scenario: Explicit path overrides auto-resolution

- **WHEN** the operator passes `--blastp-db /some/other/db` and `SSIGN_BLAST_NR` is also set
- **THEN** the runner uses `/some/other/db` unchanged and does not consult `SSIGN_BLAST_NR`

#### Scenario: No NR available still errors clearly

- **WHEN** BLASTp is enabled, `blastp_db` is unset, and NR is not discoverable via `SSIGN_BLAST_NR` or the marker root
- **THEN** the BLASTp step returns the existing definitive "BLASTp requires a local database" error rather than starting and crashing

### Requirement: Full-tier CX3 batch job wires its databases

The CX3 batched submit path SHALL support running at the `full` tier. `submit_batched_overnight.sh` SHALL accept a `--tier <base|extended|full>` flag (default `extended`) and pass the selected tier through to the run. The PBS wrapper SHALL export `SSIGN_HHSUITE_UNICLUST` (UniRef30) and `SSIGN_BLAST_NR` (NR) from the databases root when those directories are present, so that HH-suite's MSA step and BLASTp find their databases with no per-run database flags.

#### Scenario: Full-tier submit runs BLASTp and HH-suite

- **WHEN** the operator runs `submit_batched_overnight.sh --tier full ...` on a node where `blast_nr/` and `hhsuite/uniref30/` exist under the databases root
- **THEN** the job runs `ssign run --tier full`, BLASTp and HH-suite are enabled, and both resolve their databases from the exported environment variables

#### Scenario: Default tier unchanged

- **WHEN** the operator runs `submit_batched_overnight.sh` with no `--tier` flag
- **THEN** the job runs at the `extended` tier exactly as before this change
