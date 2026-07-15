# annotation-resource-sizing Specification

## Purpose
TBD - created by archiving change annotation-subsystem-cleanup. Update Purpose after archive.
## Requirements
### Requirement: Annotation tools size memory to the parallel-group share
When annotation tools run concurrently in the annotation group, each tool SHALL size its memory (and CPU) budget from its parallel-group share (`parallel_share_ram_gb()` = effective RAM / group size), not from the whole allocation, so the concurrent group fits the job's RAM allocation instead of OOM-killing on small nodes.

#### Scenario: InterProScan heap inside the group
- **WHEN** InterProScan runs inside an N-tool annotation group on a job with R GB effective RAM
- **THEN** its Java max-heap SHALL be sized from the share (R/N), clamped to [4, 64] GB, not from R

#### Scenario: HH-suite workers inside the group
- **WHEN** HH-suite runs inside the annotation group
- **THEN** its worker fan-out SHALL be clamped by the parallel-group RAM share, not by CPU count alone

#### Scenario: pLM-BLAST CPU on a GPU-less node
- **WHEN** pLM-BLAST runs in the annotation group with no GPU present
- **THEN** it SHALL be counted in the annotation CPU budget so it does not oversubscribe cores against its co-scheduled peers

#### Scenario: Standalone tool keeps the full budget
- **WHEN** a tool runs standalone (group size 1)
- **THEN** `parallel_share_ram_gb()` returns the full allocation and the tool's memory sizing is unchanged from today

### Requirement: EggNOG dbmem auto-detects from the RAM share
EggNOG `--dbmem` SHALL default to auto-detection from the parallel-group RAM share (enabled when the share is at least 50 GB, disabled otherwise), and SHALL NOT be force-enabled by a CLI default.

#### Scenario: Small node auto-disables dbmem
- **WHEN** the EggNOG RAM share is below 50 GB and the user passes no dbmem flag
- **THEN** EggNOG SHALL run without `--dbmem` (memory-mapping the on-disk DB), avoiding the ~44 GB in-RAM load and the resulting OOM

#### Scenario: Large node keeps the speedup
- **WHEN** the RAM share is at least 50 GB and the user passes no dbmem flag
- **THEN** EggNOG SHALL use `--dbmem`

#### Scenario: Explicit override wins
- **WHEN** the user passes `--eggnog-dbmem` or `--no-eggnog-dbmem`
- **THEN** that explicit choice SHALL override the autodetect

