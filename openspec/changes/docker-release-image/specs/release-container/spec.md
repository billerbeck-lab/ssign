## ADDED Requirements

### Requirement: Reproducible pinned image build
The container image SHALL build deterministically: the CUDA base image is referenced by a locked digest (not a floating tag or a placeholder), and the Python dependency set is installed from a committed, fully-resolved lock/constraints file. A rebuild from the same change SHALL resolve the identical base image and package versions.

#### Scenario: No unresolved placeholders remain
- **WHEN** the image is built from `containers/Dockerfile`
- **THEN** `ARG CUDA_BASE` carries a real `@sha256:` digest (no all-zero placeholder) and `pip install` is constrained by the committed lockfile

#### Scenario: Rebuild resolves the same stack
- **WHEN** the image is rebuilt later from the same change
- **THEN** the base image digest and the installed package versions match the recorded ones

### Requirement: Offline pipeline run with bind-mounted data
The image SHALL run the ssign pipeline offline using databases and model weights provided by the host at run time (bind-mounted onto the in-image `.ssign/databases` and `.ssign/models` paths), with no network access required for a run whose tools are all local.

#### Scenario: Golden fixture runs in-container
- **WHEN** the image runs the minimal T5aSS golden fixture with DeepLocPro bind-mounted and all other tools disabled, with no network
- **THEN** the run completes and reports the known T5aSS substrate (BIMENO_04457), matching the golden expected output

### Requirement: One image, per-tier by mounted data
The image SHALL serve every install tier without a separate build per tier: the tier is selected by which host database/weight directories are bind-mounted, since the extended and full tiers share identical Python dependencies. The documentation SHALL state exactly which host directories to mount for base, extended, and full.

#### Scenario: Extended and full use the same image
- **WHEN** a user runs full-tier vs extended-tier analyses
- **THEN** both use the same image tag and differ only in the mounted database volume

### Requirement: Licence-restricted tools are host-provided, not bundled
The image SHALL NOT bundle DTU SignalP/DeepLocPro weights or the EggNOG database. It SHALL support a host-provided local install via bind-mount (`--signalp-path` / `--deeplocpro-path`) and SHALL support the DTU webserver fallback (`--*-mode remote`) for users without a local licence.

#### Scenario: SignalP via bind-mounted local install
- **WHEN** a host SignalP install is bind-mounted and `--signalp-path` points at it
- **THEN** the in-container run uses the local SignalP without contacting the network

#### Scenario: DeepLocPro webserver fallback
- **WHEN** no local DeepLocPro is mounted and `--deeplocpro-mode remote` is set
- **THEN** the run uses the DTU webserver and does not require a bundled weight

### Requirement: HPC (Singularity/Apptainer) conversion path
The image SHALL be convertible to a Singularity/Apptainer `.sif` for HPC hosts that cannot run Docker, and the conversion + run commands SHALL be documented and verified on the target HPC.

#### Scenario: Convert and run as .sif
- **WHEN** the built image is converted to a `.sif` and run on the HPC with databases bind-mounted
- **THEN** the pipeline runs equivalently to the Docker invocation
