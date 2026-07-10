## ADDED Requirements

### Requirement: Version-locked default artifact (Singularity/Apptainer)
The default reproducible install SHALL be a Singularity/Apptainer image built from a digest-pinned base and the committed, fully-resolved `uv.lock`. A rebuild from the same change SHALL resolve the identical base image and Python package versions. The image SHALL run with no root.

#### Scenario: Rebuild resolves the same stack
- **WHEN** the `.sif` is rebuilt from the same change
- **THEN** the base image digest and the installed package versions match the recorded ones

#### Scenario: Runs unprivileged
- **WHEN** a user runs the `.sif` on a host where they have no root (e.g. an HPC login node)
- **THEN** the pipeline runs without requiring elevated privileges

### Requirement: Offline pipeline run with bind-mounted data
The image SHALL run the ssign pipeline offline using databases and model weights provided by the host at run time (bind-mounted onto the in-image `.ssign/databases` and `.ssign/models` paths), with no network access required for a run whose tools are all local.

#### Scenario: Golden fixture runs in-container
- **WHEN** the image runs the minimal T5aSS golden fixture with DeepLocPro bind-mounted, all other tools disabled, and no network
- **THEN** the run completes and reports the known T5aSS substrate (BIMENO_04457), matching the golden expected output

### Requirement: Base-only conda environment for macOS
A conda `environment.yml` SHALL install ssign at the **base** tier on macOS, pinned to the same Python versions as the lock. It SHALL NOT include any Linux-only tool (InterProScan, HH-suite, BLAST+, Bakta); base's tools (DeepLocPro, SignalP, DeepSecE, ProtParam, MacSyFinder) are Mac-compatible.

#### Scenario: Mac base install runs the pipeline
- **WHEN** a macOS user creates the env from `environment.yml` and runs the golden fixture at base tier
- **THEN** the base pipeline runs and no extended/full (Linux-only) tool is invoked

#### Scenario: No Linux-only tool in the Mac env
- **WHEN** the conda environment is resolved on macOS
- **THEN** it contains no InterProScan / HH-suite / BLAST+ / Bakta dependency

### Requirement: Licence-restricted tools are host-provided, not bundled
Neither artifact SHALL bundle DTU SignalP/DeepLocPro weights or the EggNOG database. The container SHALL support a host-provided local install via bind-mount (`--signalp-path` / `--deeplocpro-path`) and the DTU webserver fallback (`--*-mode remote`).

#### Scenario: SignalP via bind-mounted local install
- **WHEN** a host SignalP install is bind-mounted and `--signalp-path` points at it
- **THEN** the in-container run uses the local SignalP without contacting the network

### Requirement: Windows is not a target
Native Windows SHALL NOT be a supported platform (several core tools have no Windows build). The documentation SHALL direct Windows users to WSL2, which is served by the Linux path.

#### Scenario: Windows user directed to WSL2
- **WHEN** a Windows user consults the install docs
- **THEN** they are directed to WSL2 (the Linux/Singularity or conda path), not a native-Windows install
