"""Cgroup/scheduler-aware CPU and RAM detection.

`os.cpu_count()` and `psutil.virtual_memory().total` read host-wide values
and ignore the job's actual allocation under PBS/SLURM cgroups. On a shared
HPC node this leads to massive thread oversubscription (24 CPUs requested
on a 4-CPU allocation thrashes at ~2% efficiency) and OOM-by-default for
memory-heavy steps. These helpers query the real allocation.

Lookup priority for both is: explicit env override → scheduler env vars →
scheduler CLI (qstat/scontrol) → cgroup v2 → cgroup v1 → host total.
Take the smallest reasonable value, since under-using cores is far cheaper
than thrashing them.
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
from typing import cast

logger = logging.getLogger(__name__)

# cgroup v1 "no limit" sentinel: math.inf, encoded as 2**63 - 1 rounded
# down to the kernel page size. Anything within ~1 page of this is the
# "unlimited" marker, not a real limit.
_CGROUP_UNLIMITED_BYTES = 2**63 - 4096


# GPU-predictor batch sizes per VRAM tier. The breakpoints come from
# measuring peak activation memory for ProtT5 at max_length=512 plus
# tokenizer/model weight overhead, leaving ~20% headroom on top to survive
# Bakta and EggNOG running alongside. ProtT5 is larger than DeepSecE's
# ESM-1b, so these tiers are conservative for the current caller. Used by
# auto_batch_size_from_vram below; callers can pass --batch-size N to override.
_AUTO_BATCH_TIERS = (
    # (min_vram_gib_inclusive, batch_size). Thresholds sit 2-4 GiB below
    # the marketing-GB capacity because torch reports total - firmware
    # reserve in GiB (1024^3): A40's "48 GB" comes back as ~44.4 GiB,
    # RTX 4090's "24 GB" as ~22 GiB, RTX 3090's "12 GB" as ~11 GiB.
    # Top tier sized from measured A40 peak: 15 GB at batch=16 (12 GB ProtT5
    # weights + 3 GB activations). Extrapolating linearly, batch=64 lands
    # near 12 + 12 = 24 GB, half of A40's 44 GiB usable.
    (40, 64),  # A40 / L40S / A100 (40/80 GB) / H100 (80 GB)
    (20, 16),  # RTX 4090 / A5000 (24 GB nominal, ~22 GiB reported)
    (10, 8),  # RTX 3090 / RTX 4080 (12-16 GB)
    (0, 4),  # smaller GPUs (8 GB and under) and CPU fallback
)


def probe_cuda_device() -> tuple[str | None, float | None]:
    """Return (device name, total VRAM GiB) for CUDA device 0, or (None, None).

    Lazy-imports torch so callers in tiers that don't ship torch (e.g.
    base) don't pay the import cost. Swallows the usual failure modes
    (no torch, no CUDA, driver mismatch) by returning the sentinel pair.
    """
    try:
        import torch
    except ImportError:
        return (None, None)
    if not torch.cuda.is_available():
        return (None, None)
    try:
        name = torch.cuda.get_device_name(0)
        total_gib = torch.cuda.get_device_properties(0).total_memory / (1024**3)
    except Exception as e:
        logger.warning("Could not read GPU properties: %s", e)
        return (None, None)
    return (name, total_gib)


def auto_batch_size_from_vram(default_when_no_gpu: int = 4) -> int:
    """Pick a GPU-predictor batch size from the active CUDA device's total VRAM.

    Used by the ESM-based predictors (currently DeepSecE). Returns the
    smallest tier table entry's batch size when no CUDA device is visible
    or torch isn't importable; callers can override via ``--batch-size N``.
    """
    _name, total_gib = probe_cuda_device()
    if total_gib is None:
        return default_when_no_gpu
    for min_gb, batch in _AUTO_BATCH_TIERS:
        if total_gib >= min_gb:
            logger.info(
                "GPU auto-batch: detected %.1f GB VRAM, choosing batch_size=%d",
                total_gib,
                batch,
            )
            return batch
    return default_when_no_gpu


def effective_cpu_count() -> int:
    """Number of CPUs the current process is allowed to schedule on.

    Honors cgroup CPU pinning (PBS cpuset, SLURM `--cpus-per-task`). Falls
    back to host total when `sched_getaffinity` isn't available (non-Linux).
    """
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return os.cpu_count() or 4


def resolve_threads(n: int | None) -> int:
    """Pick a thread count: ``n`` when provided, else the parallel-group share.

    Standard shape for tool wrappers' ``threads`` argument so they all
    auto-scale to the scheduler's allocation when called with None.
    Routes through ``parallel_share_cpus()`` so wrappers launched inside
    a runner ``parallel_group(N)`` automatically divide the budget by N
    instead of each grabbing the whole allocation (which causes 3-4x CPU
    oversubscription on shared HPC nodes — observed live as cput=4h /
    walltime=36m on the K-12 RTX6000 annotation stage). When the env
    var is unset (standalone wrapper invocations, no parallel group),
    the share helper returns the full ``effective_cpu_count()``.
    Always returns at least 1 even if cgroup math degenerates.
    """
    if n is not None:
        return max(1, n)
    return max(1, parallel_share_cpus())


# Set by PipelineRunner before launching the DLP/DSE/SignalP parallel
# group so each wrapper subprocess can self-divide the CPU budget.
# Unset means "running standalone" — the wrapper uses the full
# effective_cpu_count.
_PARALLEL_GROUP_SIZE_ENV = "SSIGN_PARALLEL_GROUP_SIZE"


def _parallel_group_size() -> int:
    """The SSIGN_PARALLEL_GROUP_SIZE divisor: N when set to an int >= 2, else 1.

    Set by the runner's `parallel_group(N)` before it launches the prediction /
    annotation tools concurrently. 1 (env unset, non-int, or <= 1) means
    "running standalone" — the caller keeps its full CPU / RAM budget. Shared by
    parallel_share_cpus / parallel_share_ram_gb so the two can't drift.
    """
    raw = os.environ.get(_PARALLEL_GROUP_SIZE_ENV)
    if not raw:
        return 1
    try:
        n = int(raw)
    except ValueError:
        return 1
    return n if n > 1 else 1


def parallel_share_cpus() -> int:
    """Return this subprocess's CPU share when it's inside a parallel group.

    `effective_cpu_count() // N` (N = the parallel-group size), clamped to >= 1.
    Standalone (N = 1) keeps the full count so lone wrapper invocations aren't
    penalised. Stops tool wrappers 3-4x oversubscribing a shared HPC node.
    """
    return max(1, effective_cpu_count() // _parallel_group_size())


def parallel_share_ram_gb() -> float:
    """Return this subprocess's RAM share (GB) when it's inside a parallel group.

    The RAM analogue of `parallel_share_cpus()`: `effective_ram_gb() / N`. A
    wrapper that sizes memory to "fill the node" (EggNOG's `--dbmem` + DIAMOND
    `--block_size`) would otherwise overcommit when it runs concurrently with the
    other annotation tools — each grabbing the full allocation OOM-kills the
    biggest allocator (observed live: EggNOG DIAMOND OOM in a 6-tool full-tier
    group on a 120 GB node). Standalone (N = 1) returns the full RAM so
    single-tool jobs keep the speedup.
    """
    return effective_ram_gb() / _parallel_group_size()


from contextlib import contextmanager  # noqa: E402


@contextmanager
def parallel_group(n: int):
    """Mark this process region as a parallel group of size `n`.

    Sets SSIGN_PARALLEL_GROUP_SIZE for the duration of the `with` block
    so subprocess children (the prediction tool wrappers) read the
    right value via `parallel_share_cpus()`, then restores the prior
    value on exit even if the body raises. Replaces the 6-line manual
    save/set/pop dance in the runner's parallel-group launch.
    """
    prev = os.environ.get(_PARALLEL_GROUP_SIZE_ENV)
    os.environ[_PARALLEL_GROUP_SIZE_ENV] = str(n)
    try:
        yield
    finally:
        if prev is None:
            os.environ.pop(_PARALLEL_GROUP_SIZE_ENV, None)
        else:
            os.environ[_PARALLEL_GROUP_SIZE_ENV] = prev


def humanise_bytes(n: int) -> str:
    """Pretty-print bytes as B / KB / MB / GB / TB.

    Uses 1024-base (`KiB` math, `KB` label — matches `du -h` and
    `ls -lh` convention). Caps at TB; PB+ values render as e.g.
    ``"1024.0 TB"`` rather than introducing a new unit row.
    """
    f: float = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if f < 1024:
            return f"{int(f)} B" if unit == "B" else f"{f:.1f} {unit}"
        f /= 1024
    return f"{f:.1f} TB"


def _parse_size_to_gb(value: str) -> float | None:
    """Parse a size string like '32gb', '128GiB', '4096mb' → GB (decimal)."""
    m = re.match(r"\s*([\d.]+)\s*([kmgt]i?b?)?\s*$", value.lower())
    if not m:
        return None
    n = float(m.group(1))
    unit = (m.group(2) or "").rstrip("b")
    mult = {
        "": 1 / 2**30,
        "k": 1 / 2**20,
        "ki": 1 / 2**20,
        "m": 1 / 1024,
        "mi": 1 / 1024,
        "g": 1.0,
        "gi": 1.0,
        "t": 1024.0,
        "ti": 1024.0,
    }.get(unit)
    return n * mult if mult else None


def _pbs_alloc_gb() -> float | None:
    """Parse Resource_List.mem from `qstat -f $PBS_JOBID`."""
    job = os.environ.get("PBS_JOBID")
    if not job:
        return None
    try:
        out = subprocess.run(["qstat", "-f", job], capture_output=True, text=True, timeout=5)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:  # job vanished, auth error, etc.
        return None
    m = re.search(r"Resource_List\.mem\s*=\s*(\S+)", out.stdout)
    return _parse_size_to_gb(m.group(1)) if m else None


def _cgroup_mem_gb() -> float | None:
    """Read the memory cgroup limit, skipping the 'unlimited' sentinels."""
    # cgroup v2 reports the literal string "max" when no limit is set.
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(path) as fh:
                raw = fh.read().strip()
        except OSError:
            continue
        if raw == "max" or not raw.isdigit():
            continue
        n = int(raw)
        if n >= _CGROUP_UNLIMITED_BYTES:
            continue
        return n / 2**30
    return None


def host_ram_gb() -> float:
    """Total physical RAM on the host (psutil), or 0.0 if unreadable.

    Distinct from `effective_ram_gb()`, which returns the smaller of the
    scheduler allocation and the host total. doctor uses both so users
    can see the gap (the most common "scheduler is silently throttling
    me" signal on shared HPC nodes).
    """
    try:
        import psutil

        return cast(float, psutil.virtual_memory().total / 2**30)
    except Exception:
        return 0.0


def effective_ram_gb() -> float:
    """Smallest of every job-allocation signal we can find, falling back to host.

    Sources tried (returns the minimum of all that produce a positive value):
    `SSIGN_MAX_RAM_GB` env → `SLURM_MEM_PER_NODE` env → PBS qstat parse →
    cgroup v2/v1 limit → psutil host total. The min handles HPCs that
    enforce RAM externally (PBS) while the kernel cgroup still reports
    "unlimited" — the qstat parse catches what /sys misses.
    """
    candidates: list[float] = []

    override = os.environ.get("SSIGN_MAX_RAM_GB")
    if override:
        try:
            candidates.append(float(override))
        except ValueError:
            logger.warning(f"SSIGN_MAX_RAM_GB={override!r} not a number; ignored")

    slurm_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mb and slurm_mb.isdigit():
        candidates.append(int(slurm_mb) / 1024)

    pbs = _pbs_alloc_gb()
    if pbs is not None:
        candidates.append(pbs)

    cg = _cgroup_mem_gb()
    if cg is not None:
        candidates.append(cg)

    host = host_ram_gb()
    if host > 0:
        candidates.append(host)

    return min(c for c in candidates if c > 0) if candidates else 0.0


# Filesystem types whose per-read latency is high enough that database
# random-access mmap stalls for tens of minutes on multi-GB SQLite/BLAST
# files. Caching to local SSD turns the access into a local page fault.
# Conservative list: only filesystems that are *always* networked. Local
# variants of these (e.g. localhost-NFS) are rare and the false positive
# is a redundant copy, not data loss.
_REMOTE_FS_TYPES = frozenset({"nfs", "nfs4", "gpfs", "lustre", "cifs", "smb", "smb3", "sshfs", "fuse.sshfs"})


def _fs_type(path: str) -> str | None:
    """Filesystem type for `path` as recorded in /proc/mounts. None if unknown.

    Walks parent directories until a mount-point match is found, picking the
    longest matching prefix (handles nested mounts).
    """
    try:
        with open("/proc/mounts") as fh:
            mounts = [line.split() for line in fh if line.strip()]
    except OSError:
        return None
    target = os.path.realpath(path)
    best: tuple[str, str | None] = ("", None)
    for parts in mounts:
        if len(parts) < 3:
            continue
        mp, fstype = parts[1], parts[2]
        if (target == mp or target.startswith(mp.rstrip("/") + "/")) and len(mp) > len(best[0]):
            best = (mp, fstype)
    return best[1]


def is_remote_filesystem(path: str) -> bool:
    """True if `path` lives on a networked filesystem (gpfs/nfs/lustre/...).

    Used by wrappers that want to copy DB files to local SSD only when the
    source is genuinely networked, so PC users don't pay a needless copy
    cost when both source and TMPDIR are on the same local SSD.
    """
    fs = _fs_type(path)
    return fs in _REMOTE_FS_TYPES if fs else False


def stage_to_local_ssd_if_remote(
    src_dir: str,
    cache_dir: str,
    required: tuple[str, ...],
    optional: tuple[str, ...] = (),
    min_free_gb: float = 60.0,
) -> str:
    """Copy `required + optional` files from `src_dir` to `cache_dir/<basename>/` if src is on a network FS.

    Returns the path the caller should pass to the tool as `--data_dir` (or
    equivalent): the staged local dir when caching, otherwise `src_dir`
    unchanged. Idempotent — re-runs in the same job see size-equal files
    and skip the copy.

    Why: random-access mmap on multi-GB SQLite / BLAST databases stalls
    for tens of minutes on gpfs/nfs/lustre. Copying to node-local SSD
    converts subsequent lookups to local page faults. On a regular PC
    where both src and cache are on the same local disk, this is a no-op.

    Raises FileNotFoundError if any `required` file is missing from
    `src_dir` (cheap pre-flight to avoid a half-copied broken cache).
    """
    if not is_remote_filesystem(src_dir):
        logger.info(f"DB at {src_dir} is on local filesystem; skipping cache copy")
        return src_dir

    for name in required:
        if not os.path.exists(os.path.join(src_dir, name)):
            raise FileNotFoundError(f"Required DB file missing in {src_dir}: {name}")

    import shutil

    free_gb = shutil.disk_usage(cache_dir).free / 2**30
    if free_gb < min_free_gb:
        logger.warning(
            f"Skipping local DB cache: {cache_dir} has {free_gb:.1f} GB free "
            f"(need {min_free_gb}). Reads will go over the network — expect "
            f"latency penalties."
        )
        return src_dir

    target_dir = os.path.join(cache_dir, os.path.basename(src_dir.rstrip("/")) or "db")
    os.makedirs(target_dir, exist_ok=True)
    for name in required + optional:
        src = os.path.join(src_dir, name)
        if not os.path.exists(src):
            continue
        dst = os.path.join(target_dir, name)
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src):
            logger.info(f"DB file already cached: {dst}")
            continue
        size_gb = os.path.getsize(src) / 2**30
        logger.info(f"Caching {name} ({size_gb:.1f} GB) -> {target_dir}")
        tmp = f"{dst}.tmp.{os.getpid()}"
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)  # atomic; safe against concurrent stagers
    return target_dir


def stage_prefix_files_to_local_ssd_if_remote(
    prefix: str,
    cache_dir: str,
    min_free_gb: float = 10.0,
) -> str:
    """Stage files matching `glob(prefix*)` to local SSD when source is remote.

    Designed for ffindex / `<prefix>_<suffix>.ext` DBs (HH-suite uniclust,
    pfam, pdb70) where the "DB" is a prefix that expands at runtime to
    several sibling files in the same directory. Avoids the three traps
    of staging the whole parent dir:
      1. Sibling DBs sharing the same parent get copied unnecessarily.
      2. Prefix at filesystem root (`os.path.dirname` → `/rds`) would
         rsync an entire mount.
      3. Multiple per-prefix calls each apply min_free_gb independently,
         causing partial staging when the cache barely fits.

    Returns the new prefix (cache_dir/<basename>/<prefix-basename>) when
    cached, or the original prefix unchanged when src is local / not
    cacheable. Caller treats the return value as opaque.
    """
    import glob

    if not prefix or not is_remote_filesystem(prefix):
        return prefix

    matches = sorted(glob.glob(f"{prefix}*"))
    if not matches:
        logger.warning(f"Prefix-glob found no files at {prefix}*; skipping cache")
        return prefix

    import shutil

    needed_gb = sum(os.path.getsize(p) for p in matches) / 2**30
    free_gb = shutil.disk_usage(cache_dir).free / 2**30
    if free_gb < max(min_free_gb, needed_gb * 1.1):
        logger.warning(
            f"Skipping prefix-cache for {prefix}: need ~{needed_gb:.1f} GB, have {free_gb:.1f} GB free in {cache_dir}."
        )
        return prefix

    # cache_dir/db_prefixes/<basename>/ keeps multiple prefixes from
    # colliding (each gets its own subdir keyed by the prefix basename).
    target_dir = os.path.join(cache_dir, "db_prefixes", os.path.basename(prefix))
    os.makedirs(target_dir, exist_ok=True)
    logger.info(f"Staging {len(matches)} files for prefix {prefix} ({needed_gb:.1f} GB) -> {target_dir}")
    for src in matches:
        dst = os.path.join(target_dir, os.path.basename(src))
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src):
            continue
        tmp = f"{dst}.tmp.{os.getpid()}"
        shutil.copy2(src, tmp)
        os.replace(tmp, dst)
    return os.path.join(target_dir, os.path.basename(prefix))


def stage_directory_tree_to_local_ssd_if_remote(
    src_dir: str,
    cache_dir: str,
    min_free_gb: float = 60.0,
) -> str:
    """Like stage_to_local_ssd_if_remote but for whole directory trees.

    Bakta's DB has ~80 GB across many nested subdirs (amrfinderplus-db/,
    ncRNA/, rRNA/, ...); listing every file explicitly is brittle. This
    helper rsync's the entire tree to local SSD when src is on a network
    FS and skips the copy when it's local. Re-runs in the same job see
    the existing tree and short-circuit (rsync's --size-only).

    Requires `rsync` on PATH (universally present on Linux HPC).
    """
    if not is_remote_filesystem(src_dir):
        logger.info(f"DB tree at {src_dir} is on local filesystem; skipping cache copy")
        return src_dir

    import shutil
    import subprocess

    if shutil.which("rsync") is None:
        logger.warning("rsync not on PATH; cannot stage DB tree to local SSD")
        return src_dir

    free_gb = shutil.disk_usage(cache_dir).free / 2**30
    if free_gb < min_free_gb:
        logger.warning(
            f"Skipping local DB tree cache: {cache_dir} has {free_gb:.1f} GB free "
            f"(need {min_free_gb}). Reads will go over the network."
        )
        return src_dir

    target_dir = os.path.join(cache_dir, os.path.basename(src_dir.rstrip("/")) or "db_tree")
    os.makedirs(target_dir, exist_ok=True)
    # `-a` preserves perms/times so the next run's --size-only check
    # short-circuits. Trailing slash on src means "copy contents into".
    logger.info(f"Staging DB tree from {src_dir} -> {target_dir} (may take ~1-3 min)")
    try:
        subprocess.run(
            ["rsync", "-a", "--size-only", f"{src_dir.rstrip('/')}/", f"{target_dir}/"],
            check=True,
            timeout=3600,  # 1h ceiling: 84 GB at degraded GPFS speed (30 MB/s)
            # would take ~47 min, so 30 min is too tight; 1h is comfortable.
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logger.warning(f"rsync of {src_dir} failed: {e}; falling back to remote path")
        return src_dir
    return target_dir
