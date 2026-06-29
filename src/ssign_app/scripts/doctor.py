"""``ssign doctor`` — verify a fresh ssign install on the user's machine.

Three classes of dependency, three sections of report:

1. Python packages + lazy symbols (the same list the integration test asserts).
   These should always pass after ``pip install ssign[extended]``; if any fail
   here, ``pyproject.toml`` is incomplete and CI should catch it before release.
2. External binaries that pip cannot install (HH-suite, BLAST+, InterProScan,
   EggNOG-mapper, optional DTU SignalP/DeepLocPro).
3. Databases + model weights on disk (gigabytes, fetched separately).

Each failure prints the exact fix command. Exit non-zero on any required
failure so ``ssign doctor && ssign run …`` works in scripts.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import sys
from dataclasses import dataclass

from ssign_app.scripts.ssign_lib.dependency_manifest import (
    DatabasePath,
    ExternalBinary,
    ModelWeights,
    PythonDep,
    Tier,
    binaries_for_tier,
    databases_for_tier,
    deps_for_tier,
    find_db_by_env_var,
    weights_for_tier,
)

DEFAULT_DATA_ROOT = os.path.expanduser("~/.ssign")
DEFAULT_TIER: Tier = "extended"


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    fix: str = ""


# ---------------------------------------------------------------------------
# Section 1 — Python imports
# ---------------------------------------------------------------------------


def _resolve_symbol(dotted: str) -> None:
    """Import ``a.b.c.D`` by importing ``a.b.c`` then resolving attribute ``D``.

    Raises ``ImportError`` / ``AttributeError`` on failure (the caller catches).
    """
    parts = dotted.split(".")
    # walk back until we hit something importable; the rest are attributes
    for split in range(len(parts), 0, -1):
        mod_name = ".".join(parts[:split])
        try:
            mod = importlib.import_module(mod_name)
        except ImportError:
            continue
        obj = mod
        for attr in parts[split:]:
            obj = getattr(obj, attr)
        return
    raise ImportError(f"could not import any prefix of {dotted!r}")


def check_python_dep(dep: PythonDep) -> CheckResult:
    try:
        importlib.import_module(dep.module)
    except Exception as e:  # noqa: BLE001 — surface anything from native ext loaders
        return CheckResult(
            name=dep.module,
            ok=False,
            detail=f"{type(e).__name__}: {e}",
            fix=f"pip install {dep.pip_name}",
        )
    for sym in dep.symbols:
        try:
            _resolve_symbol(sym)
        except Exception as e:  # noqa: BLE001
            return CheckResult(
                name=f"{dep.module} (symbol {sym})",
                ok=False,
                detail=f"{type(e).__name__}: {e}",
                fix=f"pip install --upgrade {dep.pip_name}  # symbol moved/missing in installed version",
            )
    return CheckResult(name=dep.module, ok=True, detail="")


# ---------------------------------------------------------------------------
# Section 2 — External binaries
# ---------------------------------------------------------------------------


def check_external_binary(b: ExternalBinary, db_root: str = "") -> CheckResult:
    found = shutil.which(b.binary)
    if found:
        return CheckResult(name=b.name, ok=True, detail=found)
    # Some tools (e.g. InterProScan) ship as a tarball the user extracts
    # somewhere, then points an env var at the install dir. Honour that
    # so doctor doesn't false-flag the binary as missing.
    if b.install_dir_env:
        install_dir = os.environ.get(b.install_dir_env, "").strip()
        if install_dir:
            candidate = os.path.join(install_dir, b.binary)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return CheckResult(name=b.name, ok=True, detail=candidate)
    # Final fallback: mirror what the runner does at execution time.
    # If a DatabasePath shares this binary's install_dir_env (e.g.
    # SSIGN_INTERPROSCAN_PATH on both the InterProScan ExternalBinary
    # and the InterProScan DB DatabasePath), use the DB's resolve_path
    # to find the install via the fetch_databases.sh layout.
    if b.install_dir_env and db_root:
        dbp = find_db_by_env_var(b.install_dir_env)
        if dbp is not None:
            install_dir = dbp.resolve_path(db_root) or ""
            if install_dir:
                candidate = os.path.join(install_dir, b.binary)
                if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                    return CheckResult(name=b.name, ok=True, detail=candidate)
    detail = f"`{b.binary}` not on PATH"
    if b.install_dir_env:
        detail += f" (and ${b.install_dir_env} unset or doesn't contain it)"
    return CheckResult(
        name=b.name + (" [optional]" if b.optional else ""),
        ok=b.optional,
        detail=detail,
        fix=b.install_hint,
    )


# ---------------------------------------------------------------------------
# Section 3 — Databases + model weights
# ---------------------------------------------------------------------------


def read_db_root_marker(data_root: str) -> str | None:
    """Return the ``$TARGET`` value ``fetch_databases.sh`` recorded, or None.

    The install script writes its absolute target path to
    ``<data_root>/db_root`` at the end of a successful run, so users don't
    have to set ``SSIGN_*`` env vars for every DB just for doctor to find
    them. Returns None when the marker doesn't exist or points at a path
    that's since gone away.
    """
    marker = os.path.join(data_root, "db_root")
    if not os.path.isfile(marker):
        return None
    try:
        with open(marker, encoding="utf-8") as f:
            recorded = f.read().strip()
    except OSError:
        return None
    return recorded if recorded and os.path.isdir(recorded) else None


def resolve_db_root(data_root: str) -> str:
    """The directory where bakta/, eggnog/, etc. live as direct children.

    Precedence: ``<data_root>/db_root`` marker → ``<data_root>/databases``.
    """
    return read_db_root_marker(data_root) or os.path.join(data_root, "databases")


def check_database(d: DatabasePath, db_root: str) -> CheckResult:
    # Single source of truth: DatabasePath.resolve_path — same resolver
    # the runner consumes in its __post_init__. Doctor adds the diagnostic
    # "dir exists but sentinel missing" branch on top to help users
    # distinguish "I haven't fetched it yet" from "I fetched it but it's
    # broken / mis-layouted".
    resolved = d.resolve_path(db_root)
    if resolved:
        return CheckResult(name=d.name, ok=True, detail=resolved)
    candidate_dir = os.path.join(db_root, d.default_subpath)
    if os.path.isdir(candidate_dir):
        return CheckResult(
            name=d.name,
            ok=False,
            detail=f"{candidate_dir} exists but no {d.sentinel_file!r} matches inside",
            fix=d.install_hint,
        )
    return CheckResult(
        name=d.name,
        ok=False,
        detail=f"not found at {candidate_dir} (env: ${d.env_var})",
        fix=d.install_hint,
    )


def check_weights(w: ModelWeights, data_root: str, db_root: str) -> CheckResult:
    base = db_root if w.under_db_root else data_root
    resolved = os.path.join(base, w.default_subpath)
    if os.path.exists(resolved):
        return CheckResult(name=w.name, ok=True, detail=resolved)
    return CheckResult(
        name=w.name,
        ok=False,
        detail=f"not found at {resolved}",
        fix=w.install_hint,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def _format_gpu(name: str | None, total_gib: float | None) -> str:
    """Render `probe_cuda_device()`'s return as a one-line label."""
    if name is None or total_gib is None:
        return "none (no CUDA device or torch not installed)"
    return f"{name}, {total_gib:.1f} GiB VRAM"


def report_resources(stream) -> None:
    """Print detected CPU / RAM / GPU + the auto-scaled values ssign will use.

    Informational only — never contributes to `failures`. The aim is to
    let HPC users confirm at a glance that ssign sees the full PBS/SLURM
    allocation (one of the most common silent-throttle bugs).
    """
    # Lazy import: lets the monkeypatch.setattr in tests resolve to the
    # patched names per-call, and keeps psutil/torch out of doctor's
    # module-load path so `ssign doctor --imports-only` stays cheap.
    from ssign_app.scripts.ssign_lib.resources import (
        auto_batch_size_from_vram,
        effective_cpu_count,
        effective_ram_gb,
        host_ram_gb,
        probe_cuda_device,
    )

    cpu = effective_cpu_count()
    host_cpu = os.cpu_count() or 0
    ram = effective_ram_gb()
    host_ram = host_ram_gb()
    gpu_label = _format_gpu(*probe_cuda_device())
    plme_batch = auto_batch_size_from_vram()

    section = "Resources"
    print(f"\n{section}  (informational)", file=stream)
    print("─" * len(section), file=stream)
    print(f"  CPU      effective: {cpu}  (host: {host_cpu})", file=stream)
    if cpu < host_cpu:
        print(f"           note: scheduler is restricting ssign to {cpu}/{host_cpu} cores", file=stream)
    print(f"  RAM      effective: {ram:.1f} GB  (host: {host_ram:.1f} GB)", file=stream)
    print(f"  GPU      {gpu_label}", file=stream)
    print(f"  PLM-E    auto batch size: {plme_batch}  (override with --plme-batch-size N)", file=stream)


def _render(section: str, results: list[CheckResult], stream) -> tuple[int, int]:
    ok_count = sum(1 for r in results if r.ok)
    total = len(results)
    print(f"\n{section}  ({ok_count}/{total} OK)", file=stream)
    print("─" * len(section), file=stream)
    for r in results:
        marker = "OK  " if r.ok else "FAIL"
        print(f"  {marker}  {r.name}", file=stream)
        if r.detail:
            print(f"           {r.detail}", file=stream)
        if not r.ok and r.fix:
            print(f"           fix: {r.fix}", file=stream)
    return ok_count, total


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run(
    tier: Tier = DEFAULT_TIER,
    imports_only: bool = False,
    data_root: str = DEFAULT_DATA_ROOT,
    stream=sys.stdout,
) -> int:
    """Run all checks. Returns 0 if every required check passed, 1 otherwise.

    ``imports_only=True`` skips binaries / DBs / weights — used by CI where
    only the Python environment can be verified.
    """
    db_root = resolve_db_root(data_root)
    marker_used = read_db_root_marker(data_root) is not None

    print(f"ssign doctor — checking tier '{tier}'", file=stream)
    print(f"  data root: {data_root}", file=stream)
    if marker_used:
        print(f"  db root:   {db_root}  (from {data_root}/db_root, written by fetch_databases.sh)", file=stream)
    else:
        print(
            f"  db root:   {db_root}  (default; run fetch_databases.sh to record a custom location)",
            file=stream,
        )

    failures = 0

    py_results = [check_python_dep(d) for d in deps_for_tier(tier)]
    ok, total = _render("Python packages", py_results, stream)
    failures += total - ok

    report_resources(stream)

    if imports_only:
        print(f"\nResult: {failures} failures (imports-only).", file=stream)
        return 1 if failures else 0

    bin_results = [check_external_binary(b, db_root=db_root) for b in binaries_for_tier(tier)]
    ok, total = _render("External binaries", bin_results, stream)
    # optional-missing was already counted as ok=True in check_external_binary
    failures += sum(1 for r in bin_results if not r.ok)

    db_results = [check_database(d, db_root) for d in databases_for_tier(tier)]
    ok, total = _render("Databases", db_results, stream)
    failures += total - ok

    w_results = [check_weights(w, data_root, db_root) for w in weights_for_tier(tier)]
    ok, total = _render("Model weights", w_results, stream)
    failures += total - ok

    print(file=stream)
    if failures:
        print(f"{failures} check(s) failed.  Address the `fix:` lines above and re-run.", file=stream)
        return 1
    print("All checks passed.", file=stream)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="ssign doctor",
        description="Verify a fresh ssign install: Python packages, external binaries, databases, model weights.",
    )
    p.add_argument(
        "--tier",
        choices=("base", "extended", "full"),
        default=DEFAULT_TIER,
        help=f"Install tier to verify against (default: {DEFAULT_TIER}).",
    )
    p.add_argument(
        "--imports-only",
        action="store_true",
        help="Only check Python imports (skip binaries, DBs, weights). Used by CI.",
    )
    p.add_argument(
        "--data-root",
        default=DEFAULT_DATA_ROOT,
        help=f"Root directory for databases + models (default: {DEFAULT_DATA_ROOT}).",
    )
    args = p.parse_args(argv)
    return run(tier=args.tier, imports_only=args.imports_only, data_root=args.data_root)


if __name__ == "__main__":
    sys.exit(main())
