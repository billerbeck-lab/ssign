#!/usr/bin/env python3
"""Bootstrap pLM-BLAST's embeddings.py with shims for stripped torch wheels.

pLM-BLAST's upstream `embeddings.py` does `import torch.distributed as dist`
and `import torch.multiprocessing as mp` at module load. Both are only
exercised in `mp_process()`, which fires when `nproc > 1` (multi-GPU).
ssign always runs single-GPU (nproc == 1), so the imports do no useful
work for us — but the unconditional import statement still has to
succeed before any embedding code runs.

Some torch wheels ship without those submodules:
  - CPU-only / stripped distributions used on HPC nodes (e.g. CX3 .venv)
  - ARM builds tuned for inference
  - In-house torch builds with USE_DISTRIBUTED=0

When that happens, ProtT5 embedding crashes with
`ModuleNotFoundError: No module named 'torch.distributed'` before we
ever get to embed a single residue. This bootstrap installs an empty
stub at `sys.modules['torch.distributed']` (and likewise for
`torch.multiprocessing`) so the upstream import line resolves, then
exec's `embeddings.py` normally via `runpy`. The stub has no attributes,
so if a multi-GPU run is ever attempted on a stripped torch, the actual
call (e.g. `dist.init_process_group(...)`) raises AttributeError at the
correct failure point — not at import time.

Usage:
    python _plm_blast_embed_bootstrap.py /path/to/pLM-BLAST/embeddings.py \\
        start query.fasta output.pt -embedder pt [--cuda]
"""

from __future__ import annotations

import os
import sys
import types

_STUBBABLE_SUBMODULES = ("torch.distributed", "torch.multiprocessing")


def install_torch_submodule_stubs() -> list[str]:
    """Stub torch.distributed and torch.multiprocessing if absent.

    Returns the list of submodule names that needed stubbing — useful
    for tests and for logging when the bootstrap is doing actual work.
    Returns an empty list when the host torch already ships both
    submodules (the common case on standard pip wheels).
    """
    stubbed: list[str] = []
    try:
        __import__("torch")
    except (ImportError, ModuleNotFoundError):
        return stubbed

    for sub in _STUBBABLE_SUBMODULES:
        try:
            __import__(sub)
            continue
        except (ImportError, ModuleNotFoundError):
            pass
        sys.modules[sub] = types.ModuleType(sub)
        stubbed.append(sub)
    return stubbed


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        print(
            "Usage: _plm_blast_embed_bootstrap.py <path-to-embeddings.py> [embeddings.py args...]",
            file=sys.stderr,
        )
        return 2
    embed_script = argv[0]
    if not os.path.exists(embed_script):
        print(f"ERROR: embeddings.py not found at {embed_script}", file=sys.stderr)
        return 2

    stubbed = install_torch_submodule_stubs()
    if stubbed:
        print(
            f"INFO: stubbed missing torch submodules: {', '.join(stubbed)} "
            f"(stripped torch wheel; only valid for single-GPU embedding)",
            file=sys.stderr,
        )

    script_dir = os.path.dirname(os.path.abspath(embed_script))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    sys.argv = [embed_script] + argv[1:]

    import runpy

    runpy.run_path(embed_script, run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
