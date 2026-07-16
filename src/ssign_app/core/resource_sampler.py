"""In-process resource sampler: writes per-tool resource_samples.csv during a run.

Spawned as a daemon thread by PipelineRunner. Samples system CPU / RAM /
GPU / disk every N seconds and tags each row with the step that's
currently active. The runner sets `current_step` directly, so the sampler
needs no external log-tailing. It runs inside a normal ssign run; the user
does not launch a separate monitoring process.
"""

from __future__ import annotations

import csv
import logging
import os
import threading
import time

logger = logging.getLogger(__name__)


_CSV_HEADER = (
    "timestamp",
    "elapsed_s",
    "cpu_overall_pct",
    "cpu_per_core",
    "ram_used_gb",
    "ram_total_gb",
    "gpu_util_pct",
    "gpu_mem_gb",
    "gpu_mem_total_gb",
    "disk_read_mb_s",
    "disk_write_mb_s",
    "step",
)


class ResourceSampler:
    """Background sampling thread.

    Usage:
        sampler = ResourceSampler(out_path="runtime_data/resource_samples.csv", interval=5.0)
        sampler.start()
        sampler.set_step("parallel: dlp, dse, signalp")
        ...
        sampler.stop()

    Quietly no-ops if psutil isn't importable; GPU columns stay empty
    when pynvml/CUDA aren't available. Either gap leaves a partial CSV
    rather than crashing the pipeline.
    """

    def __init__(self, out_path: str, interval: float = 5.0):
        self.out_path = out_path
        self.interval = max(0.5, interval)
        self._current_step = ""
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()
        self._gpu = None
        self._gpu_total_gb = 0.0
        self._pynvml = None

    def set_step(self, step: str) -> None:
        with self._lock:
            self._current_step = step

    def start(self) -> None:
        try:
            import psutil  # noqa: F401
        except ImportError:
            logger.info("psutil not installed; resource sampling disabled")
            return
        self._init_gpu()
        self._thread = threading.Thread(target=self._run, name="ssign-resource-sampler", daemon=True)
        self._thread.start()

    def stop(self, timeout: float = 5.0) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=timeout)
        # NVML uses ref-counted init/shutdown; pair init with shutdown
        # so we don't leak the driver handle when the runner exits.
        if self._pynvml is not None:
            try:
                self._pynvml.nvmlShutdown()
            except Exception:
                pass
            self._pynvml = None
            self._gpu = None

    def _init_gpu(self) -> None:
        try:
            import pynvml

            pynvml.nvmlInit()
            self._gpu = pynvml.nvmlDeviceGetHandleByIndex(0)
            self._gpu_total_gb = pynvml.nvmlDeviceGetMemoryInfo(self._gpu).total / 2**30
            self._pynvml = pynvml
        except Exception:
            # No GPU, no pynvml, or driver mismatch — GPU columns stay blank.
            self._gpu = None

    def _gpu_sample(self) -> tuple[str, str]:
        if self._gpu is None or self._pynvml is None:
            return ("", "")
        try:
            util = self._pynvml.nvmlDeviceGetUtilizationRates(self._gpu).gpu
            used_gb = self._pynvml.nvmlDeviceGetMemoryInfo(self._gpu).used / 2**30
            return (f"{util:.0f}", f"{used_gb:.2f}")
        except Exception:
            return ("", "")

    def _run(self) -> None:
        import psutil

        os.makedirs(os.path.dirname(os.path.abspath(self.out_path)) or ".", exist_ok=True)
        psutil.cpu_percent(percpu=True)  # prime the per-call delta
        disk_prev = psutil.disk_io_counters()
        t_prev = time.monotonic()
        t_start = t_prev
        ram_total_gb = psutil.virtual_memory().total / 2**30

        with open(self.out_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(_CSV_HEADER)
            fh.flush()
            # Sample immediately, then on every interval thereafter. Sampling
            # first means set_step() called between samples is reliably
            # visible in the row that follows it, instead of swallowed by
            # the initial wait.
            while True:
                try:
                    now = time.monotonic()
                    dt = max(now - t_prev, 1e-6)
                    per_core = psutil.cpu_percent(percpu=True)
                    overall = sum(per_core) / len(per_core) if per_core else 0.0
                    ram_used_gb = psutil.virtual_memory().used / 2**30
                    disk_now = psutil.disk_io_counters()
                    r_mb = (disk_now.read_bytes - disk_prev.read_bytes) / 2**20 / dt
                    wr_mb = (disk_now.write_bytes - disk_prev.write_bytes) / 2**20 / dt
                    disk_prev, t_prev = disk_now, now
                    gu, gm = self._gpu_sample()
                    with self._lock:
                        step = self._current_step
                    writer.writerow(
                        [
                            time.strftime("%Y-%m-%dT%H:%M:%S"),
                            f"{now - t_start:.1f}",
                            f"{overall:.1f}",
                            ";".join(f"{c:.0f}" for c in per_core),
                            f"{ram_used_gb:.2f}",
                            f"{ram_total_gb:.2f}",
                            gu,
                            gm,
                            f"{self._gpu_total_gb:.2f}" if self._gpu is not None else "",
                            f"{r_mb:.1f}",
                            f"{wr_mb:.1f}",
                            step,
                        ]
                    )
                    fh.flush()
                except Exception as e:
                    # A sampling failure must never crash the pipeline.
                    logger.warning("resource sampler row failed: %s", e)
                if self._stop.wait(self.interval):
                    break
