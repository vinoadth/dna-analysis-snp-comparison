from __future__ import annotations

import logging
import os
import time
from contextlib import contextmanager

logger = logging.getLogger(__name__)


def analyze_timings_enabled(settings_flag: bool | None = None) -> bool:
    if settings_flag is not None:
        return bool(settings_flag)
    return os.environ.get("DNA_COMPARE_TIMINGS", "").lower() in {"1", "true", "yes"}


class AnalyzeTimings:
    def __init__(self, enabled: bool = False):
        self.enabled = enabled
        self.stages: dict[str, float] = {}
        self.native_used: dict[str, bool] = {}

    @contextmanager
    def stage(self, name: str):
        if not self.enabled:
            yield
            return
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            self.stages[name] = self.stages.get(name, 0.0) + elapsed_ms
            logger.info("analyze timing %s: %.1f ms", name, elapsed_ms)

    def mark_native(self, name: str, used: bool) -> None:
        self.native_used[name] = used

    def to_dict(self) -> dict:
        return {
            "stages_ms": {k: round(v, 2) for k, v in sorted(self.stages.items())},
            "native": dict(self.native_used),
            "total_ms": round(sum(self.stages.values()), 2),
        }
