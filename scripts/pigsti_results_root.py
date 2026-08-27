"""Resolve PIGSTI results directory (default: ./results)."""

from __future__ import annotations

import os
from pathlib import Path


def results_root_from_config(cfg: dict | None = None) -> str:
    """Read results_root from config dict, env PIGSTI_RESULTS_ROOT, or .pigsti/results_root."""
    if cfg is None:
        cfg = {}
    raw = (
        cfg.get("results_root")
        or cfg.get("results_dir")
        or os.environ.get("PIGSTI_RESULTS_ROOT")
        or "results"
    )
    root = os.path.abspath(os.path.expanduser(str(raw).strip()))
    marker = Path(".pigsti/results_root")
    if marker.is_file():
        try:
            return os.path.abspath(marker.read_text(encoding="utf-8").strip())
        except OSError:
            pass
    return root


def join_results(*parts: str) -> str:
    """Build a path under the active results root."""
    return "/".join([results_root_from_config().rstrip("/")] + [p.strip("/") for p in parts if p])
