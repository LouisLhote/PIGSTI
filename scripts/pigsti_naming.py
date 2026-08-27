"""
Canonical pathogen / sample naming helpers for PIGSTI.

Keep in sync with Snakefile safe_name() and docs/OUTPUT_SCHEMA.md.
"""

from __future__ import annotations


def safe_pathogen_name(name: str) -> str:
    """Filesystem-safe pathogen token used in results paths."""
    s = str(name).strip()
    for ch in (" ", "/", "\\", ":"):
        s = s.replace(ch, "_")
    while "__" in s:
        s = s.replace("__", "_")
    return s


def hops_species_token(name: str) -> str:
    """HOPS heatmap row labels use spaces → underscores."""
    return str(name).strip().replace(" ", "_")


def normalize_taxonomy_key(name: str) -> str:
    """Case-insensitive key for matching E-score taxonomy to spreadsheet."""
    return str(name).lower().replace("_", " ").strip()
