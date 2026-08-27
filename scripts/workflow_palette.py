#!/usr/bin/env python3
"""MetBrewer-inspired publication palettes for PIGSTI workflow figures.

Source: BlakeRMills/MetBrewer — https://github.com/BlakeRMills/MetBrewer
Colorblind-friendly defaults; switch via ACTIVE_PALETTE or set_palette().
"""

from __future__ import annotations

from typing import Any

# Authentication metric chips (shared across figures)
METRICS = [
    "ANI > 96.5%",
    "Entropy ≥ 0.9",
    "Breadth ≥ 0.8",
    "5′ C→T damage",
    "Edit distance",
    "Mapping ratio",
    "Genus rank #1",
]

# Short labels for figure callouts (no thresholds / cutoffs)
AUTH_METRIC_NAMES = [
    "ANI",
    "Entropy",
    "Breadth",
    "5′ C→T damage",
    "Edit distance",
    "Mapping ratio",
    "Genus rank",
]


def _build(prep, host, meta, detect, path, combo, opt, out, tool, ink="#1A1A2E") -> dict[str, str]:
    """Map route accents to fill tints (light washes of each accent)."""
    return {
        "prep": prep,
        "prep_fill": "#EEF4FA",
        "host": host,
        "host_fill": "#E3EEF7",
        "meta": meta,
        "meta_fill": "#EEEAF5",
        "detect": detect,
        "detect_fill": "#FBEEEC",
        "path": path,
        "path_fill": "#F9E8EB",
        "combo": combo,
        "combo_fill": "#FFF6F7",
        "combo_header": path,
        "opt": opt,
        "opt_fill": "#F3F5F8",
        "out": out,
        "out_fill": "#E8F5EC",
        "tool": tool,
        "ink": ink,
        "section_stroke": "#C8D0DC",
        "metric_chip": "#FFFFFF",
        "metric_chip_border": "#E8B4BC",
    }


# Named MetBrewer-inspired route assignments (colorblind-friendly where noted)
PALETTES: dict[str, dict[str, str]] = {
    # Original: Hiroshige blues + Isfahan teal tool
    "hiroshige": _build(
        prep="#1D335F", host="#376795", meta="#5B4A8A", detect="#D4483B",
        path="#8C1C2B", combo="#6B1420", opt="#8A94A6", out="#2E7D52", tool="#3A8F8F",
    ),
    # Egypt (4-color CB) expanded with Isfahan/Hokusai accents
    "egypt": _build(
        prep="#0F7BA2", host="#0F7BA2", meta="#43B284", detect="#DD5129",
        path="#FAB255", combo="#B03A1E", opt="#8A94A6", out="#2E7D52", tool="#178F92",
        ink="#1A1A2E",
    ),
    # Johnson (CB) — bold primaries
    "johnson": _build(
        prep="#132B69", host="#0086A8", meta="#F6C200", detect="#D04E00",
        path="#A00E00", combo="#7A0A00", opt="#8A94A6", out="#2E7D52", tool="#0086A8",
    ),
    # Derain (CB) — greens / violets / gold
    "derain": _build(
        prep="#454A74", host="#5C66A8", meta="#6F9969", detect="#EFC86E",
        path="#808FE1", combo="#454A74", opt="#8A94A6", out="#6F9969", tool="#97C684",
    ),
    # Hokusai3 (CB) — cool greens → deep blue
    "hokusai3": _build(
        prep="#0A2E57", host="#295384", meta="#5A97C1", detect="#95C36E",
        path="#74C8C3", combo="#0A2E57", opt="#8A94A6", out="#95C36E", tool="#74C8C3",
    ),
    # Java (CB) — magenta / coral / teal
    "java": _build(
        prep="#663171", host="#0C7156", meta="#CF3A36", detect="#EA7428",
        path="#E2998A", combo="#663171", opt="#8A94A6", out="#0C7156", tool="#0C7156",
    ),
    # Archambault (CB)
    "archambault": _build(
        prep="#381A61", host="#88A0DC", meta="#7C4B73", detect="#ED968C",
        path="#AB3329", combo="#381A61", opt="#8A94A6", out="#E78429", tool="#88A0DC",
    ),
    # Cassatt1 — soft rose / slate / teal
    "cassatt1": _build(
        prep="#3B4256", host="#5F6A7A", meta="#B395BD", detect="#C17E70",
        path="#8B4A3B", combo="#5C2E2A", opt="#9AA3B2", out="#4F7F6E", tool="#6B8F9A",
    ),
    # Monet — waterlily greens / blues / peach
    "monet": _build(
        prep="#1E4D5C", host="#3A7A7A", meta="#6B8F71", detect="#D17A4A",
        path="#A34E3B", combo="#6E2F28", opt="#8FA3A8", out="#4A7A5C", tool="#5A8F9E",
    ),
    # VanGogh1 — wheat / olive / sky
    "vangogh1": _build(
        prep="#2B3A4A", host="#4A6FA5", meta="#7A8F3D", detect="#E0A03A",
        path="#C45C2A", combo="#7A3018", opt="#9AA3A8", out="#5A7A3A", tool="#4A7A8A",
    ),
    # Klimt — gold / deep red / teal
    "klimt": _build(
        prep="#2C2A29", host="#1F5C6B", meta="#C4A35A", detect="#B83A2E",
        path="#7A1F1A", combo="#4A1210", opt="#8A8A84", out="#2E6B5A", tool="#1F5C6B",
    ),
    # Greek — Aegean blue / terracotta / olive
    "greek": _build(
        prep="#1B3A5F", host="#2F6FAE", meta="#6B7A3A", detect="#C45C3A",
        path="#8B2E1F", combo="#5A1A12", opt="#8A94A0", out="#3A6B4A", tool="#2F7A8A",
    ),
    # Navajo — sand / turquoise / clay
    "navajo": _build(
        prep="#3D2B1F", host="#1F6B6B", meta="#C4A06A", detect="#C45A2A",
        path="#8B3A1A", combo="#5A2410", opt="#A09080", out="#2E6B4A", tool="#1F6B6B",
    ),
    # Thomas — ink / indigo / coral
    "thomas": _build(
        prep="#1A1F36", host="#3B4F8A", meta="#6A5A8A", detect="#D95A4A",
        path="#9A2E3A", combo="#5A1520", opt="#8A90A0", out="#2E6B5A", tool="#3B6F8A",
    ),
    # Isfahan1 — teal / saffron / wine
    "isfahan1": _build(
        prep="#1A3A4A", host="#1F7A7A", meta="#4A6B8A", detect="#D48A2A",
        path="#8B2E3A", combo="#5A1520", opt="#8A94A0", out="#2E6B5A", tool="#1F7A7A",
    ),
    # Rainbow — Gilbert Baker / pride-flag inspired (6 classic stripes)
    # red → orange → yellow → green → blue → violet
    "rainbow": _build(
        prep="#750787",   # violet
        host="#004CFF",   # blue
        meta="#008026",   # green
        detect="#FFD100", # yellow (gold — readable on white)
        path="#FF8C00",   # orange
        combo="#E40303",  # red
        opt="#8A94A6",
        out="#008026",
        tool="#4C1D95",
    ),
}

ACTIVE_PALETTE = "hiroshige"
PAL: dict[str, str] = dict(PALETTES[ACTIVE_PALETTE])


def set_palette(name: str) -> dict[str, str]:
    """Activate a named palette; returns the palette dict."""
    global ACTIVE_PALETTE, PAL
    if name not in PALETTES:
        raise KeyError(f"Unknown palette {name!r}; choose from {sorted(PALETTES)}")
    ACTIVE_PALETTE = name
    PAL.clear()
    PAL.update(PALETTES[name])
    return PAL


def list_palettes() -> list[str]:
    return sorted(PALETTES)


def palette_legend_colors(pal: dict[str, str] | None = None) -> list[tuple[str, str, bool]]:
    p = pal or PAL
    return [
        (p["prep"], "Preprocess", False),
        (p["host"], "Host", False),
        (p["meta"], "Metagenomics", False),
        (p["detect"], "Detection", False),
        (p["path"], "Mapping", False),
        (p["combo"], "Authentication", False),
        (p["opt"], "Optional", True),
    ]
