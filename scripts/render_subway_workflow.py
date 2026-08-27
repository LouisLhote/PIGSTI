#!/usr/bin/env python3
"""
PIGSTI subway workflow — nf-metro **light** aesthetic.

Primary: nf-metro render (theme/mode light)
Fallback: matplotlib matching the same light tokens + PNG icons.
"""

from __future__ import annotations

import argparse
import base64
import shutil
import subprocess
import sys
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Circle, FancyBboxPatch
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
MMD = ROOT / "docs" / "images" / "pigsti_subway_workflow.mmd"
ICONS = ROOT / "docs" / "images" / "Icons" / "png"
MICROBES = ROOT / "docs" / "images" / "Icons microbes" / "png"
# All host-animal silhouettes available under docs/images/Icons/png
HOST_ANIMALS = [
    "Cattle", "Sheep", "Goat", "pig", "Dog",
    "Horse", "Cat", "Camel", "Deer", "Rat",
]
sys.path.insert(0, str(ROOT / "scripts"))
from workflow_palette import AUTH_METRIC_NAMES, METRICS, PAL, list_palettes, set_palette  # noqa: E402

# nf-metro light tokens
LIGHT = {
    "bg": "#FFFFFF",
    "section_fill": "#EDEDED",
    "section_stroke": "#E6E6E6",
    "section_label": "#666666",
    "num": "#333333",
    "ink": "#111111",
    "label": "#333333",
    "marker": "#333333",
    "tool": "#4B5563",
    "tool_bg": "#F3F4F6",
    "tool_stroke": "#D1D5DB",
    "muted": "#888888",
}

W, H = 1880, 1180
TRACK = 11.0
TRACK_OPT = 5.0  # thinner so dashes stay readable at print size
# Dash units in points — must scale with linewidth or gaps vanish under round caps
DASH_SOLID = "-"
# on/off in points — short enough to show on MALT branch, long enough to read as dashes
DASH_OPT = (0, (9.0, 8.0))

Yp, Yh, Yh2, Yh3, Ym, Yo, Yd, Yout = 130, 245, 410, 490, 580, 720, 170, 310
X = {
    "raw": 90, "trim": 240,
    # host: Damage (upper) + Mapping QC (lower) metro loop; sexing forks from Mapping QC
    "hid": 90, "hmap": 300, "mtdna": 300, "hmerge": 470,
    "hdmg": 720, "hqual": 720, "softclip": 880, "sexing": 980,
    "prinseq": 90, "unal": 260, "pool": 430, "kraken": 600, "malt": 780,
    "detect": 880, "pmap": 1120, "auth": 1340, "final": 1580,
}


def place_icon(ax, x, y, path: Path, half: float = 28, z: int = 8):
    """Place a transparent PNG silhouette centred at (x, y) via imshow (correct alpha)."""
    if not path.exists():
        return
    import numpy as np
    arr = np.asarray(Image.open(path).convert("RGBA"))
    ax.imshow(
        arr,
        extent=(x - half, x + half, y + half, y - half),
        origin="upper",
        interpolation="hanning",
        zorder=z,
        aspect="equal",
        clip_on=False,
    )


def load_icon(path: Path, zoom: float = 0.35):
    """Deprecated helper kept for compatibility — prefer place_icon(path=...)."""
    return path

def draw_track(ax, points, color, lw=TRACK, dashed=False, z=1):
    xs, ys = zip(*points)
    if dashed:
        # White underlay + coloured dashes (no path_effects — they close gaps on thick dashes)
        ax.plot(
            xs, ys, color="white", lw=lw + 3.0, solid_capstyle="butt", solid_joinstyle="miter",
            ls=DASH_OPT, zorder=z, alpha=0.95,
        )
        ax.plot(
            xs, ys, color=color, lw=lw, solid_capstyle="butt", solid_joinstyle="miter",
            ls=DASH_OPT, zorder=z + 0.1,
        )
    else:
        ax.plot(
            xs, ys, color=color, lw=lw, solid_capstyle="round", solid_joinstyle="round",
            ls=DASH_SOLID, zorder=z,
            path_effects=[pe.Stroke(linewidth=lw + 3.5, foreground="white", alpha=0.95), pe.Normal()],
        )


def section_panel(ax, x, y, w, h, title, n):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=10",
        fc=LIGHT["section_fill"], ec=LIGHT["section_stroke"], lw=1.0, zorder=0,
    ))
    ax.add_patch(Circle((x + 16, y - 14), 11, fc=LIGHT["num"], ec="none", zorder=2))
    ax.text(x + 16, y - 14, str(n), ha="center", va="center", fontsize=10,
            fontweight="bold", color="#FFF", zorder=3)
    ax.text(x + 34, y - 14, title, ha="left", va="center", fontsize=13,
            fontweight="bold", color=LIGHT["section_label"], zorder=3,
            fontfamily="sans-serif")


def _label_block(ax, x: float, y_anchor: float, purpose: str, tool: str, *, side: str, scale: float = 1.0) -> None:
    """Station label: purpose (bold) + tool (italic), no plate/box."""
    purpose_fs, tool_fs = 14.5 * scale, 12.0 * scale
    gap = 4.0 * scale

    # reading order: purpose first, tool second; side chooses above/below the station
    if side == "up":
        y_tool = y_anchor - 20 * scale
        y_purpose = y_tool - (tool_fs + gap)
    else:
        y_purpose = y_anchor + 20 * scale
        y_tool = y_purpose + (purpose_fs + gap)

    ax.text(
        x, y_purpose, purpose, ha="center", va="center",
        fontsize=purpose_fs, fontweight="bold", color="#0F172A",
        fontfamily="DejaVu Sans", zorder=14, clip_on=False,
    )
    if tool:
        ax.text(
            x, y_tool, f"({tool})", ha="center", va="center",
            fontsize=tool_fs, color="#475569",
            fontfamily="DejaVu Sans", zorder=14, clip_on=False,
        )


def station(ax, x, y, purpose, tool="", *, kind="node", side="up"):
    """Station marker (solid circle) + purpose / tool labels."""
    mk = LIGHT["marker"]
    # Terminal stations (input FASTQ / deliverables) are larger for emphasis
    if kind in {"start", "end"}:
        r, lw, scale = 18.0, 3.8, 1.5
    else:
        r, lw, scale = 9.0, 2.6, 1.0
    ax.add_patch(Circle((x, y), r, fc="#FFF", ec=mk, lw=lw, zorder=5))
    _label_block(ax, x, y, purpose, tool, side=side, scale=scale)


def draw_auth_details_box(ax, cx: float, top_y: float) -> None:
    """Nicer authentication callout: titled card with metric chips (names only)."""
    from matplotlib.patches import Rectangle

    names = AUTH_METRIC_NAMES
    cols, chip_h, gap_x, gap_y = 4, 26.0, 10.0, 10.0
    pad_x, pad_y, title_h = 16.0, 14.0, 30.0
    # estimate chip widths from label length
    chip_ws = [max(58.0, 8.2 * len(n) + 24) for n in names]
    rows = (len(names) + cols - 1) // cols
    row_ws = []
    for r in range(rows):
        chunk = chip_ws[r * cols : (r + 1) * cols]
        row_ws.append(sum(chunk) + gap_x * (len(chunk) - 1))
    inner_w = max(row_ws) if row_ws else 200
    box_w = inner_w + 2 * pad_x
    box_h = title_h + pad_y + rows * chip_h + (rows - 1) * gap_y + pad_y
    box_x = cx - box_w / 2
    box_y = top_y

    accent = PAL["combo"]
    ax.add_patch(Rectangle(
        (box_x, box_y), box_w, box_h,
        facecolor="#FFFFFF", edgecolor=accent, linewidth=2.4,
        zorder=20, clip_on=False, joinstyle="miter",
    ))
    ax.add_patch(Rectangle(
        (box_x, box_y), box_w, title_h,
        facecolor=accent, edgecolor="none", linewidth=0,
        zorder=21, clip_on=False,
    ))
    ax.text(
        cx, box_y + title_h / 2, "Authentication details",
        ha="center", va="center", fontsize=12.0, fontweight="bold",
        color="#FFFFFF", fontfamily="DejaVu Sans", zorder=22, clip_on=False,
    )

    y0 = box_y + title_h + pad_y
    for i, name in enumerate(names):
        r, c = divmod(i, cols)
        row_chunk = chip_ws[r * cols : (r + 1) * cols]
        row_w = sum(row_chunk) + gap_x * (len(row_chunk) - 1)
        x0 = cx - row_w / 2
        x = x0 + sum(row_chunk[:c]) + gap_x * c
        w = chip_ws[i]
        y = y0 + r * (chip_h + gap_y)
        ax.add_patch(FancyBboxPatch(
            (x, y), w, chip_h,
            boxstyle="round,pad=0.02,rounding_size=5",
            fc="#F8FAFC", ec="#CBD5E1", lw=1.2, zorder=21, clip_on=False,
        ))
        ax.text(
            x + w / 2, y + chip_h / 2, name, ha="center", va="center",
            fontsize=9.5, fontweight="bold", color="#0F172A",
            fontfamily="DejaVu Sans", zorder=22, clip_on=False,
        )



def render_matplotlib(outdir: Path, stem: str) -> None:
    fig, ax = plt.subplots(figsize=(19.2, 12.0), dpi=120)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor(LIGHT["bg"])
    ax.set_facecolor(LIGHT["bg"])

    ax.text(40, 36, "PIGSTI", fontsize=22, fontweight="bold", color=LIGHT["ink"],
            fontfamily="sans-serif")

    section_panel(ax, 30, 70, 280, 95, "Preprocessing", 1)
    section_panel(ax, 30, 185, 1050, 320, "Host route", 2)
    section_panel(ax, 30, 520, 820, 260, "Metagenomics screening", 3)
    section_panel(ax, 860, 70, 980, 300, "Pathogen route", 4)

    # Tracks
    draw_track(ax, [(X["raw"], Yp), (X["trim"], Yp)], PAL["prep"])
    draw_track(ax, [(X["trim"], Yp), (X["trim"], Yh), (X["hid"], Yh)], PAL["host"])
    # Host main line through Damage (upper rail of QC loop), then Soft clipping
    draw_track(
        ax,
        [
            (X["hid"], Yh), (X["hmap"], Yh), (X["hmerge"], Yh),
            (X["hdmg"], Yh), (X["softclip"], Yh),
        ],
        PAL["host"],
    )
    draw_track(ax, [(X["hid"], Yh2), (X["mtdna"], Yh2), (X["hmerge"], Yh2), (X["hmerge"], Yh)], PAL["host"])
    # Mapping QC below — lower rail of metro loop, rejoins at Soft clipping
    draw_track(
        ax,
        [
            (X["hmerge"], Yh), (X["hqual"], Yh2),
            (X["softclip"], Yh2), (X["softclip"], Yh),
        ],
        PAL["host"],
    )
    # Genetic sexing — optional fork from Mapping QC
    draw_track(
        ax,
        [(X["hqual"], Yh2), (X["hqual"], Yh3), (X["sexing"], Yh3)],
        PAL["opt"], TRACK_OPT, dashed=True,
    )
    draw_track(ax, [(X["trim"], Yp), (X["trim"], Ym), (X["prinseq"], Ym)], PAL["meta"])
    draw_track(ax, [(X["prinseq"], Ym), (X["unal"], Ym), (X["pool"], Ym), (X["kraken"], Ym)], PAL["meta"])
    # MALT optional spur starts from pool, not under Classify
    draw_track(ax, [(X["pool"], Ym), (X["malt"], Yo)], PAL["opt"], TRACK_OPT, dashed=True)
    draw_track(ax, [(X["kraken"], Ym), (X["detect"], Yd)], PAL["detect"])
    draw_track(ax, [(X["detect"], Yd), (X["pmap"], Yd)], PAL["path"])
    draw_track(ax, [(X["pmap"], Yd), (X["auth"], Yd), (X["auth"], Yout), (X["final"], Yout)], PAL["combo"])
    draw_track(ax, [(X["malt"], Yo), (X["detect"] - 50, Yo), (X["detect"], Yd)], PAL["opt"], TRACK_OPT, dashed=True)
    draw_track(ax, [(X["softclip"], Yh), (X["softclip"] + 40, Yh), (X["final"] - 80, Yout), (X["final"], Yout)], PAL["host"])

    # purpose / tool (sides alternate so plates do not overlap)
    station(ax, X["raw"], Yp, "Raw reads", "FASTQ", kind="start", side="up")
    station(ax, X["trim"], Yp, "Trim adapters", "AdapterRemoval", side="down")
    station(ax, X["hid"], Yh, "Host identification", "FastQ Screen", side="up")
    station(ax, X["hmap"], Yh, "Map host", "BWA / Bowtie2", side="down")
    station(ax, X["mtdna"], Yh2, "Map mtDNA", "BWA / Bowtie2", side="up")
    station(ax, X["hmerge"], Yh, "Merge samples", "samtools", kind="join", side="up")
    station(ax, X["hdmg"], Yh, "Damage", "DamageProfiler", side="up")
    station(ax, X["hqual"], Yh2, "Mapping QC", "Qualimap", side="down")
    station(ax, X["softclip"], Yh, "Soft clipping", "Python", side="down")
    station(ax, X["sexing"], Yh3, "Genetic sexing", "sexing_residual", kind="opt", side="down")
    station(ax, X["prinseq"], Ym, "Dedup", "PRINSEQ++", side="up")
    station(ax, X["unal"], Ym, "Mammals removal", "Bowtie2", side="down")
    station(ax, X["pool"], Ym, "Pool samples", "cat · pigz", kind="join", side="up")
    station(ax, X["kraken"], Ym, "Metagenome screening", "KrakenUniq", side="up")
    station(ax, X["malt"], Yo, "Metagenome screening", "MALT", kind="opt", side="down")
    station(ax, X["detect"], Yd, "Pathogen detection", "E-value · MaltExtract", side="up")
    station(ax, X["pmap"], Yd, "Pathogen mapping", "BWA / Bowtie2", side="up")
    station(ax, X["auth"], Yd, "Authentication", "composite score", side="down")
    station(ax, X["final"], Yout, "Deliverables", "Excel · PDF", kind="end", side="up")

    # Authentication details box below the Authentication label
    draw_auth_details_box(ax, X["auth"], Yd + 78)

    # Host species card — all animal icons (5×2)
    animals = [a for a in HOST_ANIMALS if (ICONS / f"{a}.png").exists()]
    cols, rows = 5, 2
    icon_half, gap, pad, title_h = 14, 5, 8, 13
    card_w = 2 * pad + cols * (2 * icon_half) + (cols - 1) * gap
    card_h = title_h + 2 * pad + rows * (2 * icon_half) + (rows - 1) * gap
    card_x = X["hmerge"] + 40
    card_y = Yh + 55
    ax.add_patch(FancyBboxPatch(
        (card_x, card_y), card_w, card_h,
        boxstyle="round,pad=0.02,rounding_size=10",
        fc="#FFFFFF", ec=PAL["host"], lw=1.8, zorder=6,
    ))
    ax.text(
        card_x + card_w / 2, card_y + 8, "Host species",
        ha="center", va="center", fontsize=7.0, fontweight="bold",
        color=PAL["host"], zorder=7,
    )
    for i, name in enumerate(animals):
        r, c = divmod(i, cols)
        cx = card_x + pad + icon_half + c * (2 * icon_half + gap)
        cy = card_y + title_h + pad + icon_half + r * (2 * icon_half + gap)
        place_icon(ax, cx, cy, ICONS / f"{name}.png", half=icon_half, z=8)

    place_icon(ax, 780, 555, MICROBES / "E.coli.png", half=32)
    place_icon(ax, 1720, 115, MICROBES / "virus.png", half=38)
    place_icon(ax, X["detect"] - 48, Yd - 52, MICROBES / "B.melitensis.png", half=26)
    place_icon(ax, X["detect"] + 20, Yd - 52, MICROBES / "virus.png", half=26)

    from matplotlib.patches import Rectangle

    # Route colour key — bottom-right (no title)
    legend_items = [
        (PAL["prep"], "Preprocessing", False),
        (PAL["host"], "Host", False),
        (PAL["meta"], "Metagenomics", False),
        (PAL["detect"], "Detection", False),
        (PAL["path"], "Mapping", False),
        (PAL["combo"], "Authentication", False),
        (PAL["opt"], "Optional", True),
    ]
    row_h = 20
    cols = 2
    rows = (len(legend_items) + cols - 1) // cols
    col_w = 175
    pad_top, pad_bot = 16, 14
    leg_w = 24 + cols * col_w
    leg_h = pad_top + rows * row_h + pad_bot
    leg_x = W - leg_w - 36
    leg_y = H - leg_h - 28

    ax.add_patch(Rectangle(
        (leg_x, leg_y), leg_w, leg_h,
        facecolor="#FFFFFF", edgecolor="none", linewidth=0,
        zorder=50, clip_on=False,
    ))
    for i, (col, label, dashed) in enumerate(legend_items):
        rr, cc = divmod(i, cols)
        lx = leg_x + 16 + cc * col_w
        ly = leg_y + pad_top + rr * row_h
        if dashed:
            ax.plot(
                [lx, lx + 28], [ly, ly], color=col, lw=3.2,
                ls=(0, (1.8, 1.3)), solid_capstyle="butt",
                zorder=52, clip_on=False,
            )
        else:
            ax.plot(
                [lx, lx + 28], [ly, ly], color=col, lw=7.0,
                solid_capstyle="butt", zorder=52, clip_on=False,
            )
        ax.text(
            lx + 36, ly, label, ha="left", va="center", fontsize=9.5,
            color="#0F172A", fontfamily="DejaVu Sans", zorder=52, clip_on=False,
        )
    ax.add_patch(Rectangle(
        (leg_x, leg_y), leg_w, leg_h,
        facecolor="none", edgecolor="#000000", linewidth=2.8,
        zorder=53, clip_on=False, joinstyle="miter",
    ))

    outdir.mkdir(parents=True, exist_ok=True)
    for ext, dpi in (("png", 320), ("pdf", 320), ("svg", 120)):
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path, dpi=dpi, facecolor=LIGHT["bg"], bbox_inches="tight", pad_inches=0.25)
        print(f"Wrote {path}")
    plt.close(fig)


def write_mmd(outdir: Path, stem: str) -> Path:
    text = f"""%%metro title: PIGSTI
%%metro style: light
%%metro file: fastq_in | FASTQ
%%metro file: deliverables | Reports
%%metro line: shared | Preprocessing | {PAL['prep']}
%%metro line: host | Host route | {PAL['host']}
%%metro line: meta | Metagenomics screening | {PAL['meta']}
%%metro line: detect | Pathogen detection | {PAL['detect']}
%%metro line: path | Pathogen mapping | {PAL['path']}
%%metro line: auth | Authentication | {PAL['combo']}
%%metro line: opt | Optional | {PAL['opt']} | dashed
%%metro legend: br

graph LR
 subgraph preprocessing [Preprocessing]
 %%metro exit: right | shared, host, meta
 fastq_in[ ]
 trim[Trim adapters\\nAdapterRemoval]

 fastq_in -->|shared,host,meta| trim
 end

 subgraph host_route [Host route]
 %%metro entry: left | host
 host_id[Host identification\\nFastQ Screen]
 hmap[Map host\\nBWA / Bowtie2]
 mtdna[Map mtDNA\\nBWA / Bowtie2]
 hmerge[Merge samples\\nsamtools]
 hdmg[Damage\\nDamageProfiler]
 hqual[Mapping QC\\nQualimap]
 softclip[Soft clipping\\nPython]
 sexing[Genetic sexing\\nsexing_residual]

 trim -->|host| host_id
 host_id -->|host| hmap
 host_id -->|host| mtdna
 hmap -->|host| hmerge
 mtdna -->|host| hmerge
 hmerge -->|host| hqual
 hmerge -->|host| hdmg
 hqual -->|host| softclip
 hdmg -->|host| softclip
 hqual -->|opt| sexing
 end

 subgraph metagenomics_screening [Metagenomics screening]
 %%metro entry: left | meta, opt
 %%metro exit: right | meta, detect, opt
 prinseq[Dedup\\nPRINSEQ++]
 unaligned[Mammals reads removal\\nBowtie2]
 pool[Pool samples\\ncat · pigz]
 kraken[Metagenome screening\\nKrakenUniq]
 malt[Metagenome screening\\nMALT]
 decom[Source tracking\\ndecOM]

 trim -->|meta| prinseq
 prinseq -->|meta| unaligned
 unaligned -->|meta| pool
 pool -->|meta| kraken
 pool -->|opt| malt
 kraken -->|opt| decom
 end

 subgraph pathogen_route [Pathogen route]
 %%metro direction: LR
 %%metro entry: left | detect, opt
 %%metro exit: right | auth, host
 detect[Pathogen detection\\nE-value · MaltExtract]
 pmap[Pathogen mapping\\nBWA / Bowtie2]
 auth[Authentication\\ncomposite score]
 deliverables[ ]

 kraken -->|detect| detect
 malt -->|opt| detect
 detect -->|path| pmap
 pmap -->|auth| auth
 auth -->|auth| deliverables
 softclip -->|host| deliverables
 end
"""
    path = outdir / f"{stem}.mmd"
    path.write_text(text, encoding="utf-8")
    if stem != "pigsti_subway_workflow":
        # also keep canonical when rendering default
        pass
    return path


def png_data_uri(path: Path) -> str | None:
    if not path.exists():
        return None
    b64 = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def thicken_nfmetro_strokes(
    svg_text: str,
    *,
    solid: float = 10.0,
    dashed: float = 5.0,
    dash_array: str = "10 9",
) -> str:
    """Bump nf-metro track stroke widths; force readable dash gaps on optional lines."""
    import re

    def _repl(m: re.Match[str]) -> str:
        tag = m.group(0)
        is_dashed = "stroke-dasharray" in tag
        width = dashed if is_dashed else solid
        if re.search(r'stroke-width="[0-9.]+"', tag):
            tag = re.sub(r'stroke-width="[0-9.]+"', f'stroke-width="{width}"', tag)
        else:
            tag = tag.replace("/>", f' stroke-width="{width}"/>') if tag.endswith("/>") else tag
        if is_dashed:
            tag = re.sub(r'stroke-dasharray="[^"]*"', f'stroke-dasharray="{dash_array}"', tag)
            if 'stroke-linecap=' in tag:
                tag = re.sub(r'stroke-linecap="[^"]*"', 'stroke-linecap="butt"', tag)
            else:
                tag = tag.replace("/>", ' stroke-linecap="butt"/>') if tag.endswith("/>") else tag
        return tag

    return re.sub(r"<path\b[^>]*class=\"metro-line-[^\"]*\"[^>]*/?>", _repl, svg_text)


def inject_icons_into_nfmetro_svg(svg_path: Path) -> None:
    """Thicken tracks + replace labels/legend to match the main figure."""
    import re

    text = svg_path.read_text(encoding="utf-8")
    if "pigsti-icon" in text:
        return  # already injected

    text = thicken_nfmetro_strokes(text)

    def _station_plate(x: float, y: float, purpose: str, tool: str, *, side: str) -> str:
        purpose_esc = escape(purpose)
        tool_txt = f"({tool})" if tool else ""
        tool_txt_esc = escape(tool_txt)
        purpose_fs, tool_fs, gap = 14.5, 12.0, 4.0
        if side == "up":
            y_tool = y - 20
            y_purpose = y_tool - (tool_fs + gap)
        else:
            y_purpose = y + 20
            y_tool = y_purpose + (purpose_fs + gap)
        parts = [
            f'<text class="pigsti-icon" x="{x}" y="{y_purpose}" text-anchor="middle" dominant-baseline="middle" '
            f'font-size="{purpose_fs}" font-weight="bold" fill="#0F172A" '
            f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">{purpose_esc}</text>',
        ]
        if tool_txt:
            parts.append(
                f'<text class="pigsti-icon" x="{x}" y="{y_tool}" text-anchor="middle" dominant-baseline="middle" '
                f'font-size="{tool_fs}" fill="#475569" '
                f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">{tool_txt_esc}</text>'
            )
        return "".join(parts)

    def node_xy(node_id: str) -> tuple[float, float] | None:
        m = re.search(
            rf'data-node-id="{node_id}"[^>]*data-node-cx="([0-9.]+)"[^>]*data-node-cy="([0-9.]+)"',
            text,
        )
        if not m:
            m = re.search(
                rf'data-node-id="{node_id}"[^>]*data-node-cy="([0-9.]+)"[^>]*data-node-cx="([0-9.]+)"',
                text,
            )
            if m:
                return float(m.group(2)), float(m.group(1))
            return None
        return float(m.group(1)), float(m.group(2))

    def section_box(section_id: str) -> tuple[float, float, float, float] | None:
        m = re.search(rf'<rect[^>]*data-section-id="{section_id}"[^>]*/?>', text)
        if not m:
            return None
        tag = m.group(0)
        nums = {k: float(v) for k, v in re.findall(r'(x|y|width|height)="([0-9.]+)"', tag)}
        if not {"x", "y", "width", "height"} <= nums.keys():
            return None
        return nums["x"], nums["y"], nums["width"], nums["height"]

    def hide_default_labels_and_legend(svg_text: str) -> str:
        svg_text = re.sub(
            r'<text[^>]*class="nf-metro-label-halo"[^>]*>.*?</text>\s*',
            "",
            svg_text,
            flags=re.S,
        )
        svg_text = re.sub(
            r'<text[^>]*class="nf-metro-station-label"[^>]*>.*?</text>\s*',
            "",
            svg_text,
            flags=re.S,
        )
        svg_text = re.sub(
            r'<rect[^>]*class="nf-metro-legend-bg"[^>]*/>\s*'
            r'(?:<line[^>]*>\s*</line>\s*|<line[^>]*/>\s*)*'
            r'(?:<path[^>]*>\s*</path>\s*|<path[^>]*/>\s*)*'
            r'(?:<text[^>]*class="nf-metro-legend-text"[^>]*>.*?</text>\s*)+',
            "",
            svg_text,
            flags=re.S,
        )
        return svg_text

    def circles_for_stations(svg_text: str) -> str:
        """Replace nf-metro rect/capsule stations with solid full circles."""

        def _to_circle(m: re.Match[str]) -> str:
            tag = m.group(0)
            nums = {
                k: float(v)
                for k, v in re.findall(r'(?<![-\w])(x|y|width|height)="([0-9.]+)"', tag)
            }
            if not {"x", "y", "width", "height"} <= nums.keys():
                return tag
            cx = nums["x"] + nums["width"] / 2
            cy = nums["y"] + nums["height"] / 2
            # preserve identity attributes
            attrs = []
            sid = ""
            for key in ("class", "data-station-id", "data-station-lines",
                        "data-station-label", "data-section-id", "data-section-name"):
                am = re.search(rf'\b{key}="([^"]*)"', tag, flags=re.S)
                if am:
                    attrs.append(f'{key}="{am.group(1)}"')
                    if key == "data-station-id":
                        sid = am.group(1)
            extra = (" " + " ".join(attrs)) if attrs else ""
            # Terminal FASTQ / Reports markers are larger
            if sid in {"fastq_in", "deliverables"}:
                r, lw = 12.0, 3.8
            else:
                r, lw = 6.0, 2.6
            return (
                f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#ffffff" '
                f'stroke="#333333" stroke-width="{lw}"{extra}/>'
            )

        return re.sub(
            r'<rect\b[^>]*?\bclass="[^"]*nf-metro-station[^"]*"[^>]*?/>',
            _to_circle,
            svg_text,
            flags=re.S,
        )

    def enlarge_terminal_file_badges(svg_text: str) -> str:
        """Scale FASTQ / Reports file badges and bump their label font size."""

        def _scale_badge(m: re.Match[str]) -> str:
            sid, body = m.group(1), m.group(2)
            if sid not in {"fastq_in", "deliverables"}:
                return m.group(0)
            # Scale badge geometry around the station centre
            cx_m = re.search(rf'data-node-id="{sid}"[^>]*data-node-cx="([0-9.]+)"', svg_text)
            cy_m = re.search(rf'data-node-id="{sid}"[^>]*data-node-cy="([0-9.]+)"', svg_text)
            if not cx_m or not cy_m:
                return m.group(0)
            cx, cy = float(cx_m.group(1)), float(cy_m.group(1))
            body2 = re.sub(
                r'font-size="([0-9.]+)"',
                lambda fm: f'font-size="{float(fm.group(1)) * 1.55:.2f}"',
                body,
            )
            return (
                f'<g data-station-id="{sid}" transform="translate({cx},{cy}) '
                f'scale(1.85) translate({-cx},{-cy})">{body2}</g>'
            )

        return re.sub(
            r'<g data-station-id="(fastq_in|deliverables)">(.*?)</g>',
            _scale_badge,
            svg_text,
            flags=re.S,
        )

    text = hide_default_labels_and_legend(text)
    text = circles_for_stations(text)
    # Also enlarge terminal circles when nf-metro already emitted <circle> stations
    text = re.sub(
        r'(<circle\b[^>]*\bdata-station-id="(?:fastq_in|deliverables)"[^>]*\br=")(?:6\.0|10\.0)(")',
        r'\g<1>12.0\2',
        text,
    )
    text = re.sub(
        r'(<circle\b[^>]*\bdata-station-id="(?:fastq_in|deliverables)"[^>]*\bstroke-width=")(?:2\.6|3\.4)(")',
        r'\g<1>3.8\2',
        text,
    )
    # attribute order variants (r before data-station-id)
    text = re.sub(
        r'(<circle\b[^>]*\br=")(?:6\.0|10\.0)("[^>]*\bdata-station-id="(?:fastq_in|deliverables)")',
        r'\g<1>12.0\2',
        text,
    )
    text = enlarge_terminal_file_badges(text)

    icons: list[str] = []
    station_specs = {
        "trim": ("Trim adapters", "AdapterRemoval", "up"),
        "host_id": ("Host identification", "FastQ Screen", "up"),
        "hmap": ("Map host", "BWA / Bowtie2", "down"),
        "mtdna": ("Map mtDNA", "BWA / Bowtie2", "up"),
        "hmerge": ("Merge samples", "samtools", "up"),
        "hdmg": ("Damage", "DamageProfiler", "up"),
        "hqual": ("Mapping QC", "Qualimap", "down"),
        "softclip": ("Soft clipping", "Python", "down"),
        "sexing": ("Genetic sexing", "sexing_residual", "down"),
        "prinseq": ("Dedup", "PRINSEQ++", "up"),
        "unaligned": ("Mammals reads removal", "Bowtie2", "down"),
        "pool": ("Pool samples", "cat · pigz", "up"),
        "kraken": ("Metagenome screening", "KrakenUniq", "up"),
        "malt": ("Metagenome screening", "MALT", "down"),
        "decom": ("Source tracking", "decOM", "down"),
        "detect": ("Pathogen detection", "E-value · MaltExtract", "up"),
        "pmap": ("Pathogen mapping", "BWA / Bowtie2", "up"),
        "auth": ("Authentication", "composite score", "down"),
    }
    for node_id, (purpose, tool, side) in station_specs.items():
        xy = node_xy(node_id)
        if xy:
            icons.append(_station_plate(xy[0], xy[1], purpose, tool, side=side))

    # Authentication details box under auth station (light SVG)
    auth_xy = node_xy("auth")
    if auth_xy:
        ax_, ay_ = auth_xy
        names = AUTH_METRIC_NAMES
        cols, chip_h, gap_x, gap_y = 4, 26.0, 10.0, 10.0
        pad_x, pad_y, title_h = 16.0, 14.0, 30.0
        chip_ws = [max(58.0, 8.2 * len(n) + 24) for n in names]
        nrows = (len(names) + cols - 1) // cols
        row_ws = []
        for r in range(nrows):
            chunk = chip_ws[r * cols : (r + 1) * cols]
            row_ws.append(sum(chunk) + gap_x * (len(chunk) - 1))
        inner_w = max(row_ws) if row_ws else 200
        box_w = inner_w + 2 * pad_x
        box_h = title_h + pad_y + nrows * chip_h + (nrows - 1) * gap_y + pad_y
        box_x, box_y = ax_ - box_w / 2, ay_ + 58
        accent = PAL["combo"]
        icons.extend([
            f'<rect class="pigsti-icon" x="{box_x}" y="{box_y}" width="{box_w}" height="{box_h}" '
            f'fill="#FFFFFF" stroke="{accent}" stroke-width="2.4"/>',
            f'<rect class="pigsti-icon" x="{box_x}" y="{box_y}" width="{box_w}" height="{title_h}" fill="{accent}"/>',
            f'<text class="pigsti-icon" x="{ax_}" y="{box_y + title_h / 2}" text-anchor="middle" '
            f'dominant-baseline="middle" font-size="12" font-weight="bold" fill="#FFFFFF" '
            f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">Authentication details</text>',
        ])
        y0 = box_y + title_h + pad_y
        for i, name in enumerate(names):
            r, c = divmod(i, cols)
            row_chunk = chip_ws[r * cols : (r + 1) * cols]
            row_w = sum(row_chunk) + gap_x * (len(row_chunk) - 1)
            x0 = ax_ - row_w / 2
            x = x0 + sum(row_chunk[:c]) + gap_x * c
            w = chip_ws[i]
            y = y0 + r * (chip_h + gap_y)
            icons.extend([
                f'<rect class="pigsti-icon" x="{x}" y="{y}" width="{w}" height="{chip_h}" rx="5" ry="5" '
                f'fill="#F8FAFC" stroke="#CBD5E1" stroke-width="1.2"/>',
                f'<text class="pigsti-icon" x="{x + w / 2}" y="{y + chip_h / 2}" text-anchor="middle" '
                f'dominant-baseline="middle" font-size="9.5" font-weight="bold" fill="#0F172A" '
                f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif">{escape(name)}</text>',
            ])

    animals = [(f"{name}.png", ICONS) for name in HOST_ANIMALS]
    hid = node_xy("host_id")
    if hid:
        cx, cy = hid
        animals_list = animals
        cols, rows = 5, 2
        icon_sz, gap, pad, title_h = 28, 5, 8, 14
        panel_w = pad * 2 + cols * icon_sz + (cols - 1) * gap
        panel_h = title_h + pad * 2 + rows * icon_sz + (rows - 1) * gap
        px = cx + 28
        py = cy + 14
        host_col = PAL["host"]
        icons.append(
            f'<rect class="pigsti-icon pigsti-host-card" x="{px}" y="{py}" '
            f'width="{panel_w}" height="{panel_h}" rx="10" ry="10" '
            f'fill="#FFFFFF" stroke="{host_col}" stroke-width="1.8"/>'
        )
        icons.append(
            f'<text class="pigsti-icon" x="{px + panel_w / 2}" y="{py + 12}" '
            f'text-anchor="middle" font-size="10" font-weight="bold" '
            f'fill="{host_col}" font-family="Helvetica Neue, Helvetica, Arial, sans-serif">'
            f"Host species</text>"
        )
        for i, (fname, folder) in enumerate(animals_list):
            uri = png_data_uri(folder / fname)
            if not uri:
                continue
            r, c = divmod(i, cols)
            ix = px + pad + c * (icon_sz + gap)
            iy = py + title_h + pad + r * (icon_sz + gap)
            icons.append(
                f'<image class="pigsti-icon" href="{uri}" x="{ix}" y="{iy}" '
                f'width="{icon_sz}" height="{icon_sz}" preserveAspectRatio="xMidYMid meet"/>'
            )

    # Section mascots (meta + pathogen; host covered by species card)
    placements = [
        ("metagenomics_screening", MICROBES / "E.coli.png", -54, 4, 50, 50),
        ("pathogen_route", MICROBES / "virus.png", -62, 2, 56, 56),
    ]
    for sec, path, dx, dy, w, h in placements:
        box = section_box(sec)
        uri = png_data_uri(path)
        if not box or not uri:
            continue
        x, y, bw, _bh = box
        icons.append(
            f'<image class="pigsti-icon" href="{uri}" x="{x + bw + dx}" y="{y + dy}" '
            f'width="{w}" height="{h}" preserveAspectRatio="xMidYMid meet"/>'
        )

    detect = node_xy("detect")
    if detect:
        cx, cy = detect
        for fname, ox in (("B.melitensis.png", -55), ("virus.png", 8)):
            uri = png_data_uri(MICROBES / fname)
            if not uri:
                continue
            icons.append(
                f'<image class="pigsti-icon" href="{uri}" x="{cx + ox}" y="{cy - 56}" '
                f'width="40" height="40" preserveAspectRatio="xMidYMid meet"/>'
            )

    # Route colour key for light SVG — bottom-right (no title)
    legend_items = [
        (PAL["prep"], "Preprocessing", False),
        (PAL["host"], "Host", False),
        (PAL["meta"], "Metagenomics", False),
        (PAL["detect"], "Detection", False),
        (PAL["path"], "Mapping", False),
        (PAL["combo"], "Authentication", False),
        (PAL["opt"], "Optional", True),
    ]
    row_h = 20
    cols = 2
    rows = (len(legend_items) + cols - 1) // cols
    col_w = 175
    pad_top, pad_bot = 16, 14
    leg_w = 24 + cols * col_w
    leg_h = pad_top + rows * row_h + pad_bot
    vb = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', text)
    svg_w = float(vb.group(1)) if vb else 1764.0
    svg_h = float(vb.group(2)) if vb else 626.0
    leg_x = svg_w - leg_w - 24
    leg_y = svg_h - leg_h - 18
    icons.extend([
        f'<rect class="pigsti-icon" x="{leg_x}" y="{leg_y}" width="{leg_w}" height="{leg_h}" fill="#FFFFFF"/>',
    ])
    for i, (col, label, dashed) in enumerate(legend_items):
        rr, cc = divmod(i, cols)
        lx = leg_x + 16 + cc * col_w
        ly = leg_y + pad_top + rr * row_h
        if dashed:
            icons.append(
                f'<line class="pigsti-icon" x1="{lx}" y1="{ly}" x2="{lx + 28}" y2="{ly}" '
                f'stroke="{col}" stroke-width="3.2" stroke-dasharray="1.8 1.3" stroke-linecap="butt"/>'
            )
        else:
            icons.append(
                f'<line class="pigsti-icon" x1="{lx}" y1="{ly}" x2="{lx + 28}" y2="{ly}" '
                f'stroke="{col}" stroke-width="7" stroke-linecap="butt"/>'
            )
        icons.append(
            f'<text class="pigsti-icon" x="{lx + 36}" y="{ly}" dominant-baseline="middle" '
            f'font-size="9.5" fill="#0F172A" font-family="DejaVu Sans, Helvetica, Arial, sans-serif">{escape(label)}</text>'
        )
    icons.append(
        f'<rect class="pigsti-icon" x="{leg_x}" y="{leg_y}" width="{leg_w}" height="{leg_h}" '
        f'fill="none" stroke="#000000" stroke-width="2.8"/>'
    )

    if not icons:
        # still write thickened strokes even without icons
        svg_path.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8")
        print(f"Thickened strokes in {svg_path.name} (no icons)")
        return

    block = "\n".join(icons) + "\n</svg>"
    if not text.rstrip().endswith("</svg>"):
        print(f"Unexpected SVG footer in {svg_path.name}")
        return
    svg_path.write_text(text.rstrip()[:-6] + block, encoding="utf-8")
    print(f"Beautified {svg_path.name} (thicker tracks + {len(icons)} icon elements)")


def render_nf_metro(outdir: Path, mmd: Path, stem: str) -> bool:
    if shutil.which("nf-metro") is None:
        return False
    svg_out = outdir / f"{stem}_light.svg"
    try:
        subprocess.run(
            [
                "nf-metro", "render", str(mmd),
                "-o", str(svg_out),
                "--format", "svg",
                "--theme", "light",
                "--mode", "light",
            ],
            check=True, capture_output=True, text=True,
        )
        inject_icons_into_nfmetro_svg(svg_out)
        print(f"Wrote {svg_out}")
        # High-quality PDF of the light SVG (icons preserved)
        try:
            from export_light_svg_pdf import svg_to_pdf
            pdf_out = svg_out.with_name(svg_out.stem + ".pdf")
            svg_to_pdf(svg_out, pdf_out, scale=2.5)
            print(f"Wrote {pdf_out}")
        except Exception as e:
            msg = str(e).encode("ascii", "replace").decode("ascii")
            print(f"light PDF export skipped: {msg}")
        return True
    except subprocess.CalledProcessError as e:
        err = (e.stderr or e.stdout or str(e))[:500]
        print(f"nf-metro failed: {err}")
        return False


def sync_canonical(outdir: Path, palette: str = "rainbow") -> None:
    """Mirror chosen palette to pigsti_subway_workflow.* (default: pride rainbow)."""
    pairs = [
        (f"pigsti_subway_{palette}.png", "pigsti_subway_workflow.png"),
        (f"pigsti_subway_{palette}.pdf", "pigsti_subway_workflow.pdf"),
        (f"pigsti_subway_{palette}.svg", "pigsti_subway_workflow.svg"),
        (f"pigsti_subway_{palette}_light.svg", "pigsti_subway_workflow_light.svg"),
        (f"pigsti_subway_{palette}_light.pdf", "pigsti_subway_workflow_light.pdf"),
        (f"pigsti_subway_{palette}.mmd", "pigsti_subway_workflow.mmd"),
    ]
    for src_name, dst_name in pairs:
        src, dst = outdir / src_name, outdir / dst_name
        if src.exists():
            dst.write_bytes(src.read_bytes())
            print(f"Synced {dst.name} from {src_name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--palette", default="rainbow", choices=list_palettes())
    ap.add_argument("--all-palettes", action="store_true")
    ap.add_argument("--outdir", type=Path, default=ROOT / "docs" / "images")
    args = ap.parse_args()

    names = list_palettes() if args.all_palettes else [args.palette]
    for name in names:
        set_palette(name)
        stem = f"pigsti_subway_{name}"
        mmd = write_mmd(args.outdir, stem)
        render_matplotlib(args.outdir, stem)
        render_nf_metro(args.outdir, mmd, stem)

    # Prefer rainbow as the paper canonical when available
    sync_name = "rainbow" if "rainbow" in names else (args.palette if not args.all_palettes else "hiroshige")
    if sync_name in names or (args.outdir / f"pigsti_subway_{sync_name}.png").exists():
        sync_canonical(args.outdir, sync_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
