#!/usr/bin/env python3
"""
PIGSTI workflow diagram — nf-core/rnaseq style (Function + Tool).

Grid layout: spine row (raw → trim → screen) · host row above · pathogen row below.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath

# ---------------------------------------------------------------------------
# Palette & grid
# ---------------------------------------------------------------------------
NF = {
    "brand": "#24B498",
    "brand_dark": "#21325B",
    "bg": "#FAFBFC",
    "lane_host": "#EEF4FC",
    "lane_path": "#FDF2F3",
    "lane_spine": "#F3F4F6",
    "box_fill": "#FFFFFF",
    "box_border": "#24B498",
    "tool_bar": "#F0F4F8",
    "tool_text": "#5C6B7A",
    "ink": "#21325B",
    "muted": "#8792A2",
    "join_border": "#475569",
    "checkpoint_fill": "#FFFBF0",
    "checkpoint_border": "#D97706",
    "optional_border": "#B0B8C4",
    "input_fill": "#21325B",
    "output_fill": "#24B498",
    "output_border": "#1A8A74",
    "line_host": "#2563EB",
    "line_path": "#DC2626",
    "line_shared": "#64748B",
    "bus": "#CBD5E1",
}

BOX_W = 2.28
BOX_H = 1.02
TOOL_H = 0.36
COL_GAP = 0.38
MARGIN_L = 0.55

SPINE_Y = 5.35
HOST_Y1 = 7.05
HOST_Y2 = 8.2
PATH_Y = 3.55
OPT_Y = 2.15
FINAL_COL = 12

PHASES: list[tuple[str, int, int, str]] = [
    ("Preprocess", 0, 3, NF["line_shared"]),
    ("Host · mtDNA", 4, 6, NF["line_host"]),
    ("Screening", 7, 8, NF["line_path"]),
    ("Mapping", 9, 11, NF["line_path"]),
    ("Deliverables", 12, 12, NF["brand"]),
]


def _cx(col: int) -> float:
    return MARGIN_L + col * (BOX_W + COL_GAP)


ROUTES = {
    "shared": NF["line_shared"],
    "host": NF["line_host"],
    "pathogen": NF["line_path"],
}


@dataclass(frozen=True)
class Step:
    pid: str
    col: int
    lane: str  # spine | host1 | host2 | path | opt
    function: str
    tool: str
    kind: str = "process"
    route: str = "shared"

    @property
    def x(self) -> float:
        return _cx(self.col)

    @property
    def y(self) -> float:
        return {
            "spine": SPINE_Y,
            "host1": HOST_Y1,
            "host2": HOST_Y2,
            "path": PATH_Y,
            "opt": OPT_Y,
        }[self.lane]

    @property
    def w(self) -> float:
        return BOX_W

    @property
    def h(self) -> float:
        return BOX_H


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    route: str = "shared"
    optional: bool = False
    label: str = ""
    manhattan: bool = False


def _build_steps() -> list[Step]:
    """Spine: AdapterRemoval → FastQ Screen (host); PRINSEQ++ feeds pathogen branch (matches Snakemake)."""
    return [
        Step("raw", 0, "spine", "Raw sequencing data", "FASTQ · samples.tsv", "input", "shared"),
        Step("adapt", 1, "spine", "Adapter & quality trimming", "AdapterRemoval · cutadapt", "process", "shared"),
        Step("fqs", 2, "spine", "Host & contaminant screen", "FastQ Screen", "process", "shared"),
        Step("prinseq", 3, "spine", "Complexity filter & dedup", "PRINSEQ++", "process", "shared"),
        Step("host_pcr", 4, "host1", "Host genome mapping", "BWA / Bowtie2", "process", "host"),
        Step("mtdna_pcr", 4, "host2", "Mitochondrial mapping", "BWA / Bowtie2", "process", "host"),
        Step("merge_h", 5, "host1", "Merge host libraries", "samtools", "join", "host"),
        Step("merge_m", 5, "host2", "Merge mtDNA libraries", "samtools", "join", "host"),
        Step("quali_h", 6, "host1", "Merged host QC", "Qualimap · DamageProfiler", "process", "host"),
        Step("quali_m", 6, "host2", "Merged mtDNA QC", "Qualimap", "process", "host"),
        Step("chimera", 5, "path", "Separate unaligned reads", "Bowtie2 · samtools", "process", "pathogen"),
        Step("merge_u", 6, "path", "Merge unaligned FASTQs", "cat · pigz", "join", "pathogen"),
        Step("kraken", 7, "path", "Taxonomic classification", "KrakenUniq 1.0.4", "process", "pathogen"),
        Step("decom", 7, "opt", "Source DNA screening", "decOM", "optional", "pathogen"),
        Step("evalue", 8, "path", "Kraken E-value scoring", "dExp_Escore.py · evalue/", "process", "pathogen"),
        Step("hops", 8, "opt", "aDNA pathogen screening", "HOPS", "optional", "pathogen"),
        Step("ckpt", 9, "path", "Select mapping targets", "cohort checkpoint", "checkpoint", "pathogen"),
        Step("pmap", 10, "path", "Pathogen reference mapping", "BWA / Bowtie2", "process", "pathogen"),
        Step("pqc", 11, "path", "Pathogen mapping QC", "ANI · DamageProfiler · Qualimap", "process", "pathogen"),
        Step("final", FINAL_COL, "spine", "Cohort reports & manifest", "Excel · PDF · TSV · catalog", "output", "shared"),
    ]


STEPS = _build_steps()
STEP_BY_ID = {s.pid: s for s in STEPS}

EDGES: list[Edge] = [
    Edge("raw", "adapt", label="FASTQ"),
    Edge("adapt", "fqs", label="trimmed collapsed"),
    Edge("adapt", "prinseq", label="trimmed collapsed"),
    Edge("prinseq", "chimera", route="pathogen", label="PRINSEQ-passed"),
    Edge("host_pcr", "merge_h", route="host"),
    Edge("mtdna_pcr", "merge_m", route="host"),
    Edge("merge_h", "quali_h", route="host"),
    Edge("merge_m", "quali_m", route="host"),
    Edge("chimera", "merge_u", route="pathogen"),
    Edge("merge_u", "kraken", route="pathogen"),
    Edge("kraken", "evalue", route="pathogen", label="report"),
    Edge("evalue", "ckpt", route="pathogen", label="pathogen.csv"),
    Edge("ckpt", "pmap", route="pathogen"),
    Edge("pmap", "pqc", route="pathogen"),
    Edge("merge_u", "decom", route="pathogen", optional=True, manhattan=True),
    Edge("hops", "ckpt", route="pathogen", optional=True, manhattan=True),
    Edge("quali_h", "final", route="host", manhattan=True),
    Edge("quali_m", "final", route="host", manhattan=True),
    Edge("pqc", "final", route="pathogen", manhattan=True),
]

X_MAX = _cx(FINAL_COL) + BOX_W + 0.8
Y_MAX = 9.5

MERMAID_TEMPLATE = """\
%%{{init: {{
  'theme': 'base',
  'themeVariables': {{
    'fontFamily': 'Helvetica Neue, Arial, sans-serif',
    'primaryColor': '#ffffff',
    'primaryTextColor': '#21325B',
    'primaryBorderColor': '#24B498',
    'lineColor': '#64748B',
    'clusterBkg': '#F6F7F8',
    'clusterBorder': '#21325B'
  }},
  'flowchart': {{ 'curve': 'basis', 'padding': 20, 'htmlLabels': true, 'rankSpacing': 50, 'nodeSpacing': 35 }}
}}}}%%
flowchart TB

  subgraph SPINE[" "]
    direction LR
    RAW([/"<b>Raw sequencing data</b><br/><sub>FASTQ · samples.tsv</sub>"/])
    ADAPT["<b>Adapter & quality trimming</b><br/><sub>AdapterRemoval · cutadapt</sub>"]
    FQS["<b>Host & contaminant screen</b><br/><sub>FastQ Screen · collapsed</sub>"]
    PR["<b>Filter & dedup</b><br/><sub>PRINSEQ++</sub>"]
    RAW --> ADAPT --> FQS
    ADAPT --> PR
  end

  subgraph HOST["Host & mtDNA"]
    direction LR
    HM["<b>Host genome mapping</b><br/><sub>BWA / Bowtie2 · collapsed</sub>"]
    DP["<b>DamageProfiler</b><br/><sub>per PCR BAM</sub>"]
    MT["<b>Mitochondrial mapping</b><br/><sub>BWA / Bowtie2</sub>"]
    MH{{"<b>Merge host</b><br/><sub>samtools · results/host/</sub>"}}
    SX["<b>Genetic sexing</b><br/><sub>R residual · Cow/Goat/Sheep/Dog <i>opt</i></sub>"]
    MM{{"<b>Merge mtDNA</b><br/><sub>samtools</sub>"}}
    QH["<b>Merged host QC</b><br/><sub>Qualimap · samples/</sub>"]
    QM["<b>Merged mtDNA QC</b><br/><sub>Qualimap · samples/</sub>"]
    HM --> MH --> QH
    HM -.-> DP
    MH -.-> SX
    MT --> MM --> QM
  end

  subgraph PATH["Pathogen & metagenomics"]
    direction LR
    CH["<b>Unaligned reads</b><br/><sub>Bowtie2 · samtools · PRINSEQ-passed</sub>"]
    MU{{"<b>Merge FASTQs</b><br/><sub>cat · pigz · pools/</sub>"}}
    KR["<b>Classification</b><br/><sub>KrakenUniq</sub>"]
    DC["<b>decOM</b><br/><sub>(config)</sub>"]
    HO["<b>HOPS</b><br/><sub>(config)</sub>"]
    EV["<b>Kraken E-value scoring</b><br/><sub>dExp_Escore · evalue/</sub>"]
    CK{{"<b>Checkpoint</b><br/><sub>pathogen_targets</sub>"}}
    PM["<b>Reference mapping</b><br/><sub>BWA / Bowtie2</sub>"]
    PQ["<b>Mapping QC</b><br/><sub>ANI · DamageProfiler · Qualimap</sub>"]
    CH --> MU --> KR
    MU -.-> DC
    KR -.-> HO
    KR --> EV --> CK --> PM --> PQ
    HO -.-> CK
  end

  FIN([/"<b>Cohort deliverables</b><br/><sub>final/ · manifests · catalogs</sub>"/])

  FQS --> HM
  FQS --> MT
  PR --> CH
  QH --> FIN
  QM --> FIN
  PQ --> FIN

  classDef input fill:#21325B,stroke:#21325B,color:#fff
  classDef process fill:#fff,stroke:#24B498,stroke-width:2px,color:#21325B
  classDef optional fill:#F8FAFC,stroke:#94A3B8,stroke-width:2px,stroke-dasharray:5 3,color:#21325B
  classDef custom fill:#F6F7F8,stroke:#475569,stroke-width:2px,color:#21325B
  classDef checkpoint fill:#FFFBF0,stroke:#D97706,stroke-width:2px,color:#21325B
  classDef output fill:#24B498,stroke:#1A8A74,color:#fff
  class RAW input
  class ADAPT,FQS,PR,HM,MT,DP,QH,QM,CH,KR,EV,PM,PQ process
  class MH,MM,MU,CK custom
  class DC,HO,SX optional
  class FIN output
"""


def write_mermaid(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MERMAID_TEMPLATE, encoding="utf-8")


def _center(s: Step) -> tuple[float, float]:
    return s.x + s.w / 2, s.y + s.h / 2


def _right_mid(s: Step) -> tuple[float, float]:
    return s.x + s.w, s.y + s.h / 2


def _left_mid(s: Step) -> tuple[float, float]:
    return s.x, s.y + s.h / 2


def _phase_x0_x1(col_start: int, col_end: int) -> tuple[float, float]:
    x0 = _cx(col_start) - 0.08
    x1 = _cx(col_end) + BOX_W + 0.08
    return x0, x1


def _draw_phase_headers(ax) -> None:
    y_top = 9.02
    for label, c0, c1, color in PHASES:
        x0, x1 = _phase_x0_x1(c0, c1)
        ax.add_patch(
            Rectangle(
                (x0, y_top - 0.22),
                x1 - x0,
                0.28,
                facecolor="#FFFFFF",
                edgecolor="#E2E8F0",
                linewidth=0.8,
                zorder=2,
            )
        )
        ax.text(
            (x0 + x1) / 2,
            y_top - 0.08,
            label,
            ha="center",
            va="center",
            fontsize=7,
            fontweight="bold",
            color=color,
            zorder=3,
        )


def _draw_lanes(ax) -> None:
    x0, x1 = MARGIN_L - 0.15, X_MAX - 0.35
    bands = [
        (HOST_Y1 - 0.45, HOST_Y2 + BOX_H + 0.35, NF["lane_host"], "Host & organellar genomics"),
        (SPINE_Y - 0.42, SPINE_Y + BOX_H + 0.42, NF["lane_spine"], None),
        (PATH_Y - 0.42, PATH_Y + BOX_H + 0.42, NF["lane_path"], None),
        (OPT_Y - 0.35, OPT_Y + BOX_H + 0.3, NF["lane_path"], None),
    ]
    for y0, y1, color, label in bands:
        ax.add_patch(
            Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color, edgecolor="none", zorder=0)
        )
        if label:
            ax.text(
                x0 + 0.12,
                y1 - 0.12,
                label,
                fontsize=8.5,
                fontweight="bold",
                color=NF["line_host"],
                va="top",
                zorder=1,
            )
    ax.text(
        x0 + 0.12,
        PATH_Y + BOX_H + 0.52,
        "Pathogen metagenomics & mapping",
        fontsize=8.5,
        fontweight="bold",
        color=NF["line_path"],
        va="bottom",
        zorder=1,
    )
    ax.text(
        x0 + 0.12,
        SPINE_Y + BOX_H + 0.48,
        "Shared preprocessing (all PCR libraries)",
        fontsize=8.5,
        fontweight="bold",
        color=NF["line_shared"],
        va="bottom",
        zorder=1,
    )


def _draw_title(ax) -> None:
    ax.text(
        MARGIN_L,
        9.15,
        "PIGSTI v3 workflow",
        fontsize=18,
        fontweight="bold",
        color=NF["brand_dark"],
        zorder=10,
    )
    ax.text(
        MARGIN_L,
        8.78,
        "nf-core/rnaseq style  ·  each box = Function (top) + Tool (bottom)",
        fontsize=9,
        color=NF["muted"],
        zorder=10,
    )


def _draw_step(ax, s: Step) -> None:
    x, y, w, h = s.x, s.y, s.w, s.h
    route_col = ROUTES.get(s.route, ROUTES["shared"])

    if s.kind == "input":
        face, border, fn_col, tool_col = NF["input_fill"], NF["input_fill"], "#FFFFFF", "#D1FAE8"
        ls, lw = "-", 1.8
        tool_bar = False
    elif s.kind == "output":
        face, border, fn_col, tool_col = NF["output_fill"], NF["output_border"], "#FFFFFF", "#D1FAE8"
        ls, lw = "-", 1.8
        tool_bar = False
    elif s.kind == "join":
        face, border, fn_col, tool_col = NF["box_fill"], NF["join_border"], NF["ink"], NF["tool_text"]
        ls, lw = "-", 1.6
        tool_bar = True
    elif s.kind == "checkpoint":
        face, border, fn_col, tool_col = NF["checkpoint_fill"], NF["checkpoint_border"], NF["ink"], NF["tool_text"]
        ls, lw = "-", 1.6
        tool_bar = True
    elif s.kind == "optional":
        face, border, fn_col, tool_col = NF["box_fill"], NF["optional_border"], NF["muted"], NF["tool_text"]
        ls, lw = "--", 1.2
        tool_bar = True
    else:
        face, border, fn_col, tool_col = NF["box_fill"], NF["box_border"], NF["ink"], NF["tool_text"]
        ls, lw = "-", 2.0
        tool_bar = True

    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.015,rounding_size=0.06",
            linewidth=lw,
            edgecolor=border,
            facecolor=face,
            linestyle=ls,
            zorder=4,
        )
    )

    if s.kind == "process":
        ax.add_patch(
            Rectangle((x, y + 0.04), 0.1, h - 0.08, facecolor=route_col, edgecolor="none", zorder=5)
        )

    if tool_bar and s.kind not in ("input", "output"):
        ax.add_patch(
            Rectangle((x + 0.06, y + 0.06), w - 0.12, TOOL_H, facecolor=NF["tool_bar"], edgecolor="none", zorder=5)
        )
        ax.plot([x + 0.1, x + w - 0.1], [y + 0.06 + TOOL_H, y + 0.06 + TOOL_H], color="#E2E8F0", lw=0.6, zorder=6)
        fn_y = y + 0.06 + TOOL_H + (h - 0.12 - TOOL_H) / 2
        tool_y = y + 0.06 + TOOL_H / 2
        ax.text(x + 0.14, tool_y, s.tool.replace("\n", " · "), ha="left", va="center", fontsize=6.3, color=tool_col, zorder=7)
    else:
        fn_y = y + h / 2 + (0.06 if s.tool else 0)
        tool_y = y + h / 2 - 0.18

    fs = 6.8 if len(s.function) > 26 else 7.4
    ax.text(
        x + w / 2 + (0.06 if s.kind == "process" else 0),
        fn_y,
        s.function,
        ha="center",
        va="center",
        fontsize=fs,
        fontweight="bold",
        color=fn_col,
        zorder=7,
    )
    if not tool_bar and s.tool and s.kind not in ("input", "output"):
        ax.text(
            x + w / 2,
            tool_y,
            s.tool,
            ha="center",
            va="center",
            fontsize=6.2,
            color=tool_col,
            zorder=7,
        )


def _manhattan_arrow(ax, p0: tuple[float, float], p1: tuple[float, float], color: str, dashed: bool) -> None:
    x0, y0 = p0
    x1, y1 = p1
    mid_x = (x0 + x1) / 2 if abs(x1 - x0) > 0.5 else x0 + 0.55
    verts = [(x0, y0), (mid_x, y0), (mid_x, y1), (x1, y1)]
    codes = [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO]
    path = MplPath(verts, codes)
    patch = PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        linewidth=1.6,
        linestyle=(0, (4, 3)) if dashed else "-",
        zorder=3,
    )
    ax.add_patch(patch)
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(mid_x, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6, shrinkA=0, shrinkB=2),
        zorder=3,
    )


def _draw_edge(ax, e: Edge) -> None:
    a, b = STEP_BY_ID[e.src], STEP_BY_ID[e.dst]
    color = NF["optional_border"] if e.optional else ROUTES.get(e.route, ROUTES["shared"])
    dashed = e.optional

    if e.manhattan:
        _manhattan_arrow(ax, _right_mid(a), _left_mid(b), color, dashed)
    else:
        p0, p1 = _right_mid(a), _left_mid(b)
        if abs(p0[1] - p1[1]) < 0.05:
            style = "arc3"
        else:
            style = "arc3,rad=0.08"
        ax.add_patch(
            FancyArrowPatch(
                p0,
                p1,
                arrowstyle="-|>",
                mutation_scale=10,
                linewidth=1.5,
                color=color,
                linestyle=(0, (4, 3)) if dashed else "-",
                connectionstyle=style,
                shrinkA=4,
                shrinkB=4,
                zorder=3,
            )
        )

    if e.label:
        p0, p1 = _right_mid(a), _left_mid(b)
        ax.text(
            (p0[0] + p1[0]) / 2,
            (p0[1] + p1[1]) / 2 + 0.18,
            e.label,
            ha="center",
            fontsize=6,
            color=NF["muted"],
            bbox=dict(boxstyle="round,pad=0.08", facecolor="#fff", edgecolor="#E2E8F0", linewidth=0.5),
            zorder=8,
        )


def _draw_fork_bus(ax) -> None:
    """Branch after FastQ Screen: host + mtDNA (pathogen branch from PRINSEQ)."""
    fqs = STEP_BY_ID["fqs"]
    jx = _cx(3) - 0.22
    cy = SPINE_Y + BOX_H / 2
    rx = fqs.x + fqs.w

    ax.plot([rx, jx], [cy, cy], color=ROUTES["shared"], lw=2, solid_capstyle="round", zorder=3)
    ax.plot([jx, jx], [HOST_Y1 + BOX_H / 2, HOST_Y2 + BOX_H / 2], color=NF["bus"], lw=2, zorder=2)

    for step_id, y_tgt in [
        ("host_pcr", HOST_Y1 + BOX_H / 2),
        ("mtdna_pcr", HOST_Y2 + BOX_H / 2),
    ]:
        tgt = STEP_BY_ID[step_id]
        tx = tgt.x
        col = ROUTES["host"]
        ax.plot([jx, tx], [y_tgt, y_tgt], color=col, lw=1.8, zorder=3)
        ax.annotate(
            "",
            xy=(tx, y_tgt),
            xytext=(jx, y_tgt),
            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8, shrinkA=0, shrinkB=3),
            zorder=4,
        )


def _draw_legend(ax) -> None:
    lx, ly = MARGIN_L, 0.35
    ax.add_patch(
        FancyBboxPatch(
            (lx, ly),
            4.6,
            1.15,
            boxstyle="round,pad=0.03,rounding_size=0.08",
            facecolor="#fff",
            edgecolor="#E2E8F0",
            linewidth=1,
            zorder=5,
        )
    )
    ax.text(lx + 0.15, ly + 0.88, "Box layout", fontsize=8, fontweight="bold", color=NF["ink"])
    ex = lx + 0.15
    ey = ly + 0.15
    ew, eh = 1.5, 0.65
    ax.add_patch(
        FancyBboxPatch(
            (ex, ey),
            ew,
            eh,
            boxstyle="round,pad=0.02,rounding_size=0.04",
            linewidth=2,
            edgecolor=NF["box_border"],
            facecolor="#fff",
            zorder=6,
        )
    )
    ax.add_patch(
        Rectangle((ex + 0.05, ey + 0.05), ew - 0.1, 0.22, facecolor=NF["tool_bar"], edgecolor="none", zorder=7)
    )
    ax.text(ex + ew / 2, ey + eh * 0.68, "Function", ha="center", fontsize=6.5, fontweight="bold", color=NF["ink"])
    ax.text(ex + ew / 2, ey + 0.16, "Tool", ha="center", fontsize=6, color=NF["tool_text"])
    ax.text(lx + 1.85, ly + 0.55, "→ workflow direction", fontsize=7, color=NF["muted"])
    for i, (col, lab) in enumerate(
        [(NF["line_shared"], "Shared"), (NF["line_host"], "Host"), (NF["line_path"], "Pathogen")]
    ):
        ax.plot([lx + 1.85, lx + 2.15], [ly + 0.28 - i * 0.18, ly + 0.28 - i * 0.18], color=col, lw=3)
        ax.text(lx + 2.22, ly + 0.28 - i * 0.18, lab, fontsize=6.5, va="center", color=NF["muted"])


def render(out_png: Path, out_svg: Path, mmd_path: Path | None = None) -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

    fig_h = 10
    fig, ax = plt.subplots(figsize=(X_MAX * 0.92, fig_h))
    ax.set_xlim(0, X_MAX)
    ax.set_ylim(0, Y_MAX)
    ax.axis("off")
    ax.set_facecolor(NF["bg"])

    _draw_title(ax)
    _draw_phase_headers(ax)
    _draw_lanes(ax)
    _draw_fork_bus(ax)
    for s in STEPS:
        _draw_step(ax, s)
    for e in EDGES:
        _draw_edge(ax, e)
    _draw_legend(ax)

    ax.annotate(
        "",
        xy=(X_MAX - 0.3, 0.55),
        xytext=(MARGIN_L, 0.55),
        arrowprops=dict(arrowstyle="-|>", color="#E2E8F0", lw=1.0),
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=260, bbox_inches="tight", facecolor=NF["bg"], pad_inches=0.35)
    fig.savefig(out_svg, bbox_inches="tight", facecolor=NF["bg"], pad_inches=0.35)
    plt.close(fig)
    if mmd_path:
        write_mermaid(mmd_path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render PIGSTI workflow diagram")
    parser.add_argument("--png", default="docs/images/pipeline_workflow.png")
    parser.add_argument("--svg", default="docs/images/pipeline_workflow.svg")
    parser.add_argument("--mmd", default="docs/images/pipeline_workflow.mmd")
    args = parser.parse_args()
    render(Path(args.png), Path(args.svg), Path(args.mmd))
    print(f"Wrote {args.png}")
    print(f"Wrote {args.svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
