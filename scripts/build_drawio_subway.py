#!/usr/bin/env python3
"""Subway / metro-map style PIGSTI workflow as editable draw.io.

Matches nf-metro **light** aesthetic:
  soft gray section panels, dark station markers, coloured tracks only.
"""

from __future__ import annotations

import argparse
import base64
import html
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from workflow_palette import AUTH_METRIC_NAMES, METRICS, PAL, list_palettes, set_palette  # noqa: E402

ICONS = ROOT / "docs" / "images" / "Icons"
MICROBES = ROOT / "docs" / "images" / "Icons microbes"
OUT_DIR = ROOT / "docs" / "drawio"
HEATMAP = ROOT / "docs" / "pathogen_detection_scores_heatmap.png"
if not HEATMAP.exists():
    HEATMAP = ROOT / "final" / "pathogen_detection_scores_heatmap.png"

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
}

PW, PH = 1680, 760
TRACK_W = 11.0      # bold light-metro tracks
TRACK_W_OPT = 5.0   # thinner optional so dash gaps stay visible
DASH_OPT = "10 9"


def _id() -> str:
    return uuid.uuid4().hex[:12]


def data_uri(path: Path) -> str:
    raw = path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    mime = "image/png" if path.suffix.lower() == ".png" else "image/svg+xml"
    return f"data:{mime};base64,{b64}"


def resolve_icon(name: str, folder: Path) -> Path:
    png = folder / "png" / f"{name}.png"
    if png.exists():
        return png
    for p in (folder / "png").glob("*.png") if (folder / "png").exists() else []:
        if p.stem.lower() == name.lower():
            return p
    svg = folder / f"{name}.svg"
    if svg.exists():
        return svg
    for p in folder.glob("*.svg"):
        if p.stem.lower() == name.lower():
            return p
    raise FileNotFoundError(name)


class SubwayDrawio:
    def __init__(self) -> None:
        self.cells: list[str] = []
        self._n = 2
        self.icon = {
            "cattle": data_uri(resolve_icon("Cattle", ICONS)),
            "sheep": data_uri(resolve_icon("Sheep", ICONS)),
            "goat": data_uri(resolve_icon("Goat", ICONS)),
            "pig": data_uri(resolve_icon("pig", ICONS)),
            "dog": data_uri(resolve_icon("Dog", ICONS)),
            "horse": data_uri(resolve_icon("Horse", ICONS)),
            "cat": data_uri(resolve_icon("Cat", ICONS)),
            "camel": data_uri(resolve_icon("Camel", ICONS)),
            "deer": data_uri(resolve_icon("Deer", ICONS)),
            "rat": data_uri(resolve_icon("Rat", ICONS)),
            "ecoli": data_uri(resolve_icon("E.coli", MICROBES)),
            "virus": data_uri(resolve_icon("virus", MICROBES)),
            "brucella": data_uri(resolve_icon("B.melitensis", MICROBES)),
            "heatmap": data_uri(HEATMAP) if HEATMAP.exists() else "",
        }

    def next_id(self) -> str:
        i = str(self._n)
        self._n += 1
        return i

    def image(self, x: float, y: float, w: float, h: float, uri: str) -> str:
        cid = self.next_id()
        esc = uri.replace("&", "&amp;")
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="shape=image;html=1;verticalAlign=top;'
            f'verticalLabelPosition=bottom;imageAspect=1;aspect=fixed;image={esc}" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        return cid

    def panel(self, x: float, y: float, w: float, h: float, title: str, n: int) -> None:
        cid = self.next_id()
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="rounded=1;whiteSpace=wrap;html=1;'
            f'fillColor={LIGHT["section_fill"]};strokeColor={LIGHT["section_stroke"]};'
            f'strokeWidth=1;arcSize=8;" vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        badge = self.next_id()
        self.cells.append(
            f'<mxCell id="{badge}" value="{n}" style="ellipse;whiteSpace=wrap;html=1;aspect=fixed;'
            f'fillColor={LIGHT["num"]};strokeColor=none;fontColor=#FFFFFF;fontStyle=1;fontSize=10;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x + 10}" y="{y - 12}" width="20" height="20" as="geometry"/></mxCell>'
        )
        lid = self.next_id()
        self.cells.append(
            f'<mxCell id="{lid}" value="&lt;b&gt;{html.escape(title)}&lt;/b&gt;" '
            f'style="text;html=1;strokeColor=none;fillColor=none;align=left;fontSize=13;'
            f'fontColor={LIGHT["section_label"]};fontFamily=Helvetica;" vertex="1" parent="1">'
            f'<mxGeometry x="{x + 36}" y="{y - 14}" width="{min(w, 280)}" height="20" as="geometry"/></mxCell>'
        )

    def station(
        self,
        x: float,
        y: float,
        purpose: str,
        tool: str = "",
        *,
        kind: str = "node",
        label_side: str = "up",
        r: float = 9,
    ) -> str:
        """nf-metro light station: solid circle marker."""
        cid = self.next_id()
        mk = LIGHT["marker"]
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="ellipse;whiteSpace=wrap;html=1;aspect=fixed;'
            f'fillColor=#FFFFFF;strokeColor={mk};strokeWidth=2.6;" vertex="1" parent="1">'
            f'<mxGeometry x="{x - r}" y="{y - r}" width="{2 * r}" height="{2 * r}" as="geometry"/></mxCell>'
        )

        # Purpose (bold sans process) + tool chip (mono in indigo pill)
        val = (
            f'&lt;font face=&quot;Helvetica&quot; color=&quot;#0F172A&quot;&gt;'
            f'&lt;b&gt;{html.escape(purpose)}&lt;/b&gt;&lt;/font&gt;'
        )
        tid = self.next_id()
        if label_side == "up":
            ly = y - r - 42
        else:
            ly = y + r + 4
        self.cells.append(
            f'<mxCell id="{tid}" value="{val}" style="text;html=1;align=center;verticalAlign=middle;'
            f'strokeColor=none;fillColor=none;fontSize=12;fontFamily=Helvetica;fontStyle=1;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x - 70}" y="{ly}" width="140" height="16" as="geometry"/></mxCell>'
        )
        if tool:
            py = ly + 16
            pill_w = max(48, 7 * len(tool) + 16)
            self.cells.append(
                f'<mxCell id="{self.next_id()}" value="" '
            f'style="rounded=1;whiteSpace=wrap;html=1;fillColor=#4338CA;'
            f'strokeColor=#312E81;strokeWidth=1.2;arcSize=40;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x - pill_w / 2}" y="{py}" width="{pill_w}" height="15" '
            f'as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" '
            f'value="&lt;font face=&quot;Courier New&quot; color=&quot;#FFFFFF&quot; '
            f'style=&quot;font-size:8px;&quot;&gt;&lt;b&gt;{html.escape(tool)}&lt;/b&gt;&lt;/font&gt;" '
                f'style="text;html=1;align=center;verticalAlign=middle;strokeColor=none;fillColor=none;'
                f'fontFamily=Courier New;fontSize=8;" vertex="1" parent="1">'
                f'<mxGeometry x="{x - pill_w / 2}" y="{py}" width="{pill_w}" height="14" '
                f'as="geometry"/></mxCell>'
            )
        return cid

    def track(
        self,
        points: list[tuple[float, float]],
        color: str,
        *,
        dashed: bool = False,
        width: float | None = None,
    ) -> None:
        if len(points) < 2:
            return
        eid = self.next_id()
        lw = width if width is not None else (TRACK_W_OPT if dashed else TRACK_W)
        dash = f"dashed=1;dashPattern={DASH_OPT};" if dashed else ""
        x0, y0 = points[0]
        x1, y1 = points[-1]
        extras = ""
        if len(points) > 2:
            mid = "".join(f'<mxPoint x="{x}" y="{y}"/>' for x, y in points[1:-1])
            extras = f'<Array as="points">{mid}</Array>'
        # butt-like appearance for dashed: no rounded ends (draw.io rounded=0)
        round_flag = "0" if dashed else "1"
        self.cells.append(
            f'<mxCell id="{eid}" value="" style="endArrow=none;html=1;rounded={round_flag};curved=0;'
            f'strokeColor={color};strokeWidth={lw};{dash}" edge="1" parent="1">'
            f'<mxGeometry relative="1" as="geometry">'
            f'<mxPoint x="{x0}" y="{y0}" as="sourcePoint"/>'
            f'<mxPoint x="{x1}" y="{y1}" as="targetPoint"/>'
            f"{extras}"
            f"</mxGeometry></mxCell>"
        )

    def auth_note(self, x: float, y: float, w: float, h: float) -> None:
        """Authentication details card with short metric names only."""
        note = self.next_id()
        line1 = html.escape(" · ".join(AUTH_METRIC_NAMES[:4]))
        line2 = html.escape(" · ".join(AUTH_METRIC_NAMES[4:]))
        self.cells.append(
            f'<mxCell id="{note}" value="&lt;font style=&quot;font-size:8px;&quot; '
            f'color=&quot;#0F172A&quot;&gt;&lt;b&gt;Authentication details&lt;/b&gt;&lt;br&gt;'
            f'{line1}&lt;br&gt;{line2}&lt;/font&gt;" '
            f'style="rounded=1;whiteSpace=wrap;html=1;align=center;verticalAlign=middle;'
            f'fillColor=#FFFFFF;strokeColor={PAL["combo"]};strokeWidth=2;arcSize=8;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h + 8}" as="geometry"/></mxCell>'
        )

    def host_species_card(self, x: float, y: float) -> None:
        """Host icons — all animal silhouettes in a 5×2 card."""
        keys = [
            "cattle", "sheep", "goat", "pig", "dog",
            "horse", "cat", "camel", "deer", "rat",
        ]
        cols, rows = 5, 2
        icon, gap, pad, title_h = 26, 5, 8, 14
        w = pad * 2 + cols * icon + (cols - 1) * gap
        h = title_h + pad * 2 + rows * icon + (rows - 1) * gap
        cid = self.next_id()
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="rounded=1;whiteSpace=wrap;html=1;'
            f'fillColor=#FFFFFF;strokeColor={PAL["host"]};strokeWidth=1.8;arcSize=12;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        tid = self.next_id()
        self.cells.append(
            f'<mxCell id="{tid}" value="&lt;font color=&quot;{PAL["host"]}&quot;&gt;'
            f'&lt;b&gt;Host species&lt;/b&gt;&lt;/font&gt;" '
            f'style="text;html=1;align=center;fontSize=9;fontFamily=Helvetica;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y + 2}" width="{w}" height="12" as="geometry"/></mxCell>'
        )
        for i, key in enumerate(keys):
            r, c = divmod(i, cols)
            self.image(
                x + pad + c * (icon + gap),
                y + title_h + pad + r * (icon + gap) - 2,
                icon, icon, self.icon[key],
            )
    def build(self) -> str:
        # Soft gray section panels (nf-metro light)
        self.panel(24, 70, 280, 90, "Preprocessing", 1)
        self.panel(24, 190, 560, 165, "Host route", 2)
        self.panel(24, 385, 620, 160, "Metagenomics screening", 3)
        self.panel(660, 70, 980, 250, "Pathogen route", 4)

        # Section mascots
        self.image(580, 385, 46, 46, self.icon["ecoli"])
        self.image(1370, 70, 56, 56, self.icon["virus"])

        Yp, Yh, Yh2, Ym, Yo, Yd, Yout = 115, 250, 310, 445, 505, 175, 280
        X = {
            "raw": 70, "trim": 200,
            "hid": 80, "hmap": 230, "mtdna": 230, "hmerge": 360, "hqual": 490, "hdmg": 490,
            "prinseq": 80, "unal": 200, "pool": 320, "kraken": 460, "malt": 460,
            "detect": 720, "pmap": 900, "auth": 1080, "final": 1280,
        }

        # Tracks (colour carries the route identity)
        self.track([(X["raw"], Yp), (X["trim"], Yp)], PAL["prep"])
        self.track([(X["trim"], Yp), (X["trim"], Yh), (X["hid"], Yh)], PAL["host"])
        self.track([(X["hid"], Yh), (X["hmap"], Yh), (X["hmerge"], Yh), (X["hqual"], Yh)], PAL["host"])
        self.track([(X["hid"], Yh2), (X["mtdna"], Yh2), (X["hmerge"], Yh2), (X["hmerge"], Yh)], PAL["host"])
        self.track([(X["hmerge"], Yh), (X["hdmg"], Yh2)], PAL["host"])
        self.track([(X["trim"], Yp), (X["trim"], Ym), (X["prinseq"], Ym)], PAL["meta"])
        self.track([(X["prinseq"], Ym), (X["unal"], Ym), (X["pool"], Ym), (X["kraken"], Ym)], PAL["meta"])
        self.track([(X["pool"], Ym), (X["malt"], Yo)], PAL["opt"], dashed=True)
        self.track([(X["kraken"], Ym), (X["detect"], Yd)], PAL["detect"])
        self.track([(X["detect"], Yd), (X["pmap"], Yd)], PAL["path"])
        self.track([(X["pmap"], Yd), (X["auth"], Yd), (X["auth"], Yout), (X["final"], Yout)], PAL["combo"])
        self.track([(X["malt"], Yo), (X["detect"] - 40, Yo), (X["detect"], Yd)], PAL["opt"], dashed=True)
        self.track([(X["hqual"], Yh), (X["hqual"] + 40, Yh), (X["final"] - 80, Yout), (X["final"], Yout)], PAL["host"])
        self.track(
            [(X["hdmg"], Yh2), (X["hdmg"] + 80, Yh2), (X["final"] - 40, Yout + 40), (X["final"], Yout)],
            PAL["host"], dashed=True, width=TRACK_W_OPT,
        )

        # Stations (dark markers)
        self.station(X["raw"], Yp, "Raw reads", "FASTQ", kind="start", label_side="up")
        self.station(X["trim"], Yp, "Trim adapters", "AdapterRemoval", label_side="down")
        self.station(X["hid"], Yh, "Host identification", "FastQ Screen", label_side="up")
        self.host_species_card(X["hid"] + 50, Yh + 18)
        self.station(X["hmap"], Yh, "Map host", "BWA / Bowtie2", label_side="up")
        self.station(X["mtdna"], Yh2, "Map mtDNA", "BWA / Bowtie2", label_side="down")
        self.station(X["hmerge"], Yh, "Merge samples", "samtools", kind="join", label_side="up")
        self.station(X["hqual"], Yh, "Mapping QC", "Qualimap", label_side="up")
        self.station(X["hdmg"], Yh2, "Damage", "DamageProfiler", label_side="down")
        self.station(X["prinseq"], Ym, "Dedup", "PRINSEQ++", label_side="up")
        self.station(X["unal"], Ym, "Mammals removal", "Bowtie2", label_side="up")
        self.station(X["pool"], Ym, "Pool samples", "cat · pigz", kind="join", label_side="up")
        self.station(X["kraken"], Ym, "Metagenome screening", "KrakenUniq", label_side="up")
        self.station(X["malt"], Yo, "Metagenome screening", "MALT", kind="opt", label_side="down")
        self.image(X["detect"] - 55, Yd - 62, 36, 36, self.icon["brucella"])
        self.image(X["detect"] + 5, Yd - 62, 36, 36, self.icon["virus"])
        self.station(X["detect"], Yd, "Pathogen detection", "E-value · MaltExtract", label_side="up")
        self.station(X["pmap"], Yd, "Pathogen mapping", "BWA / Bowtie2", label_side="up")
        self.station(X["auth"], Yd, "Authentication", "composite score", label_side="down")
        self.auth_note(X["auth"] - 110, Yd + 58, 220, 36)
        self.station(X["final"], Yout, "Deliverables", "Excel · heatmap · PDF", kind="end", label_side="up")

        # Detection-score heatmap inset (compact)
        if self.icon.get("heatmap"):
            hx, hy, hw, hh = X["final"] + 60, Yout - 55, 140, 120
            self.cells.append(
                f'<mxCell id="{self.next_id()}" value="" '
                f'style="rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFFFF;'
                f'strokeColor={LIGHT["section_stroke"]};strokeWidth=1;arcSize=10;" '
                f'vertex="1" parent="1">'
                f'<mxGeometry x="{hx - 6}" y="{hy - 16}" width="{hw + 12}" height="{hh + 28}" '
                f'as="geometry"/></mxCell>'
            )
            self.cells.append(
                f'<mxCell id="{self.next_id()}" '
                f'value="&lt;font color=&quot;{LIGHT["section_label"]}&quot;&gt;'
                f'&lt;b&gt;Detection score heatmap&lt;/b&gt;&lt;/font&gt;" '
                f'style="text;html=1;align=center;fontSize=7;fontFamily=Helvetica;" '
                f'vertex="1" parent="1">'
                f'<mxGeometry x="{hx - 6}" y="{hy - 14}" width="{hw + 12}" height="12" '
                f'as="geometry"/></mxCell>'
            )
            self.image(hx, hy, hw, hh, self.icon["heatmap"])

        # Framed legend with title bar
        items = [
            (PAL["prep"], "Preprocessing", False),
            (PAL["host"], "Host", False),
            (PAL["meta"], "Metagenomics", False),
            (PAL["detect"], "Detection", False),
            (PAL["path"], "Mapping", False),
            (PAL["combo"], "Authentication", False),
            (PAL["opt"], "Optional", True),
        ]
        slot = 118
        box_x, box_y = 28, PH - 78
        box_w = 24 + len(items) * slot
        box_h = 58
        title_h = 20
        self.cells.append(
            f'<mxCell id="{self.next_id()}" value="" '
            f'style="rounded=0;whiteSpace=wrap;html=1;fillColor=#FFFFFF;'
            f'strokeColor=#1F2937;strokeWidth=2.2;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{box_x}" y="{box_y}" width="{box_w}" height="{box_h}" '
            f'as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" value="" '
            f'style="rounded=0;whiteSpace=wrap;html=1;fillColor=#F3F4F6;'
            f'strokeColor=#1F2937;strokeWidth=2.2;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{box_x}" y="{box_y}" width="{box_w}" height="{title_h}" '
            f'as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" value="&lt;b&gt;Legend&lt;/b&gt;" '
            f'style="text;html=1;fontSize=11;align=center;fontColor=#111827;'
            f'fontFamily=Helvetica;" vertex="1" parent="1">'
            f'<mxGeometry x="{box_x}" y="{box_y + 2}" width="{box_w}" height="16" '
            f'as="geometry"/></mxCell>'
        )
        lx, leg_y = box_x + 14, box_y + title_h + 22
        for col, label, dashed in items:
            lid = self.next_id()
            dash = "dashed=1;dashPattern=10 9;" if dashed else ""
            self.cells.append(
                f'<mxCell id="{lid}" value="" style="endArrow=none;html=1;strokeColor={col};'
                f'strokeWidth={"4" if dashed else "8"};{dash}" edge="1" parent="1">'
                f'<mxGeometry relative="1" as="geometry">'
                f'<mxPoint x="{lx}" y="{leg_y}" as="sourcePoint"/>'
                f'<mxPoint x="{lx + 22}" y="{leg_y}" as="targetPoint"/>'
                f"</mxGeometry></mxCell>"
            )
            tid = self.next_id()
            self.cells.append(
                f'<mxCell id="{tid}" value="{html.escape(label)}" '
                f'style="text;html=1;fontSize=9;fontColor=#111827;align=left;'
                f'fontFamily=Helvetica;" vertex="1" parent="1">'
                f'<mxGeometry x="{lx + 26}" y="{leg_y - 8}" width="95" height="16" as="geometry"/></mxCell>'
            )
            lx += slot

        title = self.next_id()
        self.cells.insert(
            0,
            f'<mxCell id="{title}" value="&lt;b&gt;PIGSTI&lt;/b&gt;" '
            f'style="text;html=1;fontSize=22;align=left;fontColor={LIGHT["ink"]};'
            f'fontFamily=Helvetica;" vertex="1" parent="1">'
            f'<mxGeometry x="24" y="16" width="200" height="28" as="geometry"/></mxCell>',
        )
        # Process vs Tool key with chip example
        self.cells.append(
            f'<mxCell id="{self.next_id()}" value="" '
            f'style="rounded=0;whiteSpace=wrap;html=1;fillColor=#FFFFFF;'
            f'strokeColor=#1F2937;strokeWidth=1.6;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="1385" y="10" width="250" height="48" as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" '
            f'value="&lt;font color=&quot;#6B7280&quot; style=&quot;font-size:8px;&quot;&gt;Process&lt;/font&gt; '
            f'&lt;font face=&quot;Helvetica&quot;&gt;&lt;b&gt;Classify&lt;/b&gt;&lt;/font&gt;&lt;br&gt;'
            f'&lt;font color=&quot;#6B7280&quot; style=&quot;font-size:8px;&quot;&gt;Tool&lt;/font&gt;" '
            f'style="text;html=1;align=left;verticalAlign=middle;strokeColor=none;fillColor=none;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="1395" y="12" width="120" height="42" as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" value="" '
            f'style="rounded=1;whiteSpace=wrap;html=1;fillColor=#4338CA;'
            f'strokeColor=#312E81;strokeWidth=1.2;arcSize=40;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="1485" y="30" width="88" height="15" as="geometry"/></mxCell>'
        )
        self.cells.append(
            f'<mxCell id="{self.next_id()}" '
            f'value="&lt;font face=&quot;Courier New&quot; color=&quot;#FFFFFF&quot; '
            f'style=&quot;font-size:8px;&quot;&gt;&lt;b&gt;KrakenUniq&lt;/b&gt;&lt;/font&gt;" '
            f'style="text;html=1;align=center;verticalAlign=middle;strokeColor=none;fillColor=none;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="1485" y="30" width="88" height="15" as="geometry"/></mxCell>'
        )

        body = "\n        ".join(self.cells)
        return f"""<mxfile host="app.diagrams.net" agent="PIGSTI" version="22.1.0" type="device">
  <diagram name="PIGSTI subway light" id="{_id()}">
    <mxGraphModel dx="{PW}" dy="{PH}" grid="0" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1"
      fold="1" page="1" pageScale="1" pageWidth="{PW}" pageHeight="{PH}" math="0" shadow="0"
      background="{LIGHT["bg"]}">
      <root>
        <mxCell id="0"/>
        <mxCell id="1" parent="0"/>
        {body}
      </root>
    </mxGraphModel>
  </diagram>
</mxfile>
"""


def main() -> int:
    ap = argparse.ArgumentParser(description="Build nf-metro-light style subway draw.io")
    ap.add_argument("--palette", default="hiroshige", choices=list_palettes())
    ap.add_argument("--all-palettes", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    names = list_palettes() if args.all_palettes else [args.palette]
    for name in names:
        set_palette(name)
        if args.all_palettes:
            out = OUT_DIR / f"pigsti_subway_{name}.drawio"
        elif name == "hiroshige":
            out = OUT_DIR / "pigsti_subway_workflow.drawio"
        else:
            out = OUT_DIR / f"pigsti_subway_{name}.drawio"
        out.write_text(SubwayDrawio().build(), encoding="utf-8")
        print(f"Wrote {out} (palette={name}, style=light)")
        if name == "hiroshige" and args.all_palettes:
            canon = OUT_DIR / "pigsti_subway_workflow.drawio"
            canon.write_text(out.read_text(encoding="utf-8"), encoding="utf-8")
            print(f"Wrote {canon}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
