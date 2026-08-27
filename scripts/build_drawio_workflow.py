#!/usr/bin/env python3
"""Publication PIGSTI workflow figure (draw.io). MetBrewer palette + icon assets."""

from __future__ import annotations

import base64
import html
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from workflow_palette import METRICS, PAL  # noqa: E402

OUT = ROOT / "docs" / "drawio" / "pigsti_workflow.drawio"
ICONS = ROOT / "docs" / "images" / "Icons"
MICROBES = ROOT / "docs" / "images" / "Icons microbes"

PW, PH = 1280, 460


def _id() -> str:
    return uuid.uuid4().hex[:12]


def icon_data_uri(name: str, folder: Path) -> str:
    """Prefer PNG (visible in more export paths) then SVG."""
    png = folder / "png" / f"{name}.png"
    path = png if png.exists() else folder / f"{name}.svg"
    if not path.exists():
        for p in (folder / "png").glob("*.png") if (folder / "png").exists() else []:
            if p.stem.lower() == name.lower():
                path = p
                break
        else:
            path = next(folder.glob(f"{name}.*"), folder / f"{name}.svg")
    raw = path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    mime = "image/png" if path.suffix.lower() == ".png" else "image/svg+xml"
    return f"data:{mime};base64,{b64}"


class DrawioBuilder:
    def __init__(self) -> None:
        self.cells: list[str] = []
        self._n = 2
        self.icon = {
            "cattle": icon_data_uri("Cattle", ICONS),
            "sheep": icon_data_uri("Sheep", ICONS),
            "goat": icon_data_uri("Goat", ICONS),
            "pig": icon_data_uri("pig", ICONS),
            "dog": icon_data_uri("Dog", ICONS),
            "horse": icon_data_uri("Horse", ICONS),
            "ecoli": icon_data_uri("E.coli", MICROBES),
            "virus": icon_data_uri("virus", MICROBES),
            "brucella": icon_data_uri("B.melitensis", MICROBES),
        }

    def next_id(self) -> str:
        i = str(self._n)
        self._n += 1
        return i

    def image(self, x: float, y: float, w: float, h: float, data_uri: str) -> str:
        cid = self.next_id()
        # Escape for XML attribute
        uri = data_uri.replace("&", "&amp;")
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="shape=image;html=1;verticalAlign=top;'
            f'verticalLabelPosition=bottom;imageAspect=1;aspect=fixed;image={uri}" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        return cid

    def section(self, x: float, y: float, w: float, h: float, title: str, n: int, *, fill: str) -> None:
        cid = self.next_id()
        self.cells.append(
            f'<mxCell id="{cid}" value="" style="rounded=1;whiteSpace=wrap;html=1;'
            f'fillColor={fill};strokeColor={PAL["section_stroke"]};strokeWidth=1;arcSize=5;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        badge = self.next_id()
        self.cells.append(
            f'<mxCell id="{badge}" value="{n}" style="ellipse;whiteSpace=wrap;html=1;aspect=fixed;'
            f'fillColor={PAL["prep"]};strokeColor=none;fontColor=#FFFFFF;fontStyle=1;fontSize=9;" '
            f'vertex="1" parent="1">'
            f'<mxGeometry x="{x + 5}" y="{y - 14}" width="18" height="18" as="geometry"/></mxCell>'
        )
        lid = self.next_id()
        self.cells.append(
            f'<mxCell id="{lid}" value="&lt;b&gt;{html.escape(title)}&lt;/b&gt;" '
            f'style="text;html=1;strokeColor=none;fillColor=none;align=left;fontSize=11;'
            f'fontColor={PAL["ink"]};" vertex="1" parent="1">'
            f'<mxGeometry x="{x + 26}" y="{y - 16}" width="{w}" height="18" as="geometry"/></mxCell>'
        )

    def station(
        self,
        x: float, y: float, w: float, h: float,
        purpose: str, tool: str = "",
        *, stroke: str, fill: str = "#FFFFFF",
        optional: bool = False, purpose_color: str | None = None, extra: str = "",
    ) -> str:
        cid = self.next_id()
        pc = purpose_color or PAL["ink"]
        val = f'&lt;font color=&quot;{pc}&quot;&gt;&lt;b&gt;{html.escape(purpose)}&lt;/b&gt;&lt;/font&gt;{extra}'
        if tool:
            val += (
                f'&lt;br&gt;&lt;font style=&quot;font-size:8px;&quot; color=&quot;{PAL["tool"]}&quot;&gt;'
                f"&lt;i&gt;{html.escape(tool)}&lt;/i&gt;&lt;/font&gt;"
            )
        dash = "dashed=1;dashPattern=6 4;" if optional else ""
        self.cells.append(
            f'<mxCell id="{cid}" value="{val}" style="rounded=1;whiteSpace=wrap;html=1;align=center;'
            f'fillColor={fill};strokeColor={stroke};strokeWidth=1.8;arcSize=10;{dash}" vertex="1" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        return cid

    def fqs_mini_chart(self, x: float, y: float, w: float, h: float) -> str:
        """Fake FastQ Screen panel using animal icons + % bars."""
        gid = self.next_id()
        self.cells.append(
            f'<mxCell id="{gid}" value="" style="group" vertex="1" connectable="0" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        panel = self.next_id()
        self.cells.append(
            f'<mxCell id="{panel}" value="" style="rounded=1;whiteSpace=wrap;html=1;'
            f'fillColor=#FFFFFF;strokeColor={PAL["host"]};strokeWidth=1.5;arcSize=6;" '
            f'vertex="1" parent="{gid}">'
            f'<mxGeometry width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        hdr = self.next_id()
        self.cells.append(
            f'<mxCell id="{hdr}" value="&lt;font color=&quot;{PAL["host"]}&quot;&gt;&lt;b&gt;'
            f'FastQ Screen&lt;/b&gt; &lt;font style=&quot;font-size:7px;&quot;&gt;'
            f'(example)&lt;/font&gt;&lt;/font&gt;" style="text;html=1;align=center;" '
            f'vertex="1" parent="{gid}">'
            f'<mxGeometry y="2" width="{w}" height="14" as="geometry"/></mxCell>'
        )

        # Fake one-hit-one-genome style percentages
        rows = [
            ("cattle", "Cattle", 0.48, "#2E7D52"),
            ("sheep", "Sheep", 0.12, "#376795"),
            ("goat", "Goat", 0.07, "#5B4A8A"),
            ("pig", "Pig", 0.04, "#D4483B"),
            ("dog", "Dog", 0.02, "#8C1C2B"),
            ("horse", "Horse", 0.01, "#8A94A6"),
        ]
        bar_max = w - 78
        row_h = (h - 22) / len(rows)
        for i, (key, label, frac, color) in enumerate(rows):
            ry = 18 + i * row_h
            # icon
            ic = self.next_id()
            uri = self.icon[key].replace("&", "&amp;")
            self.cells.append(
                f'<mxCell id="{ic}" value="" style="shape=image;html=1;imageAspect=1;aspect=fixed;'
                f'image={uri}" vertex="1" parent="{gid}">'
                f'<mxGeometry x="4" y="{ry + 1}" width="{row_h - 3}" height="{row_h - 3}" as="geometry"/></mxCell>'
            )
            # bar background
            bg = self.next_id()
            self.cells.append(
                f'<mxCell id="{bg}" value="" style="rounded=0;fillColor=#EEF2F6;strokeColor=none;" '
                f'vertex="1" parent="{gid}">'
                f'<mxGeometry x="28" y="{ry + row_h * 0.28}" width="{bar_max}" height="{row_h * 0.4}" '
                f'as="geometry"/></mxCell>'
            )
            # bar fill
            bw = max(2, bar_max * frac)
            bar = self.next_id()
            self.cells.append(
                f'<mxCell id="{bar}" value="" style="rounded=0;fillColor={color};strokeColor=none;" '
                f'vertex="1" parent="{gid}">'
                f'<mxGeometry x="28" y="{ry + row_h * 0.28}" width="{bw}" height="{row_h * 0.4}" '
                f'as="geometry"/></mxCell>'
            )
            # % label
            pct = self.next_id()
            self.cells.append(
                f'<mxCell id="{pct}" value="&lt;font style=&quot;font-size:7px;&quot;&gt;'
                f'{int(frac * 100)}%&lt;/font&gt;" style="text;html=1;align=left;" '
                f'vertex="1" parent="{gid}">'
                f'<mxGeometry x="{28 + bw + 2}" y="{ry}" width="28" height="{row_h}" as="geometry"/></mxCell>'
            )
        return gid

    def auth_panel(self, x: float, y: float, w: float, h: float) -> str:
        gid = self.next_id()
        self.cells.append(
            f'<mxCell id="{gid}" value="" style="group" vertex="1" connectable="0" parent="1">'
            f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        outer = self.next_id()
        self.cells.append(
            f'<mxCell id="{outer}" value="" style="rounded=1;whiteSpace=wrap;html=1;'
            f'fillColor={PAL["combo_fill"]};strokeColor={PAL["combo"]};strokeWidth=2.2;arcSize=8;" '
            f'vertex="1" parent="{gid}">'
            f'<mxGeometry width="{w}" height="{h}" as="geometry"/></mxCell>'
        )
        header = self.next_id()
        self.cells.append(
            f'<mxCell id="{header}" value="" style="rounded=0;whiteSpace=wrap;html=1;'
            f'fillColor={PAL["combo_header"]};strokeColor=none;" vertex="1" parent="{gid}">'
            f'<mxGeometry width="{w}" height="22" as="geometry"/></mxCell>'
        )
        title = self.next_id()
        self.cells.append(
            f'<mxCell id="{title}" value="&lt;font color=&quot;#FFFFFF&quot;&gt;&lt;b&gt;'
            f'Authentication&lt;/b&gt;&lt;/font&gt;" style="text;html=1;'
            f'strokeColor=none;fillColor=none;align=center;fontSize=10;" vertex="1" parent="{gid}">'
            f'<mxGeometry y="2" width="{w}" height="18" as="geometry"/></mxCell>'
        )
        sub = self.next_id()
        self.cells.append(
            f'<mxCell id="{sub}" value="&lt;font color=&quot;{PAL["path"]}&quot;&gt;'
            f'&lt;b&gt;Composite score&lt;/b&gt;&lt;/font&gt;" '
            f'style="text;html=1;strokeColor=none;fillColor=none;align=center;fontSize=8;" '
            f'vertex="1" parent="{gid}">'
            f'<mxGeometry y="26" width="{w}" height="14" as="geometry"/></mxCell>'
        )
        chip_w, chip_h, gap = (w - 20 - 3 * 4) / 4, 18, 4
        oy = 44
        for i, metric in enumerate(METRICS):
            row, col = divmod(i, 4)
            n_cols = 4 if row == 0 else 3
            row_w = n_cols * chip_w + (n_cols - 1) * gap
            ox = 10 + (w - 20 - row_w) / 2
            cx = ox + col * (chip_w + gap)
            cy = oy + row * (chip_h + 4)
            chip = self.next_id()
            self.cells.append(
                f'<mxCell id="{chip}" value="&lt;font style=&quot;font-size:7px;&quot;&gt;'
                f'{html.escape(metric)}&lt;/font&gt;" '
                f'style="rounded=1;whiteSpace=wrap;html=1;align=center;'
                f'fillColor={PAL["metric_chip"]};strokeColor={PAL["metric_chip_border"]};'
                f'strokeWidth=1;arcSize=40;" vertex="1" parent="{gid}">'
                f'<mxGeometry x="{cx}" y="{cy}" width="{chip_w}" height="{chip_h}" as="geometry"/></mxCell>'
            )
        return gid

    def edge(self, src: str, dst: str, *, stroke: str, optional: bool = False) -> None:
        eid = self.next_id()
        dash = "dashed=1;dashPattern=6 4;" if optional else ""
        sw = 2.4 if optional else 3.4
        self.cells.append(
            f'<mxCell id="{eid}" value="" style="edgeStyle=orthogonalEdgeStyle;rounded=1;'
            f'orthogonalLoop=1;jettySize=auto;html=1;strokeColor={stroke};strokeWidth={sw};'
            f'endArrow=classic;endFill=1;{dash}" edge="1" parent="1" source="{src}" target="{dst}">'
            f'<mxGeometry relative="1" as="geometry"/></mxCell>'
        )

    def build(self) -> str:
        lx, rx = 12, 430
        lw, rw = 404, 836
        sw, sh = 86, 40

        self.section(lx, 28, lw, 72, "Preprocessing", 1, fill=PAL["prep_fill"])
        self.section(lx, 112, lw, 138, "Host route", 2, fill=PAL["host_fill"])
        self.section(lx, 262, lw, 114, "Metagenomics screening", 3, fill=PAL["meta_fill"])
        self.section(rx, 28, rw, 348, "Pathogen route", 4, fill="#FFFCFC")

        # Section mascots
        self.image(lx + lw - 54, 88, 48, 38, self.icon["cattle"])       # host
        self.image(lx + lw - 52, 236, 46, 46, self.icon["ecoli"])        # metagenomics
        self.image(rx + rw - 60, 4, 54, 54, self.icon["virus"])         # pathogen

        # 1 Preprocessing
        raw = self.station(lx + 14, 48, sw, sh, "Input reads", "FASTQ",
                           stroke=PAL["prep"], purpose_color=PAL["prep"])
        trim = self.station(lx + 116, 48, sw, sh, "Trim adapters", "AdapterRemoval",
                            stroke=PAL["prep"], purpose_color=PAL["prep"])

        # 2 Host — identification + mini FQS chart under it
        host_id = self.station(
            lx + 10, 128, 100, 40, "Host identification", "FastQ Screen",
            stroke=PAL["host"], fill=PAL["host_fill"], purpose_color=PAL["host"],
        )
        fqs = self.fqs_mini_chart(lx + 10, 172, 100, 70)

        hmap = self.station(lx + 122, 128, sw, sh, "Map host", "BWA / Bowtie2",
                            stroke=PAL["host"], purpose_color=PAL["host"])
        mtdna = self.station(lx + 122, 176, sw, sh, "Map mtDNA", "BWA / Bowtie2",
                             stroke=PAL["host"], purpose_color=PAL["host"])
        hmerge = self.station(lx + 220, 152, sw, sh, "Merge", "samtools",
                              stroke=PAL["host"], purpose_color=PAL["host"])
        hqual = self.station(lx + 318, 128, sw, sh, "Mapping QC", "Qualimap",
                             stroke=PAL["host"], purpose_color=PAL["host"])
        hdmg = self.station(lx + 318, 176, sw, sh, "Damage", "DamageProfiler",
                            stroke=PAL["host"], purpose_color=PAL["host"])
        endog = self.station(lx + 220, 210, sw, 32, "Endogenous DNA", "host fraction",
                             stroke=PAL["host"], fill=PAL["host_fill"], purpose_color=PAL["host"])
        sexing = self.station(lx + 122, 210, sw, 32, "Genetic sexing", "sexing_residual",
                              stroke=PAL["opt"], fill=PAL["opt_fill"], optional=True)

        # 3 Metagenomics
        prinseq = self.station(lx + 10, 282, sw, sh, "Filter & dedup", "PRINSEQ++",
                               stroke=PAL["meta"], purpose_color=PAL["meta"])
        unaligned = self.station(lx + 108, 282, sw, sh, "Mammals reads removal", "Bowtie2",
                                 stroke=PAL["meta"], purpose_color=PAL["meta"])
        pool = self.station(lx + 206, 282, sw, sh, "Pool", "cat · pigz",
                            stroke=PAL["meta"], purpose_color=PAL["meta"])
        kraken = self.station(lx + 304, 282, 96, sh, "Classify", "KrakenUniq",
                              stroke=PAL["meta"], purpose_color=PAL["meta"])
        malt = self.station(lx + 304, 330, 96, 32, "Metagenome screen", "MALT",
                            stroke=PAL["opt"], fill=PAL["opt_fill"], optional=True)
        decom = self.station(lx + 206, 330, sw, 32, "Source screen", "decOM",
                             stroke=PAL["opt"], fill=PAL["opt_fill"], optional=True)

        # 4 Pathogen sequence labels
        for x_off, label, col in (
            (20, "1 · Detection", PAL["detect"]),
            (170, "2 · Mapping", PAL["path"]),
            (310, "3 · Authentication", PAL["combo"]),
            (620, "4 · Outputs", PAL["out"]),
        ):
            z = self.next_id()
            self.cells.append(
                f'<mxCell id="{z}" value="&lt;font color=&quot;{col}&quot;&gt;&lt;b&gt;'
                f'{html.escape(label)}&lt;/b&gt;&lt;/font&gt;" style="text;html=1;strokeColor=none;'
                f'fillColor=none;fontSize=9;" vertex="1" parent="1">'
                f'<mxGeometry x="{rx + x_off}" y="42" width="140" height="14" as="geometry"/></mxCell>'
            )

        detect = self.station(
            rx + 16, 150, 130, 58,
            "Pathogen detection",
            "E-value · MaltExtract",
            stroke=PAL["detect"], fill=PAL["detect_fill"], purpose_color=PAL["detect"],
        )
        # small pathogen icons near detection
        self.image(rx + 24, 88, 38, 38, self.icon["brucella"])
        self.image(rx + 68, 86, 40, 40, self.icon["virus"])

        pmap = self.station(
            rx + 168, 150, 110, 58,
            "Pathogen mapping",
            "BWA / Bowtie2",
            stroke=PAL["path"], fill=PAL["path_fill"], purpose_color=PAL["path"],
        )
        auth = self.auth_panel(rx + 300, 100, 290, 110)
        final = self.station(
            rx + 620, 158, 100, 42,
            "Cohort reports", "Excel · PDF",
            stroke=PAL["out"], fill=PAL["out_fill"], purpose_color=PAL["out"],
        )

        # edges
        self.edge(raw, trim, stroke=PAL["prep"])
        self.edge(trim, host_id, stroke=PAL["host"])
        self.edge(trim, prinseq, stroke=PAL["meta"])

        self.edge(host_id, hmap, stroke=PAL["host"])
        self.edge(host_id, mtdna, stroke=PAL["host"])
        self.edge(hmap, hmerge, stroke=PAL["host"])
        self.edge(mtdna, hmerge, stroke=PAL["host"])
        self.edge(hmerge, hqual, stroke=PAL["host"])
        self.edge(hmerge, hdmg, stroke=PAL["host"])
        self.edge(hmerge, endog, stroke=PAL["host"])
        self.edge(hqual, sexing, stroke=PAL["opt"], optional=True)

        self.edge(prinseq, unaligned, stroke=PAL["meta"])
        self.edge(unaligned, pool, stroke=PAL["meta"])
        self.edge(pool, kraken, stroke=PAL["meta"])
        self.edge(pool, malt, stroke=PAL["opt"], optional=True)
        self.edge(kraken, decom, stroke=PAL["opt"], optional=True)

        self.edge(kraken, detect, stroke=PAL["detect"])
        self.edge(malt, detect, stroke=PAL["opt"], optional=True)
        self.edge(detect, pmap, stroke=PAL["path"])
        self.edge(pmap, auth, stroke=PAL["path"])
        self.edge(auth, final, stroke=PAL["out"])
        self.edge(hqual, final, stroke=PAL["host"])
        self.edge(endog, final, stroke=PAL["host"], optional=True)

        # legend
        leg_y = PH - 28
        items = [
            (PAL["prep"], "Preprocess", False),
            (PAL["host"], "Host", False),
            (PAL["meta"], "Metagenomics", False),
            (PAL["detect"], "Detection", False),
            (PAL["path"], "Mapping", False),
            (PAL["combo"], "Authentication", False),
            (PAL["opt"], "Optional", True),
        ]
        lx_leg = 12
        for col, label, dashed in items:
            lid = self.next_id()
            dash = "dashed=1;dashPattern=6 4;" if dashed else ""
            self.cells.append(
                f'<mxCell id="{lid}" value="" style="endArrow=none;html=1;strokeColor={col};'
                f'strokeWidth=4;{dash}" edge="1" parent="1">'
                f'<mxGeometry relative="1" as="geometry">'
                f'<mxPoint x="{lx_leg}" y="{leg_y}" as="sourcePoint"/>'
                f'<mxPoint x="{lx_leg + 18}" y="{leg_y}" as="targetPoint"/>'
                f"</mxGeometry></mxCell>"
            )
            tid = self.next_id()
            self.cells.append(
                f'<mxCell id="{tid}" value="{html.escape(label)}" '
                f'style="text;html=1;fontSize=8;fontColor=#5C6578;align=left;" vertex="1" parent="1">'
                f'<mxGeometry x="{lx_leg + 22}" y="{leg_y - 6}" width="90" height="12" as="geometry"/></mxCell>'
            )
            lx_leg += 108

        title = self.next_id()
        self.cells.insert(
            0,
            f'<mxCell id="{title}" value="&lt;b&gt;PIGSTI&lt;/b&gt; '
            f'&lt;font style=&quot;font-size:9px;color:#5C6578;&quot;&gt;'
            f'aDNA host · metagenomics · pathogen detection → mapping → authentication&lt;/font&gt;" '
            f'style="text;html=1;fontSize=14;align=left;" vertex="1" parent="1">'
            f'<mxGeometry x="12" y="4" width="700" height="20" as="geometry"/></mxCell>',
        )

        body = "\n        ".join(self.cells)
        return f"""<mxfile host="app.diagrams.net" agent="PIGSTI" version="22.1.0" type="device">
  <diagram name="PIGSTI workflow" id="{_id()}">
    <mxGraphModel dx="{PW}" dy="{PH}" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1"
      fold="1" page="1" pageScale="1" pageWidth="{PW}" pageHeight="{PH}" math="0" shadow="0">
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
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(DrawioBuilder().build(), encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"Icons from: {ICONS} and {MICROBES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
