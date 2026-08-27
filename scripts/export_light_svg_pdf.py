#!/usr/bin/env python3
"""Convert nf-metro light SVGs to high-quality vector-friendly PDFs via Chromium."""

from __future__ import annotations

import argparse
import base64
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
IMG = ROOT / "docs" / "images"


def svg_dimensions(svg_text: str) -> tuple[float, float]:
    m = re.search(r'viewBox="([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)"', svg_text)
    if m:
        return float(m.group(3)), float(m.group(4))
    wm = re.search(r'\bwidth="([0-9.]+)"', svg_text)
    hm = re.search(r'\bheight="([0-9.]+)"', svg_text)
    if wm and hm:
        return float(wm.group(1)), float(hm.group(1))
    return 2000.0, 650.0


def svg_to_pdf(svg_path: Path, pdf_path: Path, *, scale: float = 2.0) -> None:
    """Render SVG in Chromium and print to PDF (preserves embedded PNG icons)."""
    from playwright.sync_api import sync_playwright

    raw = svg_path.read_bytes()
    text = raw.decode("utf-8")
    w, h = svg_dimensions(text)
    # Scale page so print is sharp (CSS px → PDF points via deviceScaleFactor)
    page_w = max(1, int(w * scale))
    page_h = max(1, int(h * scale))

    b64 = base64.b64encode(raw).decode("ascii")
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/>
<style>
  @page {{ size: {page_w}px {page_h}px; margin: 0; }}
  html, body {{ margin: 0; padding: 0; background: #ffffff; width: {page_w}px; height: {page_h}px; }}
  img {{ width: {page_w}px; height: {page_h}px; display: block; }}
</style></head>
<body><img id="s" src="data:image/svg+xml;base64,{b64}"/></body></html>
"""
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(
            viewport={"width": page_w, "height": page_h},
            device_scale_factor=2,
        )
        page.set_content(html, wait_until="networkidle")
        page.wait_for_timeout(300)
        page.pdf(
            path=str(pdf_path),
            width=f"{page_w}px",
            height=f"{page_h}px",
            print_background=True,
            margin={"top": "0", "right": "0", "bottom": "0", "left": "0"},
            prefer_css_page_size=True,
        )
        browser.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Export *_light.svg → high-quality PDF")
    ap.add_argument("--indir", type=Path, default=IMG)
    ap.add_argument("--scale", type=float, default=2.0, help="Render scale (2 = high-res)")
    ap.add_argument("files", nargs="*", help="Specific SVG paths (default: all *_light.svg)")
    args = ap.parse_args()

    files = [Path(f) for f in args.files] if args.files else sorted(args.indir.glob("pigsti_subway*_light.svg"))
    if not files:
        print("No light SVGs found")
        return 1

    for svg in files:
        if not svg.exists():
            print(f"missing {svg}")
            continue
        pdf = svg.with_name(svg.stem + ".pdf")  # e.g. pigsti_subway_workflow_light.pdf
        svg_to_pdf(svg, pdf, scale=args.scale)
        print(f"Wrote {pdf} ({pdf.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
