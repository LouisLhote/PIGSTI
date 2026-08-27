#!/usr/bin/env python3
"""Rasterize SVG icons to PNG via Playwright Chromium (no Cairo required).

Output PNGs are tightly cropped silhouettes on a transparent background.
"""

from __future__ import annotations

import base64
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
ICON_DIRS = [
    ROOT / "docs" / "images" / "Icons",
    ROOT / "docs" / "images" / "Icons microbes",
]


def make_transparent_silhouette(png_path: Path, size: int = 160) -> None:
    """Key out light background; keep dark ink; crop & pad on transparent canvas."""
    im = Image.open(png_path).convert("RGBA")
    a = np.asarray(im).copy()
    r, g, b, al = (a[:, :, i].astype(np.float32) for i in range(4))
    lum = (r + g + b) / 3.0
    out_a = al.copy()
    out_a[lum >= 240] = 0
    soft = (lum >= 200) & (lum < 240) & (al > 0)
    out_a[soft] = al[soft] * ((240.0 - lum[soft]) / 40.0)
    # Never invent opacity on already-transparent pixels (RGB often 0,0,0).
    drawn = (al > 8) & (lum < 200)
    out_a[drawn] = 255
    a[:, :, 0:3][drawn] = 0
    a[:, :, 3] = np.clip(out_a, 0, 255).astype(np.uint8)
    a[a[:, :, 3] == 0, 0:3] = 0
    # extra pass: any remaining light ink → transparent
    lum2 = a[:, :, 0:3].astype(np.float32).mean(axis=2)
    a[(lum2 > 180) & (a[:, :, 3] > 0), 3] = 0

    out = Image.fromarray(a)
    bbox = out.getbbox()
    if bbox is None:
        out.save(png_path)
        return
    cropped = out.crop(bbox)
    side = int(max(cropped.size) * 1.08) or size
    canvas = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    ox = (side - cropped.size[0]) // 2
    oy = (side - cropped.size[1]) // 2
    canvas.paste(cropped, (ox, oy), cropped)
    canvas.resize((size, size), Image.Resampling.LANCZOS).save(png_path)


def svg_to_png(svg_path: Path, png_path: Path, size: int = 160) -> None:
    from playwright.sync_api import sync_playwright

    raw = svg_path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    data_uri = f"data:image/svg+xml;base64,{b64}"
    html = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
        "<body style='margin:0;background:transparent;'>"
        f"<img id='i' src='{data_uri}' "
        f"style='width:{size}px;height:{size}px;object-fit:contain;display:block;'/>"
        "</body></html>"
    )
    png_path.parent.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": size + 8, "height": size + 8})
        page.set_content(html, wait_until="load")
        page.wait_for_timeout(150)
        page.locator("#i").screenshot(path=str(png_path), omit_background=True)
        browser.close()
    make_transparent_silhouette(png_path, size=size)


def main() -> int:
    size = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    n_ok = 0
    for d in ICON_DIRS:
        out = d / "png"
        for svg in sorted(d.glob("*.svg")):
            if svg.stat().st_size > 400_000 and svg.stem.lower() == "pox":
                print(f"skip {svg.name} (too large)")
                continue
            png = out / f"{svg.stem}.png"
            try:
                svg_to_png(svg, png, size=size)
                print(f"ok {png.relative_to(ROOT)} ({png.stat().st_size} B)")
                n_ok += 1
            except Exception as e:
                msg = str(e).encode("ascii", "replace").decode("ascii")
                print(f"err {svg.name}: {type(e).__name__}: {msg}")
    print(f"Rasterized {n_ok} icons (transparent)")
    return 0 if n_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
