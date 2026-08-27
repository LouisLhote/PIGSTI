# PIGSTI workflow figure (nf-core metro-map style)

Publication workflow schematics can use the **official nf-core component library** (CC0 / public domain). This repo ships the library and two ways to build the figure.

## Option A — nf-metro (recommended for papers)

[nf-metro](https://seqeralabs.github.io/nf-metro/) renders true nf-core metro maps from a Mermaid source, including **section boxes** (Preprocessing, Host route, Metagenomics route, Pathogen route).

```bash
pip install nf-metro
python scripts/render_subway_workflow.py          # matplotlib fallback + nf-metro if installed
# or directly:
nf-metro render docs/images/pigsti_subway_workflow.mmd -o docs/images/pigsti_subway_workflow.svg --theme nfcore
nf-metro render docs/images/pigsti_subway_workflow.mmd -o docs/images/pigsti_subway_workflow.html --format html --theme nfcore
```

**Outputs**

| File | Use |
|------|-----|
| `docs/images/pigsti_subway_workflow.svg` | Dark nf-core theme (screen / slides) |
| `docs/images/pigsti_subway_workflow_light.svg` | Light theme (print / Word / Google Docs) |
| `docs/images/pigsti_subway_workflow.png` | Raster preview |
| `docs/images/pigsti_subway_workflow.pdf` | Print / submission |
| `docs/images/pigsti_subway_workflow.html` | Interactive preview (pan, zoom, line legend) |
| `docs/images/pigsti_subway_workflow.mmd` | Editable source |

Edit stations, colours, and section layout in `pigsti_subway_workflow.mmd`, then re-run the render command.

## Option B — draw.io with nf-core shapes (manual, pixel-perfect)

The official shape library is already in this repo:

**[`docs/nf-core_components.xml`](nf-core_components.xml)**

### Load the library in draw.io

1. Open [https://app.diagrams.net/](https://app.diagrams.net/)
2. **File → Open Library from → Device**
3. Select `PIGSTI_publication/docs/nf-core_components.xml`
4. nf-core stations, tracks, file icons, and elbows appear in the left sidebar

You can also download the latest library from nf-core:

<https://raw.githubusercontent.com/nf-core/website/refs/heads/main/sites/docs/src/assets/images/graphic_design_assets/workflow_schematics_components/generic/nf-core_components.xml>

### Build the diagram

1. Open **`docs/drawio/pigsti_workflow.drawio`** in draw.io (or regenerate with `python scripts/build_drawio_workflow.py`)
2. Load **`docs/nf-core_components.xml`** for official nf-core station shapes
3. Replace rounded boxes with nf-core shapes; add animal icons at **Host identification**
4. **File → Export as → SVG** (or PNG/PDF) into `docs/images/`

### Paste into Google Docs

1. Export **SVG** or **PNG** from draw.io (SVG scales better)
2. In Google Docs: **Insert → Image → Upload from computer**
3. Or copy the PNG from the export folder and paste (Ctrl+V)

For Word / PowerPoint, prefer **SVG** or high-DPI **PNG** (300 dpi export in draw.io: **File → Export → Zoom 300%**).

## Attribution

nf-core metro-map components: James A. Fellows Yates, Maxime Garcia, Louis Le Nézet & nf-core — [CC0](https://creativecommons.org/publicdomain/zero/1.0/).

When adapting nf-core example schematics (e.g. nf-core/eager), check the specific licence on [nf-co.re workflow schematics](https://nf-co.re/docs/community/brand/workflow-schematics).
