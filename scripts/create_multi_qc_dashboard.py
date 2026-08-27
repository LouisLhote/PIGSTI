#!/usr/bin/env python3
"""
Create a lightweight multi-QC HTML dashboard (per biological sample).

The dashboard combines:
  - Host / mtDNA QC metrics from results/final/host_mtdna_summary_all_samples.xlsx
    (uses the "Sample_level" sheet)
  - Pathogen QC / detection metrics from results/final/pathogen_summary_all_samples.xlsx

It also links to per-pathogen PDF reports and (optionally) edit-distance plots.
"""

import os
import sys
import html
import math
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pigsti_naming import safe_pathogen_name as pathogen_safe


def _safe_str(x):
    if x is None:
        return "NA"
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return "NA"
    if str(x).strip() == "" or str(x).strip().lower() in {"nan", "none", "na"}:
        return "NA"
    return str(x)


def _parse_score_num(score):
    """
    Score is stored as a string like "3/9". Return numerator as float for sorting.
    """
    if score is None:
        return -1.0
    s = str(score).strip()
    if not s or s.lower() in {"na", "nan", "none"}:
        return -1.0
    if "/" in s:
        try:
            return float(s.split("/", 1)[0])
        except Exception:
            return -1.0
    try:
        return float(s)
    except Exception:
        return -1.0


def load_host_sample_level(host_excel: str) -> pd.DataFrame:
    # summarize_host_mtdna.py writes both PCR_level and Sample_level
    df = pd.read_excel(host_excel, sheet_name="Sample_level")
    # Normalize column names just in case
    df.columns = [str(c).strip() for c in df.columns]
    return df


def load_pathogen_table(pathogen_excel: str) -> pd.DataFrame:
    try:
        df = pd.read_excel(pathogen_excel, sheet_name="pathogen_summary")
    except ValueError:
        df = pd.read_excel(pathogen_excel, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def main():
    try:
        host_excel = snakemake.input.host_mtdna
        pathogen_excel = snakemake.input.pathogen
        out_html = snakemake.output.html
        top_n = int(snakemake.params.top_pathogens)
        results_dir = snakemake.params.results_dir
    except Exception:
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--host-mtdna", required=True)
        parser.add_argument("--pathogen", required=True)
        parser.add_argument("--output-html", required=True)
        parser.add_argument("--top-pathogens", type=int, default=5)
        parser.add_argument("--results-dir", default="results")
        args = parser.parse_args()
        host_excel = args.host_mtdna
        pathogen_excel = args.pathogen
        out_html = args.output_html
        top_n = args.top_pathogens
        results_dir = args.results_dir

    host_df = load_host_sample_level(host_excel)
    pathogen_df = load_pathogen_table(pathogen_excel)

    # Expected columns:
    # - host_df: bio_sample, species, raw_reads, collapsed_reads, host_q30_reads, ...
    # - pathogen_df: Sample, Pathogen, Score, Coverage, Mapped_reads, ANI, ...
    if "bio_sample" not in host_df.columns:
        raise ValueError(f"Expected host Sample_level sheet to include 'bio_sample'. Columns: {list(host_df.columns)}")
    if "Sample" not in pathogen_df.columns or "Pathogen" not in pathogen_df.columns:
        raise ValueError(
            "Expected pathogen summary Excel to include 'Sample' and 'Pathogen'. "
            f"Columns: {list(pathogen_df.columns)}"
        )

    # Basic normalization
    host_df["bio_sample"] = host_df["bio_sample"].astype(str)
    pathogen_df["Sample"] = pathogen_df["Sample"].astype(str)
    pathogen_df["Pathogen"] = pathogen_df["Pathogen"].astype(str)

    # Keep only columns we might show (safe even if absent)
    host_cols = [
        "species",
        "raw_reads",
        "collapsed_reads",
        "host_q30_reads",
        "host_dedup_reads",
        "host_coverage",
        "host_endogenous_pct",
        "duplication_rate",
        "mt_q30_reads",
        "mt_dedup_reads",
        "mtdna_coverage",
        "damageprofiler_5pCtoT_mean",
        "human_contamination_pct",
        "other_warnings",
        "host_mean_read_length",
        "host_median_read_length",
        "mt_mean_read_length",
        "mt_median_read_length",
    ]
    host_cols = [c for c in host_cols if c in host_df.columns]

    pathogen_cols_preferred = [
        "Pathogen",
        "Score",
        "Coverage",
        "Mapped_reads",
        "Read_mapping_ratio",
        "ANI",
        "Relative_entropy",
        "Breadth_ratio",
        "Edit_distance_decay_quality",
        "Edit_distance_decay_quality",
        "Criteria_passed",
        "Damage_5p_CtoT",
        "Genus_ranking",
    ]
    pathogen_cols = [c for c in pathogen_cols_preferred if c in pathogen_df.columns]
    if "Score" not in pathogen_cols and "Score" in pathogen_df.columns:
        pathogen_cols.append("Score")

    # Build per-sample cards
    bios = sorted(host_df["bio_sample"].unique().tolist())

    cards = []
    for bio in bios:
        host_row = host_df[host_df["bio_sample"] == bio]
        if host_row.empty:
            continue
        host_row = host_row.iloc[0].to_dict()

        sample_species = _safe_str(host_row.get("species", "NA"))

        # Pathogens for this sample
        sub = pathogen_df[pathogen_df["Sample"] == bio].copy()
        if not sub.empty and "Score" in sub.columns:
            sub["_score_num"] = sub["Score"].apply(_parse_score_num)
            sub = sub.sort_values("_score_num", ascending=False)
        sub_top = sub.head(top_n) if not sub.empty else sub

        # Links & plot paths for each pathogen
        pathogen_rows_html = []
        for _, prow in sub_top.iterrows():
            pname = str(prow["Pathogen"])
            psafe = pathogen_safe(pname)

            pdf = os.path.join(results_dir, bio, "summary", f"{bio}_{psafe}_pathogen_report.pdf")
            edit_damage_png = os.path.join(results_dir, bio, "pathogen_mapping", f"{bio}_{psafe}.edit_distance_logr2_damaged.png")
            edit_default_png = os.path.join(
                results_dir, bio, "pathogen_mapping", f"{bio}_{psafe}.edit_distance_logr2_default.png"
            )
            edit_no_png = edit_default_png
            if not os.path.isfile(edit_no_png):
                edit_no_png = os.path.join(
                    results_dir, bio, "pathogen_mapping", f"{bio}_{psafe}.edit_distance_logr2_no_damage.png"
                )

            pdf_rel = pdf.replace("\\", "/")
            img_src = None
            # Prefer damage image if available; else no-damage; else legacy
            if os.path.exists(edit_damage_png):
                img_src = edit_damage_png.replace("\\", "/")
            elif os.path.exists(edit_no_png):
                img_src = edit_no_png.replace("\\", "/")
            else:
                legacy = os.path.join(results_dir, bio, "pathogen_mapping", f"{bio}_{psafe}.edit_distance_logr2.png")
                if os.path.exists(legacy):
                    img_src = legacy.replace("\\", "/")

            # Convert to relative to HTML location (results/final/multi_qc_dashboard.html)
            # i.e., strip leading "results/" and prefix "./"
            def to_rel(path):
                p = str(path).replace("\\", "/")
                if p.startswith(results_dir + "/"):
                    p = "./" + p[len(results_dir) + 1:]
                elif p == results_dir:
                    p = "./"
                return p

            row_cells = []
            # Build a small subset row: Pathogen | Score | ANI | Edit
            for col in pathogen_cols:
                if col == "Pathogen":
                    cell = f"<b>{html.escape(pname)}</b>"
                    links = [f'<a href="{html.escape(to_rel(pdf))}">PDF</a>'] if os.path.exists(pdf) else []
                    if links:
                        cell += "<br/>" + " | ".join(links)
                    if img_src is not None and os.path.exists(img_src):
                        cell += f'<br/><img src="{html.escape(to_rel(img_src))}" style="max-width:240px; height:auto; border:1px solid #ddd; border-radius:4px;" />'
                    row_cells.append(f"<td>{cell}</td>")
                else:
                    row_cells.append(f"<td>{html.escape(_safe_str(prow.get(col, 'NA')))}</td>")

            pathogen_rows_html.append("<tr>" + "".join(row_cells) + "</tr>")

        # Host table
        host_kvs = []
        for col in host_cols:
            host_kvs.append(f"<tr><th>{html.escape(col)}</th><td>{html.escape(_safe_str(host_row.get(col)))}</td></tr>")

        host_table_html = "<table class='table host-table'>" + "".join(host_kvs) + "</table>"

        pathogens_table = ""
        if sub_top is not None and len(sub_top) > 0 and pathogen_cols:
            # Create header from pathogen_cols
            header = "".join([f"<th>{html.escape(c)}</th>" for c in pathogen_cols if c != "Pathogen"])
            # But we handle Pathogen col specially by injecting it first.
            pathogens_table = (
                "<table class='table pathogen-table'>"
                "<thead><tr>"
                "<th>Pathogen</th>"
                + "".join([f"<th>{html.escape(c)}</th>" for c in pathogen_cols if c != "Pathogen"])
                + "</tr></thead><tbody>"
                + "".join(pathogen_rows_html)
                + "</tbody></table>"
            )
        else:
            pathogens_table = "<p><i>No pathogen records in this summary for this sample.</i></p>"

        card = f"""
        <section class="sample-card" id="{html.escape(bio)}">
          <h2>Sample: {html.escape(bio)} <span class="muted">({html.escape(sample_species)})</span></h2>
          {host_table_html}
          <h3>Top pathogens (by Score)</h3>
          {pathogens_table}
        </section>
        """
        cards.append(card)

    html_doc = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>PIGSTI Multi-QC Dashboard (per biological sample)</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 30px; color: #1f2937; }}
    .muted {{ color: #6b7280; font-weight: normal; font-size: 0.9em; }}
    .sample-card {{ border: 1px solid #e5e7eb; border-radius: 10px; padding: 18px; margin-bottom: 22px; }}
    h1 {{ margin-top: 0; }}
    h2 {{ margin-bottom: 10px; }}
    h3 {{ margin-top: 16px; margin-bottom: 10px; }}
    .table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
    .table th, .table td {{ border: 1px solid #e5e7eb; padding: 8px; vertical-align: top; }}
    .table th {{ background: #f9fafb; text-align: left; width: 240px; }}
    .host-table th {{ width: 260px; }}
    .host-table td {{ width: auto; }}
    .pathogen-table th {{ background: #f3f4f6; }}
    .pathogen-table img {{ display:block; margin-top: 6px; }}
    a {{ color: #2563eb; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
  </style>
</head>
<body>
  <h1>PIGSTI Multi-QC Dashboard</h1>
  <p><b>Scope:</b> per biological sample (<code>BIO_SAMPLES</code>) with pathogen info from pathogen summary.</p>
  <p><b>Generated:</b> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

  {''.join(cards)}
</body>
</html>
"""

    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html_doc)

    print(f"[OK] Multi-QC dashboard written to: {out_html}")


if __name__ == "__main__":
    main()

