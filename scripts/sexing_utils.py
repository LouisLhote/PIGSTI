"""Helpers for reading and aggregating genetic sexing results."""

from __future__ import annotations

import csv
import os
from pathlib import Path

SUPPORTED_SEXING_SPECIES = {"Cow", "Goat", "Sheep", "Dog"}


def read_sexing_tsv(tsv_path: str) -> dict[str, str | float]:
    """Read a per-library sexing TSV written by run_pigsti_sexing.py / R scripts."""
    defaults: dict[str, str | float] = {
        "sexing_call": "NA",
        "sexing_female_prob": "NA",
        "sexing_likelihood_ratio": "NA",
        "sexing_status": "missing",
        "sexing_note": "",
    }
    if not tsv_path or not os.path.isfile(tsv_path):
        return defaults

    try:
        with open(tsv_path, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError as exc:
        defaults["sexing_status"] = "read_error"
        defaults["sexing_note"] = str(exc)
        return defaults
    except csv.Error as exc:
        defaults["sexing_status"] = "read_error"
        defaults["sexing_note"] = str(exc)
        return defaults

    if not rows:
        defaults["sexing_status"] = "empty"
        return defaults

    row = rows[0]
    status = (row.get("status") or "ok").strip()
    call = (row.get("call") or "NA").strip()

    out = dict(defaults)
    out["sexing_status"] = status
    out["sexing_call"] = call
    if row.get("female_prob") not in (None, "", "NA"):
        try:
            out["sexing_female_prob"] = float(row["female_prob"])
        except ValueError:
            out["sexing_female_prob"] = row["female_prob"]
    if row.get("likelihood_ratio") not in (None, "", "NA"):
        try:
            out["sexing_likelihood_ratio"] = float(row["likelihood_ratio"])
        except ValueError:
            out["sexing_likelihood_ratio"] = row["likelihood_ratio"]
    if row.get("note"):
        out["sexing_note"] = str(row["note"])
    elif status == "skipped" and call == "NA":
        out["sexing_note"] = "not_applicable"
    return out


def sexing_plot_path(bio_sample: str) -> str:
    """Return relative path to biological-sample sexing PDF under results/samples/."""
    return f"results/samples/{bio_sample}/sexing/{bio_sample}_sexing.pdf"


def resolve_bio_host_species(species_files: list[str], mismatch_warning: str | None = None) -> tuple[str, str]:
    """
    Pick one host species for a biological sample from PCR-level FastQ Screen results.

    Returns (species, note). Empty note means OK to proceed.
    """
    if mismatch_warning and os.path.isfile(mismatch_warning) and os.path.getsize(mismatch_warning) > 0:
        return "No species found", "host_species_mismatch"

    values: list[str] = []
    for path in species_files:
        if not path or not os.path.isfile(path):
            continue
        s = Path(path).read_text(encoding="utf-8").strip()
        if s and s not in ("NA", "No species found"):
            values.append(s)

    if not values:
        return "No species found", "no_species"

    unique = sorted(set(values))
    if len(unique) > 1:
        return unique[0], f"fastq_screen_species_mismatch:{','.join(unique)}"

    return unique[0], ""
