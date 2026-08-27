#!/usr/bin/env python3
"""Parse FastQ Screen collapsed reports and select the best host species."""

import argparse
import os
import sys


def parse_headers_and_data(file_handle):
    # Skip metadata lines starting with #
    for line in file_handle:
        if not line.startswith("#"):
            header_line = line.strip()
            break
    else:
        raise ValueError("No header line found")

    headers = header_line.split()
    required_columns = [
        "Genome",
        "#Unmapped",
        "#One_hit_one_genome",
        "%Unmapped",
        "%One_hit_one_genome",
        "%Multiple_hits_one_genome",
    ]

    col_indices = {}
    for col in required_columns:
        if col not in headers:
            raise ValueError(f"Missing required column: {col}")
        col_indices[col] = headers.index(col)

    data = []
    for line in file_handle:
        if not line.strip():
            continue
        cols = line.strip().split()
        try:
            genome = cols[col_indices["Genome"]]
            unmapped = float(cols[col_indices["%Unmapped"]])
            one_hit = int(cols[col_indices["#One_hit_one_genome"]])
            one_hit_pct = float(cols[col_indices["%One_hit_one_genome"]])
            multi_hit_pct = float(cols[col_indices["%Multiple_hits_one_genome"]])
            data.append(
                {
                    "genome": genome,
                    "unmapped": unmapped,
                    "one_hit": one_hit,
                    "one_hit_pct": one_hit_pct,
                    "multi_hit_pct": multi_hit_pct,
                }
            )
        except (IndexError, ValueError):
            continue
    return data


def write_warning(path, message):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", encoding="utf-8") as warn:
        warn.write(message + "\n")


def resolve_best_species(species_data, exclude_human=True):
    """
    Pick best species using the same rules as the legacy PIGSTI parser.

    Returns dict with keys: best_species, best_one_hit, best_hit_genome, warnings (list).
    """
    warnings: list[str] = []

    if not species_data:
        return {
            "best_species": "No species found",
            "best_one_hit": 0,
            "best_hit_genome": None,
            "warnings": warnings,
        }

    human_entry = next((s for s in species_data if s["genome"].lower() == "human"), None)
    if human_entry:
        one_hit_pct = human_entry.get("one_hit_pct")
        multi_hit_pct = human_entry.get("multi_hit_pct")
        if one_hit_pct is not None and multi_hit_pct is not None:
            human_contam = one_hit_pct + multi_hit_pct
            warnings.append(f"Human contamination load: {human_contam:.2f}%")
        else:
            warnings.append("Could not compute human contamination load")

    filtered_data = [
        s for s in species_data if not (exclude_human and s["genome"].lower() == "human")
    ]

    if not filtered_data:
        return {
            "best_species": "No species found",
            "best_one_hit": 0,
            "best_hit_genome": None,
            "warnings": warnings,
        }

    best_hit = max(filtered_data, key=lambda x: x["one_hit"])
    best_unmapped = min(filtered_data, key=lambda x: x["unmapped"])

    if best_hit["genome"] != best_unmapped["genome"]:
        best_species = best_unmapped["genome"]
        warnings.append(
            "Best one-hit ({}) and lowest unmapped ({}) do not match. Using {}.".format(
                best_hit["genome"], best_unmapped["genome"], best_species
            )
        )
    else:
        best_species = best_hit["genome"]

    non_best = [s for s in filtered_data if s["genome"].lower() != best_hit["genome"].lower()]
    if non_best:
        second_best = max(non_best, key=lambda x: x["one_hit"])
        if second_best["one_hit"] > 0 and (best_hit["one_hit"] / second_best["one_hit"]) < 1.5:
            warnings.append(
                "Probable {} contamination: {} only {} vs {} {}.".format(
                    second_best["genome"],
                    best_hit["genome"],
                    best_hit["one_hit"],
                    second_best["genome"],
                    second_best["one_hit"],
                )
            )

    return {
        "best_species": best_species,
        "best_one_hit": int(best_hit["one_hit"]),
        "best_hit_genome": best_hit["genome"],
        "warnings": warnings,
    }


def load_species_data(input_file):
    with open(input_file, encoding="utf-8") as f_in:
        return parse_headers_and_data(f_in)


def main_cli():
    parser = argparse.ArgumentParser(description="Parse FastQ Screen collapsed_screen.txt")
    parser.add_argument("screen_txt", nargs="?", help="collapsed_screen.txt (Snakemake mode if omitted)")
    parser.add_argument("best_species_out", nargs="?", help="Output best_species.txt")
    parser.add_argument("--screen-txt", dest="screen_txt_flag", help="collapsed_screen.txt path")
    parser.add_argument("--output", dest="output_flag", help="best_species.txt path")
    parser.add_argument("--sample", default="manual_run", help="Sample ID for warning logs")
    parser.add_argument("--exclude-human", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--print-best-one-hit",
        action="store_true",
        help="Print #One_hit_one_genome for the selected best species (stdout only)",
    )
    parser.add_argument(
        "--min-one-hit",
        type=int,
        default=None,
        help="With --check-rescreen: exit 0 if best_one_hit >= min, else exit 1",
    )
    parser.add_argument(
        "--check-rescreen",
        action="store_true",
        help="Exit 1 when best #One_hit_one_genome is below --min-one-hit (for Snakemake shell)",
    )
    args = parser.parse_args()

    screen_txt = args.screen_txt_flag or args.screen_txt
    if not screen_txt:
        parser.error("Provide collapsed_screen.txt path")

    warning_file = f"logs/contamination_warnings/{args.sample}.txt"

    try:
        species_data = load_species_data(screen_txt)
    except (OSError, ValueError) as e:
        if args.print_best_one_hit or args.check_rescreen:
            print(0)
            return 1 if args.check_rescreen else 0
        msg = (
            f"FastQ Screen parsing issue for {args.sample} ({screen_txt}): {e}. "
            "Treating as 'No species found'."
        )
        write_warning(warning_file, msg)
        species_data = []

    result = resolve_best_species(species_data, exclude_human=args.exclude_human)

    if args.print_best_one_hit:
        print(result["best_one_hit"])
        return 0

    if args.check_rescreen:
        if args.min_one_hit is None:
            parser.error("--check-rescreen requires --min-one-hit")
        print(result["best_one_hit"])
        return 0 if result["best_one_hit"] >= args.min_one_hit else 1

    output_file = args.output_flag or args.best_species_out
    if not output_file:
        parser.error("Provide output path for best_species.txt")

    for msg in result["warnings"]:
        write_warning(warning_file, f"Warning: {msg}" if not msg.startswith("Warning") else msg)

    with open(output_file, "w", encoding="utf-8") as f_out:
        f_out.write(result["best_species"] + "\n")
    return 0


if __name__ == "__main__":
    if "snakemake" in globals():
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        sample = snakemake.wildcards.sample
        warning_file = f"logs/contamination_warnings/{sample}.txt"
        os.makedirs(os.path.dirname(warning_file), exist_ok=True)
        exclude_human = snakemake.params.get("exclude_human", True)

        try:
            species_data = load_species_data(input_file)
        except (OSError, ValueError) as e:
            msg = (
                f"FastQ Screen parsing issue for {sample} ({input_file}): {e}. "
                "Treating as 'No species found'."
            )
            write_warning(warning_file, msg)
            species_data = []

        result = resolve_best_species(species_data, exclude_human=exclude_human)
        for msg in result["warnings"]:
            prefixed = msg if msg.startswith("Warning") or msg.startswith("Human") else f"Warning: {msg}"
            write_warning(warning_file, prefixed.replace("Warning: Warning:", "Warning:"))

        with open(output_file, "w", encoding="utf-8") as f_out:
            f_out.write(result["best_species"] + "\n")
    else:
        raise SystemExit(main_cli())
