import sys
import os

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
        "Genome", "#Unmapped", "#One_hit_one_genome", "%Unmapped",
        "%One_hit_one_genome", "%Multiple_hits_one_genome"
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
            data.append({
                "genome": genome,
                "unmapped": unmapped,
                "one_hit": one_hit,
                "one_hit_pct": one_hit_pct,
                "multi_hit_pct": multi_hit_pct
            })
        except (IndexError, ValueError):
            continue
    return data

def write_warning(path, message):
    with open(path, "a") as warn:
        warn.write(message + "\n")

if __name__ == "__main__":
    # Detect Snakemake or CLI mode
    if "snakemake" in globals():
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        sample = snakemake.wildcards.sample
        warning_file = f"logs/contamination_warnings/{sample}.txt"
        os.makedirs(os.path.dirname(warning_file), exist_ok=True)
        exclude_human = snakemake.params.get("exclude_human", True)
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        sample = "manual_run"
        warning_file = "warning.txt"
        exclude_human = True

    with open(input_file) as f_in:
        species_data = parse_headers_and_data(f_in)

    if not species_data:
        raise ValueError(f"No valid species entries in input file: {input_file}")

    # Handle human contamination logging
    human_entry = next((s for s in species_data if s["genome"].lower() == "human"), None)
    if human_entry:
        one_hit_pct = human_entry.get("one_hit_pct")
        multi_hit_pct = human_entry.get("multi_hit_pct")
        if one_hit_pct is not None and multi_hit_pct is not None:
            # Corrected calculation:
            human_contam = one_hit_pct + multi_hit_pct
            write_warning(warning_file, f"Human contamination load for {sample}: {human_contam:.2f}%")
        else:
            write_warning(warning_file, f"Warning: Could not compute human contamination load for {sample}")

    # Filter out human if needed
    filtered_data = [
        s for s in species_data if not (exclude_human and s["genome"].lower() == "human")
    ]

    if not filtered_data:
        best_species = "No species found"
    else:
        best_hit = max(filtered_data, key=lambda x: x["one_hit"])
        best_unmapped = min(filtered_data, key=lambda x: x["unmapped"])

        if best_hit["genome"] != best_unmapped["genome"]:
            best_species = best_unmapped["genome"]
            write_warning(
                warning_file,
                f"Warning: Best one-hit ({best_hit['genome']}) and lowest unmapped ({best_unmapped['genome']}) do not match. Using {best_species}."
            )
        else:
            best_species = best_hit["genome"]

        # Check for contamination (best one_hit vs second best)
        non_best = [s for s in filtered_data if s["genome"].lower() != best_hit["genome"].lower()]
        if non_best:
            second_best = max(non_best, key=lambda x: x["one_hit"])
            if second_best["one_hit"] > 0 and (best_hit["one_hit"] / second_best["one_hit"]) < 1.5:
                write_warning(
                    warning_file,
                    f"Warning: Probable {second_best['genome']} contamination in sample {sample}. "
                    f"{best_hit['genome']} only {best_hit['one_hit']} vs {second_best['genome']} {second_best['one_hit']}."
                )

    # Write result
    with open(output_file, "w") as f_out:
        f_out.write(best_species + "\n")
