import pandas as pd
import sys
import traceback
import seaborn as sns
import matplotlib.pyplot as plt
import os

def normalize_name(name):
    return str(name).lower().replace("_", " ").strip()

def main():
    escore_file = sys.argv[1]
    heatmap_file = sys.argv[2]
    sample_name = sys.argv[3]
    spreadsheet_path = sys.argv[4]

    # Load KrakenUniq E-score output (with lowercase headers)
    df_escore = pd.read_csv(escore_file)
    df_escore.columns = [col.strip().lower() for col in df_escore.columns]
    kraken_set = set(df_escore["taxonomy"].map(normalize_name))

    # Load HOPS heatmap
    df_hops = pd.read_csv(heatmap_file, sep="\t")
    df_hops.columns = [c.strip().strip('"') for c in df_hops.columns]
    df_hops.set_index(df_hops.columns[0], inplace=True)

    sample_col = f"{sample_name}_unaligned.rma6"
    if sample_col not in df_hops.columns:
        print(f"❌ Sample column '{sample_col}' not found in HOPS heatmap!")
        print("✔️ Available columns:", ", ".join(df_hops.columns))
        print(f"⚠️ Skipping sample {sample_name}. No output will be written.")
        
        html_output = f"results/comparison/{sample_name}_heatmap.html"
        tsv_output = f"results/comparison/{sample_name}_comparison.tsv"
        os.makedirs("results/comparison", exist_ok=True)
        with open(html_output, "w") as f:
            f.write(f"<h3>Sample column '{sample_col}' not found. Skipping.</h3>\n")
        with open(tsv_output, "w") as f:
            f.write("taxid\tKrakenuniq name\tIn_Krakenuniq\tIn_HOPS\n")
        return

    hops_present = set(
        df_hops[sample_col][df_hops[sample_col] >= 2]
        .index.map(normalize_name)
    )

    # Load pathogen spreadsheet
    df_sheet = pd.read_csv(spreadsheet_path)
    col_map = {col.lower().strip(): col for col in df_sheet.columns}
    taxid_col = col_map.get("taxid")
    kraken_col = col_map.get("krakenuniq name")
    hops_col = col_map.get("hops name")

    if not taxid_col or not kraken_col or not hops_col:
        raise KeyError(f"Missing required columns. Found: {list(df_sheet.columns)}")

    df_sheet["krakenuniq_norm"] = df_sheet[kraken_col].map(normalize_name)
    df_sheet["hops_norm"] = df_sheet[hops_col].map(normalize_name)

    df_merged = pd.DataFrame({
        "taxid": df_sheet[taxid_col],
        "Krakenuniq name": df_sheet[kraken_col],
        "In_Krakenuniq": df_sheet["krakenuniq_norm"].isin(kraken_set),
        "In_HOPS": df_sheet["hops_norm"].isin(hops_present)
    })

    df_merged = df_merged[(df_merged["In_Krakenuniq"]) | (df_merged["In_HOPS"])]

    os.makedirs("results/comparison", exist_ok=True)
    tsv_output = f"results/comparison/{sample_name}_comparison.tsv"
    df_merged.to_csv(tsv_output, sep="\t", index=False)

    if df_merged.empty:
        print(f"ℹ️ No pathogens found for sample {sample_name}, skipping heatmap.")
        html_output = f"results/comparison/{sample_name}_heatmap.html"
        with open(html_output, "w") as f:
            f.write(f"<h3>No pathogens found in {sample_name}, so no heatmap was generated.</h3>\n")
        return

    heatmap_data = df_merged.set_index("Krakenuniq name")[["In_Krakenuniq", "In_HOPS"]].astype(int)
    plt.figure(figsize=(5, max(2, 0.4 * len(heatmap_data))))
    sns.heatmap(
        heatmap_data,
        cmap=["#FFFFFF", "#3b82f6"],
        linewidths=0.5,
        linecolor="gray",
        cbar=False,
        annot=True,
        fmt="d"
    )
    plt.title(f"Presence of Pathogens — {sample_name}")
    plt.tight_layout()

    img_path = f"results/comparison/{sample_name}_heatmap.png"
    html_output = f"results/comparison/{sample_name}_heatmap.html"
    plt.savefig(img_path, dpi=150)
    plt.close()

    with open(html_output, "w") as f:
        f.write(f"<h2>Sample: {sample_name}</h2>\n")
        f.write(f"<img src='{os.path.basename(img_path)}'>\n")

if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
