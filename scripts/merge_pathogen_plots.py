#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfMerger

# ------------------ Plot Functions ------------------


def plot_escore_reads(ax, escore_val, read_count, escore_thresh):
    ax.scatter(escore_val, read_count, color='blue')
    ax.axvline(escore_thresh, linestyle='--', color='green', label='Escore threshold')
    ax.set_xlabel("Escore")
    ax.set_ylabel("Read Count")
    ax.set_title("Krakenuniq: Escore vs Read Count")
    ax.legend()

def plot_breadth_evenness(ax, breadth, evenness):
    ax.bar(['Breadth', 'Evenness'], [breadth, evenness], color=['teal', 'orange'])
    ax.set_title("Breadth & Evenness of Coverage")
    ax.set_ylim(0, 1)

def plot_entropy(ax, entropy_100, entropy_1000):
    ax.bar(['Entropy_100bp', 'Entropy_1000bp'], [entropy_100, entropy_1000], color='purple')
    ax.set_title("Relative Entropy of Read Start Positions")
    ax.set_ylim(0, 1)

def plot_damage(ax, damage_df):
    if damage_df is not None and "Position" in damage_df.columns and "C>T" in damage_df.columns:
        ax.plot(damage_df["Position"], damage_df["C>T"], label='C>T')
        ax.set_title("5' Damage (C>T)")
        ax.set_xlabel("Position")
        ax.set_ylabel("Frequency")
    else:
        ax.set_title("Damage data unavailable")

def plot_length_dist(ax, length_df):
    if length_df is not None and "Length" in length_df.columns and "Occurrences" in length_df.columns:
        ax.bar(length_df["Length"], length_df["Occurrences"], color='gray')
        ax.set_title("Read Length Distribution")
        ax.set_xlabel("Length")
        ax.set_ylabel("Occurrences")
    else:
        ax.set_title("Length distribution unavailable")

def plot_ani(ax, ani_val):
    ax.text(0.1, 0.5, f"ANI: {ani_val}", fontsize=14)
    ax.axis('off')

# ------------------ Main Script ------------------

def main(args):
    if not os.path.isfile(args.summary_csv):
        print(f"Error: summary CSV file does not exist: {args.summary_csv}")
        sys.exit(1)

    try:
        summary_df = pd.read_csv(args.summary_csv)
    except pd.errors.EmptyDataError:
        print(f"Warning: summary CSV is empty: {args.summary_csv}. Skipping processing.")
        sys.exit(0)

    if summary_df.empty:
        print(f"Warning: summary CSV is empty: {args.summary_csv}. Skipping processing.")
        sys.exit(0)

    escore_thresh = args.escore_threshold
    

    base_pdf_path = args.output_pdf
    os.makedirs(os.path.dirname(base_pdf_path), exist_ok=True)

    pdf_files = []
    summary_records = []

    for _, row in summary_df.iterrows():
        pathogen = str(row["Pathogen"]).replace(" ", "_")
        sample = str(row["Sample"])

        pathogen_pdf = base_pdf_path.replace(".pdf", f"_{pathogen}_report.pdf")
        pdf_files.append(pathogen_pdf)

        with PdfPages(pathogen_pdf) as pdf_out:
            fig, axs = plt.subplots(3, 2, figsize=(12, 14))
            axs = axs.flatten()

            # Escore vs Read Count
            try:
                escore_val = float(row.get("Escore", 0))
                read_count = int(row.get("Krakenuniq_reads", 0))
                plot_escore_reads(axs[0], escore_val, read_count, escore_thresh)
            except Exception:
                axs[0].set_title("Escore/Read Count Plot Error")

            # Breadth & Evenness
            # Breadth & Evenness
            breadth = 0
            evenness = 0
            try:
                breadth = float(row.get("Breadth", 0))
                evenness = float(row.get("Evenness", 0))
                plot_breadth_evenness(axs[1], breadth, evenness)
            except Exception:
                axs[1].set_title("Breadth/Evenness Plot Error")# Relative Entropy
            try:
                entropy_100 = float(row.get("Entropy100", 0))
                entropy_1000 = float(row.get("Entropy1000", 0))
                plot_entropy(axs[2], entropy_100, entropy_1000)
            except Exception:
                axs[2].set_title("Entropy Plot Error")

            # 5' Damage plot - from DamagePlot.pdf or damageprofiler 5p_freq_misincorporations.txt fallback
            damage_plot_path = os.path.join(args.damage_dir, f"damageprofiler_{pathogen}", "DamagePlot.pdf")
            damage_txt_path = os.path.join(args.damage_dir, f"damageprofiler_{pathogen}", "5p_freq_misincorporations.txt")

            if os.path.exists(damage_plot_path):
                # We cannot embed PDF page directly in matplotlib, so mark info for later merging
                axs[3].text(0.5, 0.5, f"DamagePlot.pdf included separately\n({os.path.basename(damage_plot_path)})",
                            ha='center', va='center')
                axs[3].axis('off')
            elif os.path.exists(damage_txt_path):
                try:
                    damage_df = pd.read_csv(damage_txt_path, sep="\t", comment="#")
                    plot_damage(axs[3], damage_df)
                except Exception:
                    axs[3].set_title("Damage Plot Error")
            else:
                axs[3].set_title("Damage file missing")

            # Read Length Distribution
            lgfile = os.path.join(args.damage_dir, f"damageprofiler_{pathogen}", "lgdistribution.txt")
            if os.path.exists(lgfile):
                try:
                    length_df = pd.read_csv(lgfile, sep="\t", comment="#")
                    plot_length_dist(axs[4], length_df)
                except Exception:
                    axs[4].set_title("Read Length Plot Error")
            else:
                axs[4].set_title("Length file missing")

            # ANI plot
            # ANI plot
            ani_val = "NA"
            ani_file = os.path.join(args.damage_dir, f"{sample}_{pathogen}.ani.txt")
            if os.path.exists(ani_file):
                try:
                    with open(ani_file, "r") as f:
                        line = f.readline().strip()
                        import re
                        match = re.search(r"(\d+\.?\d*)%", line)
                        if match:
                            ani_val = float(match.group(1))
                except Exception:
                    pass

            plot_ani(axs[5], ani_val)

            plt.suptitle(f"Sample: {sample} | Pathogen: {row['Pathogen']}", fontsize=16)
            pdf_out.savefig(fig)

        # Append extra PDFs for merging later
        extra_pdfs = []

        # aDNAPlotter.pdf located directly in bwa_pathogen directory
        aDNA_pdf_path = os.path.join(args.damage_dir, "aDNAPlotter.pdf")
        if os.path.exists(aDNA_pdf_path):
            extra_pdfs.append(aDNA_pdf_path)

        # DamagePlot.pdf from damageprofiler folder
        if os.path.exists(damage_plot_path):
            extra_pdfs.append(damage_plot_path)

        pdf_files.extend(extra_pdfs)

        # Collect summary info for CSV output
        try:
            mean_read_len = None
            if os.path.exists(lgfile):
                length_df = pd.read_csv(lgfile, sep="\t", comment="#")
                # Compute mean read length
                if all(x in length_df.columns for x in ["Length", "Occurrences"]):
                    mean_read_len = np.average(length_df["Length"], weights=length_df["Occurrences"])
        except Exception:
            mean_read_len = None

        summary_records.append({
            "Sample": sample,
            "Pathogen": row["Pathogen"],
            "ANI": ani_val,
            "Breadth": breadth,
            "Evenness": evenness,
            "Krakenuniq_reads": row.get("Krakenuniq_reads", None),
            "Mean_read_length": mean_read_len
        })

    # Merge all PDFs
    if args.merge_pdfs and pdf_files:
        merged_pdf = base_pdf_path.replace(".pdf", "_merged.pdf")
        merger = PdfMerger()
        for pdf_file in pdf_files:
            if os.path.exists(pdf_file):
                merger.append(pdf_file)
            else:
                print(f"Warning: PDF file missing, skipping merge: {pdf_file}")
        merger.write(merged_pdf)
        merger.close()
        print(f"Merged PDF saved to: {merged_pdf}")

    # Write summary CSV file
    summary_out_csv = base_pdf_path.replace(".pdf", "_summary.csv")
    summary_df_out = pd.DataFrame(summary_records)
    summary_df_out.to_csv(summary_out_csv, index=False)
    print(f"Summary CSV saved to: {summary_out_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary_csv", required=True)
    parser.add_argument("--damage_dir", required=True)
    parser.add_argument("--output_pdf", required=True)
    parser.add_argument("--escore_threshold", type=float, default=0.05)
    parser.add_argument("--merge_pdfs", action="store_true")
    args = parser.parse_args()
    main(args)
