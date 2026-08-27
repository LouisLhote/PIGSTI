#!/usr/bin/env python
"""
Calculate edit distance distribution and log(R²) from pathogen BAM file.
The edit distance (NM tag) distribution should follow a pattern where
log(R²) close to 1 indicates good alignment quality.
"""

import pysam
import numpy as np
import sys
from scipy import stats
import os


def _named_output(name, default=None):
    """Resolve an optional named Snakemake output without raising."""
    try:
        out = snakemake.output
    except NameError:
        return default
    try:
        keys = list(out.keys())
        if name in keys:
            return out[name]
    except Exception:
        pass
    try:
        return getattr(out, name)
    except (AttributeError, KeyError):
        return default


def _ensure_parent(path):
    if path:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)


def _write_empty_ani(path):
    if not path:
        return
    _ensure_parent(path)
    open(path, "w").close()


def _write_ed_distribution(path, full_counts=None):
    """Always write edit-distance bins 0-5 (zeros if no data)."""
    if not path:
        return
    _ensure_parent(path)
    if full_counts is None:
        full_counts = np.zeros(6, dtype=float)
    with open(path, "w") as f:
        for i in range(6):
            count = int(full_counts[i]) if i < len(full_counts) else 0
            f.write(f"{i}\t{count}\n")


def _write_empty_plot(path, message="No edit distance data available"):
    _ensure_parent(path)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=14)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


# Get inputs from Snakemake
try:
    bam = snakemake.input.bam
    output_r2 = snakemake.output.r2
    output_plot = snakemake.output.plot
    ani_dist_file = _named_output("ani_distribution")
    ed_dist_file = _named_output("edit_distance_distribution")
    # Fallbacks only for the main EditDistanceR2 BAM (*.dedup.bam), where these
    # files are rule outputs. Damage/default subset rules declare only r2+plot.
    if ".dedup.bam" in bam:
        if ani_dist_file is None:
            ani_dist_file = bam.replace(".dedup.bam", ".ani_distribution.txt")
        if ed_dist_file is None:
            ed_dist_file = (
                output_r2[: -len(".txt")] + "_distribution.txt"
                if output_r2.endswith(".txt")
                else output_r2 + "_distribution.txt"
            )
    metric_start_ed = 0
    try:
        metric_start_ed = int(getattr(snakemake.params, "metric_start_ed", 0))
    except Exception:
        metric_start_ed = 0
except (AttributeError, NameError):
    # Fallback for command-line usage
    import argparse

    parser = argparse.ArgumentParser(description="Calculate edit distance R² from BAM file")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--output-r2", required=True, help="Output R² text file")
    parser.add_argument("--output-plot", required=True, help="Output plot PNG file")
    parser.add_argument(
        "--ani-distribution",
        required=False,
        default=None,
        help="Optional per-read ANI distribution output",
    )
    parser.add_argument(
        "--edit-distance-distribution",
        required=False,
        default=None,
        help="Optional edit-distance count distribution output",
    )
    parser.add_argument(
        "--metric-start-ed",
        required=False,
        type=int,
        default=0,
        help="Edit-distance bin to start the decay-quality metric from (default: 0)",
    )
    args = parser.parse_args()
    bam = args.bam
    output_r2 = args.output_r2
    output_plot = args.output_plot
    metric_start_ed = int(getattr(args, "metric_start_ed", 0))
    ani_dist_file = args.ani_distribution
    ed_dist_file = args.edit_distance_distribution
    if ".dedup.bam" in bam:
        if ani_dist_file is None:
            ani_dist_file = bam.replace(".dedup.bam", ".ani_distribution.txt")
        if ed_dist_file is None:
            ed_dist_file = (
                output_r2[: -len(".txt")] + "_distribution.txt"
                if output_r2.endswith(".txt")
                else output_r2 + "_distribution.txt"
            )

# Touch required distribution outputs early so Snakemake never sees a silent miss
# if a later step fails after writing r2/plot.
_write_empty_ani(ani_dist_file)
_write_ed_distribution(ed_dist_file)

# Check if BAM index exists
bam_index = bam + ".bai"
if not os.path.exists(bam_index):
    sys.stderr.write(f"WARNING: Missing BAM index file: {bam_index}. Creating index...\n")
    pysam.index(bam)

# Open BAM file
bamfile = pysam.AlignmentFile(bam, "rb")

# Extract edit distances (NM tag) and per-read ANI from all mapped reads
edit_distances = []
ani_values = []
total_reads = 0
reads_with_nm = 0

for read in bamfile.fetch():
    if not read.is_unmapped:
        total_reads += 1
        try:
            nm = read.get_tag("NM")
        except KeyError:
            # Some reads might not have NM tag
            continue

        edit_distances.append(nm)
        reads_with_nm += 1

        # Approximate per-read ANI from NM and read length:
        # ANI (%) ≈ 100 * (1 - NM / read_length)
        read_len = read.query_length
        if read_len and read_len > 0:
            ani = 100.0 * (1.0 - float(nm) / float(read_len))
            # Clamp to [0, 100] to avoid numerical oddities
            ani = max(0.0, min(100.0, ani))
            ani_values.append(ani)

bamfile.close()

# If no reads with NM tag, return default values (still write every declared output).
if reads_with_nm == 0 or len(edit_distances) == 0:
    sys.stderr.write(f"WARNING: No reads with NM tag found in {bam}. Writing default values.\n")
    _ensure_parent(output_r2)
    with open(output_r2, "w") as f:
        f.write("NA\n")
    _write_empty_ani(ani_dist_file)
    _write_ed_distribution(ed_dist_file)
    _write_empty_plot(output_plot)
    sys.exit(0)

edit_distances = np.array(edit_distances)

# Calculate distribution (counts per edit distance value)
unique_ed, counts = np.unique(edit_distances, return_counts=True)
max_ed = min(int(max(unique_ed)), 20)  # Cap at 20 for visualization

# --- Regression counts: EXACTLY as in your example ---
# We build a count vector for all distances 0..max_ed (including zeros)
full_counts = np.zeros(max_ed + 1, dtype=float)
for ed, count in zip(unique_ed, counts):
    if 0 <= ed <= max_ed:
        full_counts[int(ed)] = count

# Initialize metrics
slope_linear = 0.0
intercept_linear = 0.0
r2_linear = 0.0
r2_log = 0.0
edit_r2_value = 0.0
decay_quality_score = 0.0  # New composite metric

# Prepare x and y vectors. Optionally skip edit-distance bins starting at metric_start_ed
# (e.g., for damage subsets exclude edit distance 0).
if metric_start_ed and metric_start_ed > 0 and metric_start_ed < len(full_counts):
    y_vals = full_counts[metric_start_ed:].astype(float)
    x_vals = np.arange(len(y_vals), dtype=float)
else:
    x_vals = np.arange(len(full_counts), dtype=float)
    y_vals = full_counts.astype(float)

if len(x_vals) >= 2:
    # Linear regression: counts vs edit distance
    slope_linear, intercept_linear, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
    r2_linear = r_value**2

    # Log-linear regression: log1p(counts) vs edit distance
    log_y = np.log1p(y_vals)  # log(y + 1) to avoid log(0)
    slope_log, intercept_log, r_log_value, _, _ = stats.linregress(x_vals, log_y)
    r2_log = r_log_value**2

    # --- NEW: Decay Quality Score (composite metric for descending pattern) ---
    # This metric combines multiple signals to detect good descending patterns
    # like 0:500, 1:350, 2:283 (good) vs 0:150, 1:350, 2:120 (bad)

    # 1. Monotonicity check: Is the distribution strictly decreasing?
    # Check if count[i] >= count[i+1] for all i (allowing for small noise)
    is_monotonic = True
    tolerance = 0.05  # Allow 5% tolerance for noise
    for i in range(len(y_vals) - 1):
        if y_vals[i] > 0 and y_vals[i + 1] > 0:
            # Check if next value is significantly larger (bad)
            if y_vals[i + 1] > y_vals[i] * (1 + tolerance):
                is_monotonic = False
                break

    monotonicity_score = 1.0 if is_monotonic else 0.0

    # 2. Dominance ratio: How much of the signal is at 0-mismatch?
    # Higher ratio = better (more reads with perfect matches)
    total_reads = np.sum(y_vals)
    count_0 = y_vals[0] if len(y_vals) > 0 else 0
    count_1 = y_vals[1] if len(y_vals) > 1 else 0

    if total_reads > 0:
        # Ratio of 0-mismatch to 1-mismatch (higher is better)
        if count_1 > 0:
            dominance_ratio = count_0 / count_1
        else:
            dominance_ratio = 10.0  # Perfect: all reads at 0-mismatch
        # Normalize to 0-1 scale (assuming good cases have ratio > 1.0)
        # Typical good: ratio ~1.4-2.0, perfect: ratio >> 1
        dominance_score = min(1.0, dominance_ratio / 2.0)  # Cap at 1.0 for ratio >= 2.0

        # Also compute fraction of reads at 0-mismatch
        fraction_0 = count_0 / total_reads
        fraction_0_score = fraction_0  # Direct: higher is better
    else:
        dominance_score = 0.0
        fraction_0_score = 0.0

    # 3. Peak detection: Penalize if peak is not at 0
    # Find the position with maximum count
    max_idx = np.argmax(y_vals)
    peak_penalty = 0.0 if max_idx == 0 else 0.5  # Heavy penalty if peak not at 0

    # 4. Exponential decay fit quality
    # Fit exponential: count(ed) = A * exp(-λ * ed)
    # Good patterns should fit well to exponential decay
    try:
        # Only fit on non-zero values
        non_zero_mask = y_vals > 0
        if np.sum(non_zero_mask) >= 2:
            x_fit = x_vals[non_zero_mask]
            y_fit = y_vals[non_zero_mask]
            log_y_fit = np.log(y_fit + 1e-10)  # log(y) for exponential fit

            # Fit: log(y) = log(A) - λ*x
            slope_exp, intercept_exp, r_exp, _, _ = stats.linregress(x_fit, log_y_fit)
            r2_exp = r_exp**2

            # Decay rate (lambda): more negative is better (faster decay)
            decay_rate = -slope_exp
            # Normalize decay rate score (typical good: 0.3-0.7, perfect: >1.0)
            decay_rate_score = min(1.0, decay_rate / 1.0)

            # R² of exponential fit (higher is better)
            r2_exp_score = r2_exp
        else:
            decay_rate_score = 0.0
            r2_exp_score = 0.0
    except Exception:
        decay_rate_score = 0.0
        r2_exp_score = 0.0

    # 5. Slope-based score (from linear regression)
    # More negative slope = better (faster decline)
    if slope_linear < 0:
        # Normalize: typical good slopes are -50 to -200, perfect: < -200
        slope_score = min(1.0, abs(slope_linear) / 200.0)
    else:
        slope_score = 0.0  # Positive slope is very bad

    # Combine all components into composite score (weighted average)
    # Weights emphasize monotonicity and dominance
    weights = {
        "monotonicity": 0.35,  # Most important: must be descending
        "dominance": 0.25,  # High 0-mismatch ratio
        "fraction_0": 0.15,  # Fraction at 0-mismatch
        "peak_penalty": -0.10,  # Penalty if peak not at 0
        "decay_rate": 0.10,  # Exponential decay rate
        "r2_exp": 0.05,  # Quality of exponential fit
    }

    decay_quality_score = (
        weights["monotonicity"] * monotonicity_score
        + weights["dominance"] * dominance_score
        + weights["fraction_0"] * fraction_0_score
        + weights["peak_penalty"] * peak_penalty
        + weights["decay_rate"] * decay_rate_score
        + weights["r2_exp"] * r2_exp_score
    )

    # Ensure score is in [0, 1] range
    decay_quality_score = max(0.0, min(1.0, decay_quality_score))

else:
    sys.stderr.write(
        "WARNING: Insufficient data points for regression. Setting edit distance metric = 0.0\n"
    )

# Choose metric for thresholding: use new decay_quality_score as primary
# Keep old edit_r2_value for backward compatibility
base_r2_for_metric = r2_log if r2_log > 0 else r2_linear
abs_r2_for_metric = abs(base_r2_for_metric)
log_abs_r2 = np.log10(abs_r2_for_metric + 1e-10)  # log10(|R²|)
edit_r2_value = max(0.0, min(1.0, 1.0 + log_abs_r2))

# Write metrics to text file for downstream use (R plotting and summaries)
_ensure_parent(output_r2)
with open(output_r2, "w") as f:
    # First line: NEW primary metric - Decay Quality Score (0-1 scale, higher is better)
    f.write(f"{decay_quality_score:.4f}\n")
    # Second line: Legacy transformed metric (for backward compatibility)
    f.write(f"Legacy_metric={edit_r2_value:.4f}\n")
    # Core regression metrics (matching your regression.py output)
    f.write(f"Slope_linear={slope_linear:.4f}\n")
    f.write(f"Decay_quality_score={decay_quality_score:.4f}\n")

# Per-read ANI distribution (required by EditDistanceR2 rule when declared)
if ani_dist_file:
    _ensure_parent(ani_dist_file)
    with open(ani_dist_file, "w") as f_ani:
        for v in ani_values:
            f_ani.write(f"{v:.4f}\n")

# Edit distance distribution bins 0-5 (required by EditDistanceR2 rule when declared)
_write_ed_distribution(ed_dist_file, full_counts)

# Create plot
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Edit distance distribution (0..max_ed for context)
ed_range_full = np.arange(0, max_ed + 1)
ax1.bar(ed_range_full, full_counts, alpha=0.7, color="steelblue", edgecolor="black")
ax1.set_xlabel("Edit Distance (NM)", fontsize=12, fontweight="bold")
ax1.set_ylabel("Read Count", fontsize=12, fontweight="bold")
ax1.set_title("Edit Distance Distribution", fontsize=14, fontweight="bold")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-0.5, max_ed + 0.5)

# Plot 2: Occurrences vs Edit Distance (for R²_linear visualization)
if len(x_vals) >= 2:
    ax2.scatter(
        x_vals,
        y_vals,
        s=50,
        alpha=0.7,
        color="crimson",
        edgecolor="black",
    )
    # Plot linear regression line
    fitted_counts = slope_linear * x_vals + intercept_linear
    ax2.plot(
        x_vals,
        fitted_counts,
        "r--",
        linewidth=2,
        label="Linear fit",
    )
    ax2.legend(fontsize=10)
    ax2.set_xlabel("Edit Distance (NM)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Read Count", fontsize=12, fontweight="bold")
    ax2.set_title(
        f"Occurrences vs Edit Distance\nDecay quality score = {decay_quality_score:.4f}",
        fontsize=14,
        fontweight="bold",
    )
    ax2.grid(True, alpha=0.3)
else:
    ax2.text(
        0.5,
        0.5,
        "Insufficient data for regression",
        ha="center",
        va="center",
        fontsize=12,
        transform=ax2.transAxes,
    )
    ax2.set_xlabel("Edit Distance (NM)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Read Count", fontsize=12, fontweight="bold")
    ax2.set_title(
        "Occurrences vs Edit Distance\nInsufficient data",
        fontsize=14,
        fontweight="bold",
    )
    ax2.grid(True, alpha=0.3)

plt.tight_layout()
_ensure_parent(output_plot)
plt.savefig(output_plot, dpi=300, bbox_inches="tight")
plt.close()

print(f"Calculated edit distance metrics for {bam}:")
print(f"   Decay quality score = {decay_quality_score:.4f} (0-1, higher is better)")
