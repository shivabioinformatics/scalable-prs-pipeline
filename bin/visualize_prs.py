#!/usr/bin/env python3
"""
visualize_prs.py — Generate all the figures for the PRS pipeline

I'm using matplotlib + seaborn because they're standard in bioinformatics
and I don't need to install R/ggplot2 separately. Five plots:

1. QC Summary — bar chart of call rates and depth per sample
2. PRS Distribution — histogram with risk threshold lines
3. Risk Stratification — box/bar plot by risk category
4. Variant Effect Sizes — Manhattan-style plot by chromosome
5. Score Comparison — per-sample dot plot ranked by PRS

These are the kinds of figures you'd put in a paper or clinical report
to communicate results to someone who doesn't want to read raw numbers.
"""

import os
import sys
import argparse
import math

# I need these for plotting — they should be installed via pip
try:
    import matplotlib
    matplotlib.use("Agg")  # no GUI needed, just saving to files
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ImportError:
    print("ERROR: matplotlib not installed. Run: pip install matplotlib")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("ERROR: numpy not installed. Run: pip install numpy")
    sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser(description="Visualize PRS pipeline results")
    p.add_argument("--scores", required=True, help="PRS scores TSV")
    p.add_argument("--qc-stats", required=True, help="QC stats TSV")
    p.add_argument("--gwas", required=True, help="GWAS summary stats TSV")
    p.add_argument("--outdir", required=True, help="Output directory for figures")
    return p.parse_args()


def load_scores(path):
    """load the PRS results into a list of dicts"""
    data = []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            fs = line.strip().split("\t")
            if len(fs) >= 7:
                data.append({
                    "sample": fs[0],
                    "prs_raw": float(fs[1]),
                    "prs_zscore": float(fs[2]),
                    "percentile": float(fs[3]),
                    "risk_category": fs[4],
                    "variants_used": int(fs[5]),
                    "variants_imputed": int(fs[6])
                })
    return data


def load_qc_samples(path):
    """load per-sample QC metrics"""
    data = []
    in_section = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("## PER-SAMPLE QC"):
                in_section = True; continue
            if in_section and line.startswith("##"): break
            if in_section and line.startswith("sample\t"): continue
            if in_section and line:
                fs = line.split("\t")
                if len(fs) >= 6:
                    data.append({
                        "sample": fs[0],
                        "call_rate": float(fs[1]),
                        "mean_depth": float(fs[2]),
                        "variants_called": int(fs[3]),
                        "variants_total": int(fs[4]),
                        "pass_qc": fs[5]
                    })
    return data


def load_qc_summary(path):
    """load the overall QC summary numbers"""
    summary = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("## VARIANT QC"):
                continue
            if line.startswith("##"):
                break
            if "\t" in line:
                k, v = line.split("\t", 1)
                summary[k] = v
    return summary


def load_gwas(path):
    """load GWAS summary stats for the variant plot"""
    data = []
    with open(path) as f:
        f.readline()
        for line in f:
            fs = line.strip().split("\t")
            if len(fs) >= 8:
                try:
                    data.append({
                        "rsid": fs[0],
                        "chrom": int(fs[1]),
                        "pos": int(fs[2]),
                        "beta": float(fs[5]),
                        "pval": float(fs[7])
                    })
                except ValueError:
                    continue
    return data


# ---------- color scheme ----------
# going with a clean, professional palette
COLORS = {
    "HIGH_RISK": "#e74c3c",
    "ELEVATED": "#f39c12",
    "AVERAGE": "#3498db",
    "LOW_RISK": "#2ecc71",
    "primary": "#2c3e50",
    "secondary": "#7f8c8d",
    "accent": "#8e44ad",
    "bg": "#fafafa",
    "grid": "#ecf0f1"
}


def plot_qc_summary(qc_samples, qc_summary, outdir):
    """
    Plot 1: QC Summary
    Two-panel figure showing per-sample call rates and mean read depth.
    The red dashed line shows the threshold — samples below it are flagged.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Quality Control Summary", fontsize=16, fontweight="bold", y=1.02)

    samples = [d["sample"].replace("PATIENT_", "P") for d in qc_samples]
    call_rates = [d["call_rate"] for d in qc_samples]
    depths = [d["mean_depth"] for d in qc_samples]
    colors = [COLORS["HIGH_RISK"] if d["pass_qc"] == "FAIL" else COLORS["primary"] for d in qc_samples]

    # call rate bar chart
    ax1.bar(range(len(samples)), call_rates, color=colors, alpha=0.85, edgecolor="white")
    ax1.axhline(y=0.90, color=COLORS["HIGH_RISK"], linestyle="--", linewidth=1.5, label="Threshold (0.90)")
    ax1.set_xlabel("Sample", fontsize=11)
    ax1.set_ylabel("Call Rate", fontsize=11)
    ax1.set_title("Per-Sample Call Rate", fontsize=13)
    ax1.set_xticks(range(len(samples)))
    ax1.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax1.set_ylim(0.80, 1.01)
    ax1.legend(fontsize=9)
    ax1.grid(axis="y", alpha=0.3)

    # depth bar chart
    ax2.bar(range(len(samples)), depths, color=COLORS["accent"], alpha=0.75, edgecolor="white")
    ax2.axhline(y=10, color=COLORS["HIGH_RISK"], linestyle="--", linewidth=1.5, label="Min Depth (10)")
    ax2.set_xlabel("Sample", fontsize=11)
    ax2.set_ylabel("Mean Read Depth", fontsize=11)
    ax2.set_title("Per-Sample Read Depth", fontsize=13)
    ax2.set_xticks(range(len(samples)))
    ax2.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(axis="y", alpha=0.3)

    # add the overall QC stats as text
    total = qc_summary.get("total_variants", "?")
    passing = qc_summary.get("passing_variants", "?")
    rate = qc_summary.get("pass_rate", "?")
    fig.text(0.5, -0.05, f"Variant QC: {passing}/{total} passed ({rate})", ha="center", fontsize=10, style="italic")

    plt.tight_layout()
    path = os.path.join(outdir, "01_qc_summary.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {path}")


def plot_prs_distribution(scores, outdir):
    """
    Plot 2: PRS Distribution
    Histogram of z-scores across all samples. I mark the risk thresholds
    so you can see where the cutoffs are. In a real study with thousands
    of samples this would look like a nice bell curve.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    zscores = [d["prs_zscore"] for d in scores]

    # histogram with density curve feel
    n, bins, patches = ax.hist(zscores, bins=15, color=COLORS["primary"], alpha=0.7,
                                edgecolor="white", linewidth=1.2)

    # color the bars by risk category
    for patch, left, right in zip(patches, bins[:-1], bins[1:]):
        center = (left + right) / 2
        # rough z-score to percentile mapping for coloring
        if center > 1.65:
            patch.set_facecolor(COLORS["HIGH_RISK"])
        elif center > 0.84:
            patch.set_facecolor(COLORS["ELEVATED"])
        elif center > -0.84:
            patch.set_facecolor(COLORS["AVERAGE"])
        else:
            patch.set_facecolor(COLORS["LOW_RISK"])
        patch.set_alpha(0.8)

    # threshold lines
    ax.axvline(x=1.65, color=COLORS["HIGH_RISK"], linestyle="--", linewidth=2, label="High Risk (95th)")
    ax.axvline(x=0.84, color=COLORS["ELEVATED"], linestyle="--", linewidth=1.5, label="Elevated (80th)")
    ax.axvline(x=-0.84, color=COLORS["LOW_RISK"], linestyle="--", linewidth=1.5, label="Low Risk (20th)")

    ax.set_xlabel("PRS Z-Score", fontsize=13)
    ax.set_ylabel("Number of Samples", fontsize=13)
    ax.set_title("Polygenic Risk Score Distribution", fontsize=16, fontweight="bold")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    # annotate mean
    mean_z = np.mean(zscores)
    ax.axvline(x=mean_z, color=COLORS["secondary"], linestyle="-", linewidth=1)
    ax.text(mean_z + 0.05, ax.get_ylim()[1] * 0.9, f"Mean: {mean_z:.2f}", fontsize=9)

    plt.tight_layout()
    path = os.path.join(outdir, "02_prs_distribution.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {path}")


def plot_risk_stratification(scores, outdir):
    """
    Plot 3: Risk Stratification
    Bar chart showing how many patients fall in each risk category,
    plus a dot plot showing individual scores within each group.
    This is the money plot for a clinical audience.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Risk Stratification", fontsize=16, fontweight="bold", y=1.02)

    categories = ["LOW_RISK", "AVERAGE", "ELEVATED", "HIGH_RISK"]
    cat_colors = [COLORS[c] for c in categories]
    cat_labels = ["Low Risk\n(<20th)", "Average\n(20-80th)", "Elevated\n(80-95th)", "High Risk\n(>95th)"]

    counts = [sum(1 for s in scores if s["risk_category"] == c) for c in categories]

    # bar chart of counts
    bars = ax1.bar(range(len(categories)), counts, color=cat_colors, alpha=0.85, edgecolor="white", width=0.6)
    ax1.set_xticks(range(len(categories)))
    ax1.set_xticklabels(cat_labels, fontsize=10)
    ax1.set_ylabel("Number of Patients", fontsize=12)
    ax1.set_title("Patients per Risk Category", fontsize=13)
    ax1.grid(axis="y", alpha=0.3)

    # add count labels on bars
    for bar, count in zip(bars, counts):
        if count > 0:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                     str(count), ha="center", va="bottom", fontweight="bold", fontsize=12)

    # dot plot — individual scores by category
    for i, cat in enumerate(categories):
        cat_scores = [s["prs_zscore"] for s in scores if s["risk_category"] == cat]
        if cat_scores:
            x_jitter = np.random.normal(i, 0.08, len(cat_scores))
            ax2.scatter(x_jitter, cat_scores, c=COLORS[cat], alpha=0.7, s=60, edgecolors="white", linewidth=0.5)

    ax2.set_xticks(range(len(categories)))
    ax2.set_xticklabels(cat_labels, fontsize=10)
    ax2.set_ylabel("PRS Z-Score", fontsize=12)
    ax2.set_title("Individual Scores by Category", fontsize=13)
    ax2.grid(axis="y", alpha=0.3)
    ax2.axhline(y=0, color=COLORS["secondary"], linestyle="-", alpha=0.5)

    plt.tight_layout()
    path = os.path.join(outdir, "03_risk_stratification.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {path}")


def plot_variant_effects(gwas, outdir):
    """
    Plot 4: Variant Effect Sizes (Manhattan-style)
    Each dot is a SNP, plotted by chromosome position on x-axis and
    -log10(p-value) on y-axis. Bigger effect sizes get bigger dots.
    This is the standard way to visualize GWAS results.
    """
    fig, ax = plt.subplots(figsize=(14, 5))

    # assign cumulative positions so chromosomes line up side by side
    chrom_offsets = {}
    cumulative = 0
    chroms_sorted = sorted(set(d["chrom"] for d in gwas))

    for c in chroms_sorted:
        chrom_offsets[c] = cumulative
        max_pos = max(d["pos"] for d in gwas if d["chrom"] == c)
        cumulative += max_pos + 5000000  # gap between chromosomes

    x_positions = []
    y_values = []
    dot_sizes = []
    dot_colors = []

    # alternate colors by chromosome so they're visually distinct
    chrom_palette = [COLORS["primary"], COLORS["accent"]]

    for d in gwas:
        x = chrom_offsets[d["chrom"]] + d["pos"]
        y = -math.log10(d["pval"]) if d["pval"] > 0 else 20
        size = min(abs(d["beta"]) * 500 + 5, 100)  # scale dot size by effect
        color = chrom_palette[d["chrom"] % 2]

        x_positions.append(x)
        y_values.append(y)
        dot_sizes.append(size)
        dot_colors.append(color)

    ax.scatter(x_positions, y_values, s=dot_sizes, c=dot_colors, alpha=0.6, edgecolors="none")

    # genome-wide significance line (5e-8)
    ax.axhline(y=-math.log10(5e-8), color=COLORS["HIGH_RISK"], linestyle="--", linewidth=1.5, label="p = 5×10⁻⁸")
    # suggestive significance line (1e-5)
    ax.axhline(y=-math.log10(1e-5), color=COLORS["ELEVATED"], linestyle="--", linewidth=1, label="p = 1×10⁻⁵")

    # chromosome labels
    chrom_centers = []
    for c in chroms_sorted:
        positions = [chrom_offsets[c] + d["pos"] for d in gwas if d["chrom"] == c]
        if positions:
            chrom_centers.append((c, np.mean(positions)))

    ax.set_xticks([center for _, center in chrom_centers])
    ax.set_xticklabels([str(c) for c, _ in chrom_centers], fontsize=8)

    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=12)
    ax.set_title("Variant Effect Sizes Across the Genome", fontsize=16, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.2)

    plt.tight_layout()
    path = os.path.join(outdir, "04_variant_effects.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {path}")


def plot_score_comparison(scores, outdir):
    """
    Plot 5: Per-Sample Score Comparison
    Horizontal lollipop chart ranking all patients by their PRS.
    Color-coded by risk category. Easy way to see who's at either extreme.
    """
    fig, ax = plt.subplots(figsize=(10, max(6, len(scores) * 0.35)))

    # sort by PRS
    sorted_scores = sorted(scores, key=lambda x: x["prs_zscore"])
    samples = [d["sample"].replace("PATIENT_", "P") for d in sorted_scores]
    zscores = [d["prs_zscore"] for d in sorted_scores]
    colors = [COLORS[d["risk_category"]] for d in sorted_scores]

    y_pos = range(len(samples))

    # lollipop chart — line from 0 to the score, dot at the end
    for i, (z, c) in enumerate(zip(zscores, colors)):
        ax.plot([0, z], [i, i], color=c, linewidth=2, alpha=0.7)
        ax.scatter(z, i, color=c, s=80, zorder=5, edgecolors="white", linewidth=0.5)

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(samples, fontsize=9)
    ax.set_xlabel("PRS Z-Score", fontsize=12)
    ax.set_title("Individual Risk Scores (Ranked)", fontsize=16, fontweight="bold")
    ax.axvline(x=0, color=COLORS["secondary"], linestyle="-", linewidth=0.8, alpha=0.5)
    ax.grid(axis="x", alpha=0.3)

    # legend
    legend_patches = [
        mpatches.Patch(color=COLORS["HIGH_RISK"], label="High Risk"),
        mpatches.Patch(color=COLORS["ELEVATED"], label="Elevated"),
        mpatches.Patch(color=COLORS["AVERAGE"], label="Average"),
        mpatches.Patch(color=COLORS["LOW_RISK"], label="Low Risk"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)

    plt.tight_layout()
    path = os.path.join(outdir, "05_score_comparison.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [OK] {path}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 60)
    print("Generating PRS Visualizations")
    print("=" * 60)

    scores = load_scores(args.scores)
    qc_samples = load_qc_samples(args.qc_stats)
    qc_summary = load_qc_summary(args.qc_stats)
    gwas = load_gwas(args.gwas)

    print(f"\nLoaded: {len(scores)} samples, {len(qc_samples)} QC records, {len(gwas)} variants")
    print(f"Output dir: {args.outdir}\n")

    plot_qc_summary(qc_samples, qc_summary, args.outdir)
    plot_prs_distribution(scores, args.outdir)
    plot_risk_stratification(scores, args.outdir)
    plot_variant_effects(gwas, args.outdir)
    plot_score_comparison(scores, args.outdir)

    print(f"\nAll 5 figures saved to {args.outdir}")


if __name__ == "__main__":
    main()
