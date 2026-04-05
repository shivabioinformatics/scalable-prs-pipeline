#!/usr/bin/env python3
"""
qc_filter.py
------------
Quality control is literally the first thing you do in any genomics pipeline.
If my input data is garbage, everything downstream is garbage too.

This script reads a VCF and checks:
  1. Per-variant call rate — if a variant is missing in too many samples, drop it
  2. Minor allele frequency (MAF) — super rare variants aren't useful for PRS
  3. Per-sample call rate — if a patient has too many missing calls, flag them
  4. Read depth — low depth means low confidence genotypes

I output a filtered VCF and a QC summary report (TSV) that the visualization
step will pick up later.
"""

import sys
import os
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="QC filter for VCF files")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--out-vcf", required=True, help="Output filtered VCF")
    parser.add_argument("--out-stats", required=True, help="Output QC stats TSV")
    parser.add_argument("--min-call-rate", type=float, default=0.95,
                        help="Minimum variant call rate (default: 0.95)")
    parser.add_argument("--min-maf", type=float, default=0.01,
                        help="Minimum minor allele frequency (default: 0.01)")
    parser.add_argument("--min-depth", type=float, default=10.0,
                        help="Minimum mean depth per variant (default: 10)")
    parser.add_argument("--min-sample-call-rate", type=float, default=0.90,
                        help="Minimum per-sample call rate (default: 0.90)")
    return parser.parse_args()


def parse_vcf_line(line):
    """pull apart a VCF data line into its components"""
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[1]
    rsid = fields[2]
    ref = fields[3]
    alt = fields[4]
    qual = fields[5]
    filt = fields[6]
    info = fields[7]
    fmt = fields[8]
    genotypes = fields[9:]
    return chrom, pos, rsid, ref, alt, qual, filt, info, fmt, genotypes


def compute_variant_stats(genotypes, fmt_fields):
    """
    For a single variant across all samples, figure out:
    - call rate: what fraction of samples actually have a genotype?
    - MAF: minor allele frequency
    - mean depth: average read depth across called samples

    If the call rate or MAF is too low, this variant gets tossed.
    """
    gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0
    dp_idx = fmt_fields.index("DP") if "DP" in fmt_fields else None

    n_samples = len(genotypes)
    n_called = 0
    n_alt = 0
    n_alleles = 0
    depths = []

    for geno_str in genotypes:
        geno_parts = geno_str.split(":")
        gt = geno_parts[gt_idx] if gt_idx < len(geno_parts) else "./."

        if gt in ("./.", ".|."):
            # missing call — this sample has no data for this variant
            continue

        n_called += 1
        alleles = gt.replace("|", "/").split("/")
        for a in alleles:
            if a.isdigit():
                n_alleles += 1
                if int(a) > 0:
                    n_alt += 1

        # grab read depth if available
        if dp_idx is not None and dp_idx < len(geno_parts):
            try:
                depths.append(int(geno_parts[dp_idx]))
            except ValueError:
                pass

    call_rate = n_called / n_samples if n_samples > 0 else 0.0
    af = n_alt / n_alleles if n_alleles > 0 else 0.0
    maf = min(af, 1.0 - af)  # MAF is the MINOR allele freq
    mean_depth = sum(depths) / len(depths) if depths else 0.0

    return call_rate, maf, mean_depth, n_called


def compute_sample_stats(sample_names, all_genotypes_by_sample):
    """
    Per-sample QC: how many variants does each sample have data for?
    If a sample has a really low call rate, something probably went wrong
    with their sequencing — library prep issues, contamination, whatever.
    """
    stats = {}
    for i, name in enumerate(sample_names):
        genos = all_genotypes_by_sample[i]
        total = len(genos)
        called = sum(1 for g in genos if g.split(":")[0] not in ("./.", ".|."))
        call_rate = called / total if total > 0 else 0.0
        depths = []
        for g in genos:
            parts = g.split(":")
            if len(parts) > 1:
                try:
                    depths.append(int(parts[1]))
                except (ValueError, IndexError):
                    pass
        mean_dp = sum(depths) / len(depths) if depths else 0.0
        stats[name] = {"call_rate": call_rate, "mean_depth": mean_dp, "variants_called": called, "variants_total": total}
    return stats


def main():
    args = parse_args()

    header_lines = []
    sample_names = []
    passing_lines = []
    all_genotypes_by_sample = defaultdict(list)

    # counters for the summary
    total_variants = 0
    dropped_call_rate = 0
    dropped_maf = 0
    dropped_depth = 0
    passing_variants = 0

    variant_stats_rows = []

    print(f"QC Filter — reading {args.vcf}")
    print(f"  Thresholds: call_rate >= {args.min_call_rate}, MAF >= {args.min_maf}, depth >= {args.min_depth}")

    with open(args.vcf, "r") as f:
        for line in f:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                header_lines.append(line)
                fields = line.strip().split("\t")
                sample_names = fields[9:]
                # initialize storage for per-sample tracking
                for i in range(len(sample_names)):
                    all_genotypes_by_sample[i] = []
                continue

            # data line — run QC checks
            total_variants += 1
            chrom, pos, rsid, ref, alt, qual, filt, info, fmt, genotypes = parse_vcf_line(line)
            fmt_fields = fmt.split(":")

            call_rate, maf, mean_depth, n_called = compute_variant_stats(genotypes, fmt_fields)

            # store genotypes for per-sample stats (even if variant gets dropped)
            for i, g in enumerate(genotypes):
                all_genotypes_by_sample[i].append(g)

            # apply filters — any failure and this variant is out
            if call_rate < args.min_call_rate:
                dropped_call_rate += 1
                continue
            if maf < args.min_maf:
                dropped_maf += 1
                continue
            if mean_depth < args.min_depth:
                dropped_depth += 1
                continue

            passing_variants += 1
            passing_lines.append(line)

            variant_stats_rows.append({
                "rsid": rsid, "chrom": chrom, "pos": pos,
                "call_rate": call_rate, "maf": maf, "mean_depth": mean_depth
            })

    # write filtered VCF
    with open(args.out_vcf, "w") as f:
        for h in header_lines:
            f.write(h)
        for line in passing_lines:
            f.write(line)

    # compute per-sample stats
    sample_stats = compute_sample_stats(sample_names, all_genotypes_by_sample)

    # write QC summary report
    with open(args.out_stats, "w") as f:
        # section 1: overall summary
        f.write("## VARIANT QC SUMMARY\n")
        f.write(f"total_variants\t{total_variants}\n")
        f.write(f"passing_variants\t{passing_variants}\n")
        f.write(f"dropped_low_call_rate\t{dropped_call_rate}\n")
        f.write(f"dropped_low_maf\t{dropped_maf}\n")
        f.write(f"dropped_low_depth\t{dropped_depth}\n")
        f.write(f"pass_rate\t{passing_variants/total_variants:.4f}\n")
        f.write("\n")

        # section 2: per-sample stats
        f.write("## PER-SAMPLE QC\n")
        f.write("sample\tcall_rate\tmean_depth\tvariants_called\tvariants_total\tpass_qc\n")
        for name in sample_names:
            s = sample_stats[name]
            passed = "PASS" if s["call_rate"] >= args.min_sample_call_rate else "FAIL"
            f.write(f"{name}\t{s['call_rate']:.4f}\t{s['mean_depth']:.1f}\t{s['variants_called']}\t{s['variants_total']}\t{passed}\n")

        f.write("\n")

        # section 3: per-variant stats (just the passing ones)
        f.write("## PER-VARIANT QC (passing only)\n")
        f.write("rsid\tchrom\tpos\tcall_rate\tmaf\tmean_depth\n")
        for row in variant_stats_rows:
            f.write(f"{row['rsid']}\t{row['chrom']}\t{row['pos']}\t{row['call_rate']:.4f}\t{row['maf']:.4f}\t{row['mean_depth']:.1f}\n")

    print(f"\nResults:")
    print(f"  {total_variants} total variants")
    print(f"  {passing_variants} passed QC ({passing_variants/total_variants*100:.1f}%)")
    print(f"  {dropped_call_rate} dropped (call rate)")
    print(f"  {dropped_maf} dropped (MAF)")
    print(f"  {dropped_depth} dropped (depth)")
    print(f"  Filtered VCF: {args.out_vcf}")
    print(f"  QC stats: {args.out_stats}")


if __name__ == "__main__":
    main()
